from __future__ import annotations
import re
import base64
import json
import os
import shutil
import subprocess
import tarfile
import time
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Literal
from fastapi import HTTPException
from pydantic import (
    BaseModel,
    Field,
    constr,
    conint,
    confloat,
    field_validator,
    model_validator,
)

from comps import CustomLogger, OpeaComponent, OpeaComponentRegistry
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")
# ---------------------------------------------------------------------------
# Globals / constants
# ---------------------------------------------------------------------------

logger = CustomLogger("relion_microservice")

BASE_WORKDIR = Path("/tmp/workspaces")
# Limit is ONLY for the input tar.gz from S3 (not for RELION outputs)
MAX_INPUT_TAR_BYTES = 50 * 1024 * 1024 * 1024  # 50 GB
RELION_EXE = "/opt/relion_5.0_cpu_benchmark_prefix/bin/relion_refine_mpi"

# ---------------------------------------------------------------------------
# Helper: encoding / metadata
# ---------------------------------------------------------------------------

def encode_base64_json(obj: Any) -> str:
    """
    Encode a JSON-serialisable Python object into a base64 string.
    Used to send metadata.json content inline to the client.
    """
    data = json.dumps(obj, default=str).encode("utf-8")
    return base64.b64encode(data).decode("ascii")


# ---------------------------------------------------------------------------
# Workspace management – Pattern C (always cleanup)
# ---------------------------------------------------------------------------

def create_workspace(user_id: str, workspace_id: Optional[str]) -> tuple[Path, str]:
    """
    Create a *unique* per-request workspace.

    - workspace_id from the client is treated as a logical label.
    - Physical directory always includes a random suffix to avoid collisions:
        /tmp/workspaces/<user_id>/<workspace_id>_<uuid>
      or
        /tmp/workspaces/<user_id>/<uuid>
      if workspace_id is None.

    Returns:
        (workdir_path, job_id_string)
    """
    user_dir = BASE_WORKDIR / user_id

    if workspace_id:
        job_id = f"{workspace_id}_{uuid.uuid4().hex}"
    else:
        job_id = uuid.uuid4().hex

    workdir = user_dir / job_id
    workdir.mkdir(parents=True, exist_ok=True)
    logger.info(f"[RELION] Workspace created: {workdir} (job_id={job_id})")
    return workdir, job_id
def make_tar_gz(src_dir: Path, tar_path: Path) -> None:
    """
    Create a .tar.gz archive from a directory.
    """
    if not src_dir.exists() or not src_dir.is_dir():
        raise HTTPException(
            status_code=500,
            detail=f"Output directory not found for tar: {src_dir}",
        )

    logger.info(f"[RELION] Creating tar.gz: {tar_path} from {src_dir}")

    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(src_dir, arcname=src_dir.name)

    logger.info(f"[RELION] Tar created: {tar_path} ({tar_path.stat().st_size / (1024**2):.2f} MB)")
def upload_file_to_s3(local_file: Path, s3_path: str) -> None:
    """
    Upload a single file to S3 using awscli.
    """
    if not local_file.exists():
        raise HTTPException(
            status_code=500,
            detail=f"Local file not found for upload: {local_file}",
        )

    if not s3_path.startswith("s3://"):
        raise HTTPException(
            status_code=400,
            detail="output_s3_prefix must start with 's3://'.",
        )

    cmd = [
        "aws",
        "s3",
        "cp",
        str(local_file),
        s3_path,
        "--no-progress",
    ]

    logger.info(f"[RELION] Uploading tar.gz to S3: {' '.join(cmd)}")

    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if proc.returncode != 0:
        logger.error(f"[RELION] aws s3 cp failed: {proc.stderr}")
        raise HTTPException(
            status_code=500,
            detail=f"Failed to upload tar.gz to S3: {proc.stderr}",
        )

    logger.info(f"[RELION] Uploaded tar.gz to {s3_path}")


def cleanup_workspace(workdir: Path) -> None:
    """
    Remove the workspace directory recursively (Pattern C).
    Called in a finally block to ensure cleanup on success or failure.
    """
    try:
        if workdir.exists():
            shutil.rmtree(workdir)
            logger.info(f"[RELION] Cleaned workspace: {workdir}")
    except Exception as e:
        logger.error(f"[RELION] Failed to cleanup workspace {workdir}: {e}")


# ---------------------------------------------------------------------------
# S3 download + extraction
# ---------------------------------------------------------------------------

def download_from_s3(s3_url: str, dest_path: Path) -> None:
    """
    Download a single file from S3 using awscli.

    Assumes correct IAM credentials/role are configured inside the container.
    """
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["aws", "s3", "cp", s3_url, str(dest_path), "--no-progress"]
    logger.info(f"[RELION] Downloading dataset from S3: {' '.join(cmd)}")

    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if proc.returncode != 0:
        logger.error(f"[RELION] aws s3 cp failed: {proc.stderr}")
        raise HTTPException(
            status_code=500,
            detail=f"Failed to download dataset from S3: {proc.stderr}",
        )

def extract_tar(tar_path: Path, dest_dir: Path) -> Path:
    """
    Extract a .tar.gz into dest_dir and return the detected project directory.

    Heuristic:
      - If exactly one subdirectory exists after extraction, use that as project dir.
      - Otherwise, use dest_dir itself.

    Only the input tar.gz is size-limited by MAX_INPUT_TAR_BYTES.
    RELION outputs can exceed this size (no limit on outputs).
    """
    size_bytes = tar_path.stat().st_size
    if size_bytes > MAX_INPUT_TAR_BYTES:
        raise HTTPException(
            status_code=413,
            detail=(
                f"Dataset too large (> {MAX_INPUT_TAR_BYTES // (1024**3)} GB). "
                "Increase MAX_INPUT_TAR_BYTES if this is expected."
            ),
        )

    size_gb = size_bytes / (1024**3)
    logger.info(f"[RELION] Extracting {tar_path} ({size_gb:.2f} GB) into {dest_dir}")
    start = time.time()
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(dest_dir)
    logger.info(f"[RELION] Extraction finished in {time.time() - start:.2f} s")

    # Remove tar to save space
    tar_path.unlink(missing_ok=True)

    # Detect project dir
    subdirs = [p for p in dest_dir.iterdir() if p.is_dir()]
    if len(subdirs) == 1:
        project_dir = subdirs[0]
    else:
        project_dir = dest_dir

    logger.info(f"[RELION] Project dir detected: {project_dir}")
    return project_dir


# ---------------------------------------------------------------------------
# Path helpers (prevent path traversal)
# ---------------------------------------------------------------------------

def _ensure_relative_safe(rel_path: str, field_name: str) -> str:
    """
    Ensure a user-provided path is relative and does not contain '..'.
    """
    if rel_path.startswith("/"):
        raise HTTPException(
            status_code=400,
            detail=f"{field_name} must be a relative path, not absolute.",
        )
    parts = [p for p in rel_path.split("/") if p not in ("", ".")]
    if any(p == ".." for p in parts):
        raise HTTPException(
            status_code=400,
            detail=f"{field_name} must not contain '..' path segments.",
        )
    return "/".join(parts)


def resolve_in_project(project_dir: Path, rel_path: str, field_name: str, must_exist: bool) -> Path:
    """
    Resolve a relative path inside the extracted project directory, ensuring
    it does not escape the project root. Optionally enforce existence.
    """
    rel_clean = _ensure_relative_safe(rel_path, field_name)
    full = (project_dir / rel_clean).resolve()
    project_root = project_dir.resolve()

    try:
        is_inside = full.is_relative_to(project_root)  # type: ignore[attr-defined]
    except AttributeError:
        is_inside = str(full).startswith(str(project_root))

    if not is_inside:
        raise HTTPException(
            status_code=400,
            detail=f"{field_name} resolves outside of project directory.",
        )

    if must_exist and not full.exists():
        raise HTTPException(
            status_code=400,
            detail=f"{field_name} path not found inside dataset: {rel_path}",
        )

    return full


# ---------------------------------------------------------------------------
# Pydantic models: Flags + Input + Output
# ---------------------------------------------------------------------------

class RelionRefineFlags(BaseModel):
    """
    Full set of relion_refine_mpi flags relevant for:
      - 2D classification
      - 3D classification
      - 3D auto-refine

    This microservice is command-mode:
      - Client supplies *all* flags in a structured way.
      - Server validates semantics + dataset paths.
    """

    # Paths (relative to project root)
    i: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description="Input particles star file relative to project root, e.g. 'Particles/shiny_2sets.star'. (Required)",
    )
    o: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description="Output rootname relative to project root, e.g. '2D/2D', '3D/3D', or 'Refine3D/job019/run'. (Required)",
    )
    ref: Optional[constr(strip_whitespace=True, min_length=1)] = Field(
        None,
        description=(
            "Optional reference map path relative to project root. "
            "May include a RELION suffix, e.g. 'emd_2660.map:mrc' or "
            "'Class3D/job016/run_it025_class004_box256.mrc'."
        ),
    )

    # Classification / refinement semantics
    K: Optional[conint(ge=1)] = Field(
        None,
        description="--K Number of classes/references. Required for 2D/3D classification, not for pure auto-refine.",
    )
    auto_refine: bool = Field(
        False,
        description="--auto_refine for 3D auto-refinement mode.",
    )

    # Boolean flags
    dont_combine_weights_via_disc: bool = Field(
        False,
        description="--dont_combine_weights_via_disc",
    )
    ctf: bool = Field(
        False,
        description="--ctf (perform CTF correction).",
    )
    zero_mask: bool = Field(
        False,
        description="--zero_mask",
    )
    norm: bool = Field(
        False,
        description="--norm",
    )
    scale: bool = Field(
        False,
        description="--scale",
    )
    firstiter_cc: bool = Field(
        False,
        description="--firstiter_cc",
    )
    flatten_solvent: bool = Field(
        False,
        description="--flatten_solvent",
    )
    split_random_halves: bool = Field(
        False,
        description="--split_random_halves.",
    )
    cpu: bool = Field(
        True,
        description="--cpu (Intel CPU vectorisation). Recommended True inside this microservice.",
    )
    @field_validator("cpu")
    @classmethod
    def _force_cpu(cls, v: bool) -> bool:
        """
        This microservice is CPU-only.
        flags.cpu must always be True; GPU mode is not supported.
        """
        if not v:
            raise ValueError(
                "This RELION microservice is CPU-only; flags.cpu must be True."
            )
        return True

    preread_images: bool = Field(
        False,
        description="--preread_images (pre-read particle images into memory).",
    )
    center_classes: bool = Field(
        False,
        description="--center_classes (re-center 2D classes).",
    )
    blush: bool = Field(
        False,
        description="--blush (3D classification option used in some pipelines).",
    )
    auto_ignore_angles: bool = Field(
        False,
        description="--auto_ignore_angles (auto-refine control flag).",
    )
    auto_resol_angles: bool = Field(
        False,
        description="--auto_resol_angles (auto-refine control flag).",
    )

    # Numeric flags
    tau2_fudge: Optional[confloat(gt=0)] = Field(
        None,
        description="--tau2_fudge regularisation parameter (>0).",
    )
    particle_diameter: Optional[confloat(gt=0)] = Field(
        None,
        description="--particle_diameter in Angstroms (>0).",
    )
    oversampling: Optional[conint(ge=0)] = Field(
        None,
        description="--oversampling (0,1,2,...).",
    )
    psi_step: Optional[confloat(gt=0)] = Field(
        None,
        description="--psi_step in degrees (>0).",
    )
    offset_range: Optional[confloat(gt=0)] = Field(
        None,
        description="--offset_range in pixels (>0).",
    )
    offset_step: Optional[confloat(gt=0)] = Field(
        None,
        description="--offset_step in pixels (>0).",
    )
    random_seed: Optional[int] = Field(
        None,
        description="--random_seed (int).",
    )
    pad: Optional[conint(ge=1)] = Field(
        None,
        description="--pad (Must be >=1).",
    )
    pool: conint(ge=1) = Field(
        ...,
        description="--pool. Must be between 1 and 2*j.",
    )
    j: conint(ge=1) = Field(
        ...,
        description="--j (threads per MPI rank). Must be >= 1 and <= total_cpus/num_mpi.",
    )
    iter: Optional[conint(ge=1)] = Field(
        None,
        description="--iter (>=1). For auto-refine, RELION controls iterations internally.",
    )
    ini_high: Optional[confloat(gt=0)] = Field(
        None,
        description="--ini_high (Angstrom, >0).",
    )
    healpix_order: Optional[conint(ge=0)] = Field(
        None,
        description="--healpix_order (>=0).",
    )
    auto_local_healpix_order: Optional[conint(ge=0)] = Field(
        None,
        description="--auto_local_healpix_order (>=0).",
    )
    low_resol_join_halves: Optional[confloat(gt=0)] = Field(
        None,
        description="--low_resol_join_halves (Angstrom, >0).",
    )
    sym: Optional[constr(strip_whitespace=True, min_length=1)] = Field(
        None,
        description="--sym symmetry group, e.g. 'C1', 'D2'.",
    )

    # Extra CLI args (rare / escape hatch)
    extra_args: Optional[List[str]] = Field(
        None,
        description="Optional extra CLI args to append verbatim to relion_refine_mpi.",
    )

class RelionInput(BaseModel):
    """
    Request schema for RELION classification + auto-refine microservice.

    All relion_refine_mpi flags are part of `flags` and come from the client.

    - `mode` is used ONLY for semantic validation (2d / 3d / autorefine).
    - Server:
        * downloads & extracts dataset from S3,
        * validates paths inside the dataset,
        * validates flags according to mode,
        * constructs relion_refine_mpi command,
        * runs RELION,
        * returns metadata.json as base64.
    """

    # ---- Identity / workspace isolation ----
    user_id: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description="Logical user id used as top-level directory under /tmp/workspaces.",
    )

    workspace_id: Optional[constr(strip_whitespace=True, min_length=1)] = Field(
        None,
        description=(
            "Optional workspace/job id. "
            "If not provided, the server will generate an ID."
        ),
    )

    # ---- dataset location in S3 ----
    dataset_s3_url: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description="S3 URL of a .tar.gz RELION dataset, e.g. s3://bucket/relion_benchmark.tar.gz",
    )
    output_s3_prefix: Optional[constr(strip_whitespace=True, min_length=1)] = Field(
            None,
            description=(
                "Optional S3 prefix where the parent folder of flags.o will be uploaded.\n"
                "Example: s3://my-bucket/relion_outputs/job_2d_001\n"
                "If flags.o = '2D/2D', then the local '2D' folder will be copied to "
                "s3://.../job_2d_001/2D"
                ),
            )
    @field_validator("output_s3_prefix")
    @classmethod
    def _check_output_s3_prefix(cls, v: Optional[str]) -> Optional[str]:
        if v is not None and not v.startswith("s3://"):
            raise ValueError("output_s3_prefix must start with 's3://'.")
        return v


    # ---- high-level mode (for semantics only) ----
    mode: Literal["2d", "3d", "autorefine"] = Field(
        "2d",
        description="RELION mode for semantic validation: '2d', '3d', or 'autorefine'.",
    )

    # ---- MPI ----
    num_mpi: int = Field(
        2,
        description="mpirun -n value (number of MPI ranks, >= 1).",
    )

    # ---- relion_refine_mpi flags ----
    flags: RelionRefineFlags = Field(
        ...,
        description="Structured relion_refine_mpi flags.",
    )
    @field_validator("user_id")
    @classmethod
    def _check_user_id(cls, v: str) -> str:
        if not SAFE_ID_RE.match(v):
           raise ValueError(
             "user_id may only contain letters, digits, '.', '_' or '-'."
              )
        return v
    @field_validator("workspace_id")
    @classmethod
    def _check_workspace_id(cls, v: Optional[str]) -> Optional[str]:
        if v is not None and not SAFE_ID_RE.match(v):
           raise ValueError(
               "workspace_id may only contain letters, digits, '.', '_' or '-'."
            )
        return v

    @field_validator("dataset_s3_url")
    @classmethod
    def _check_s3_url(cls, v: str) -> str:
        if not v.startswith("s3://"):
            raise ValueError("dataset_s3_url must start with 's3://'.")
        if not v.endswith(".tar.gz"):
            raise ValueError("dataset_s3_url must point to a .tar.gz file.")
        return v

    @field_validator("mode")
    @classmethod
    def _normalize_mode(cls, v: str) -> str:
        return v.lower()

    @model_validator(mode="after")
    def _check_mode_semantics(self) -> "RelionInput":
        """
        Mode-level semantic validation (dataset-agnostic):

        Common:
          - flags.i and flags.o must be set.

        mode == '2d':
          - flags.ref MUST be None.
          - flags.K and flags.tau2_fudge MUST be provided.

        mode == '3d':
          - flags.ref MUST be provided.
          - flags.K, flags.tau2_fudge, flags.healpix_order,
            flags.ini_high, flags.sym MUST be provided.

        mode == 'autorefine':
          - flags.ref MUST be provided.
          - flags.auto_refine MUST be True.
          - flags.split_random_halves MUST be True.
          - healpix_order, ini_high, auto_local_healpix_order,
            low_resol_join_halves, sym MUST be provided.
        """
        mode = self.mode.lower()
        f = self.flags

        # Common:
        if not f.i:
            raise HTTPException(status_code=400, detail="flags.i (input star file) is required.")
        if not f.o:
            raise HTTPException(status_code=400, detail="flags.o (output root) is required.")

        if mode == "2d":
            if f.ref is not None:
                raise HTTPException(
                status_code=400,
                detail="For 2D classification, flags.ref must NOT be set.",
            )
            if f.K is None:
                raise HTTPException(status_code=400, detail="For 2D, flags.K is required.")
            if f.tau2_fudge is None:
                raise HTTPException(status_code=400, detail="For 2D, flags.tau2_fudge is required.")

        elif mode == "3d":
            if not f.ref:
                raise HTTPException(status_code=400, detail="For 3D, flags.ref is required.")
            if f.K is None:
                raise HTTPException(status_code=400, detail="For 3D, flags.K is required.")
            if f.tau2_fudge is None:
                raise HTTPException(status_code=400, detail="For 3D, flags.tau2_fudge is required.")
            if f.healpix_order is None:
                raise HTTPException(status_code=400, detail="For 3D, flags.healpix_order is required.")
            if f.ini_high is None:
                raise HTTPException(status_code=400, detail="For 3D, flags.ini_high is required.")
            if f.sym is None:
                raise HTTPException(status_code=400, detail="For 3D, flags.sym is required.")

        elif mode == "autorefine":
            if not f.ref:
                raise HTTPException(status_code=400, detail="For autorefine, flags.ref is required.")
            if not f.auto_refine:
                raise HTTPException(status_code=400, detail="For autorefine, flags.auto_refine must be True.")
            if not f.split_random_halves:
                raise HTTPException(
                status_code=400,
                detail="For autorefine, flags.split_random_halves must be True.",
            )
            if f.iter is not None:
                raise HTTPException(
                        status_code=400,
                        detail="In auto-refine mode, --iter must NOT be provided. RELION controls iterations internally."
                        )
            if f.healpix_order is None:
                raise HTTPException(status_code=400, detail="For autorefine, flags.healpix_order is required.")
            if f.ini_high is None:
                raise HTTPException(status_code=400, detail="For autorefine, flags.ini_high is required.")
            if f.auto_local_healpix_order is None:
                raise HTTPException(
                status_code=400,
                detail="For autorefine, flags.auto_local_healpix_order is required.",
            )
            if f.low_resol_join_halves is None:
                raise HTTPException(
                status_code=400,
                detail="For autorefine, flags.low_resol_join_halves is required.",
            )
            if f.sym is None:
                raise HTTPException(status_code=400, detail="For autorefine, flags.sym is required.")

        else:
            raise HTTPException(status_code=400, detail="mode must be '2d', '3d', or 'autorefine'.")
        if f.pad is not None and f.pad <= 0:
            raise HTTPException(
                    status_code=400,
                    detail=f"flags.pad must be >=1 (got {f.pad}).",
                    )
        total_cpus = os.cpu_count() or 1
        if self.num_mpi <= 1:
            raise HTTPException(
                    status_code=400,
                    detail=f"num_mpi must be > 1 (got {self.num_mpi}).",
                    )
        if f.j is None or f.j < 1:
            raise HTTPException(
                    status_code=400,
                    detail=f"flags.j must be >= 1 (got {f.j}).",
                    )
        max_j = total_cpus // self.num_mpi
        if max_j < 1:
            raise HTTPException(
                    status_code=400,
                    detail=(
                        f"num_mpi={self.num_mpi} is too large for this machine "
                        f"(total_cpus={total_cpus})."
                        ),
                    )
        if f.j > max_j:
            raise HTTPException(
                    status_code=400,
                    detail=(
                        f"Invalid j: j={f.j}, but total_cpus/num_mpi={max_j}. "
                        "Constraint: j <= total_cpus / num_mpi."
                        ),
                    )
        if f.pool is None or f.pool < 1:
            raise HTTPException(
                    status_code=400,
                    detail=f"flags.pool must be >= 1 (got {f.pool}).",
                    )
        max_pool = 2 * f.j
        if f.pool > max_pool:
            raise HTTPException(
                    status_code=400,
                    detail=(
                        f"Invalid pool: pool={f.pool}, j={f.j}, max allowed=2*j={max_pool}."
                        ),
                    )
        return self


class RelionOutput(BaseModel):
    """
    Response schema for RELION microservice.
    Only metadata is returned (no project tarball).
    """

    status: str
    message: Optional[str] = None
    metadata_b64: Optional[str] = None  # base64-encoded metadata.json contents
    job_id: Optional[str] = None        # physical job id / directory name


# ---------------------------------------------------------------------------
# Command builder for relion_refine_mpi
# ---------------------------------------------------------------------------

def build_relion_cmd(
    project_dir: Path,
    num_mpi: int,
    flags: RelionRefineFlags,
) -> Dict[str, Any]:
    """
    Build the relion_refine_mpi command line from validated flags.

    - Resolves relative paths for --i, --o, --ref inside project_dir.
    - Ensures paths stay within project_dir and exist where required.
    - For flags.ref, supports both:
        * 'emd_2660.map:mrc'  (benchmark style)
        * 'Class3D/.../run_it025_class004_box256.mrc' (plain .mrc)
    """
    if num_mpi <= 1:
        raise HTTPException(
            status_code=400,
            detail=f"num_mpi must be > 1 (got {num_mpi})",
        )

    # Resolve input and output paths
    i_path = resolve_in_project(
        project_dir, flags.i, field_name="flags.i (input star)", must_exist=True
    )
    o_path = resolve_in_project(
        project_dir, flags.o, field_name="flags.o (output root)", must_exist=False
    )

    # ---- REF HANDLING (GENERIC) ----
    ref_path: Optional[Path] = None
    ref_suffix: Optional[str] = None

    if flags.ref:
        ref_str = flags.ref
        if ":" in ref_str:
            base_rel, ref_suffix = ref_str.split(":", 1)
        else:
            base_rel, ref_suffix = ref_str, None

        ref_path = resolve_in_project(
            project_dir,
            base_rel,
            field_name="flags.ref (reference map)",
            must_exist=True,
        )

    # Ensure parent dir of output root exists
    o_parent = o_path.parent
    o_parent.mkdir(parents=True, exist_ok=True)

    cmd: List[str] = ["mpirun", "-n", str(num_mpi), RELION_EXE]

    # Required paths
    cmd.extend(["--i", str(i_path)])
    cmd.extend(["--o", str(o_path)])

    if ref_path is not None:
        if ref_suffix:
            cmd.extend(["--ref", f"{str(ref_path)}:{ref_suffix}"])
        else:
            cmd.extend(["--ref", str(ref_path)])

    # Boolean flags
    if flags.dont_combine_weights_via_disc:
        cmd.append("--dont_combine_weights_via_disc")
    if flags.ctf:
        cmd.append("--ctf")
    if flags.zero_mask:
        cmd.append("--zero_mask")
    if flags.norm:
        cmd.append("--norm")
    if flags.scale:
        cmd.append("--scale")
    if flags.firstiter_cc:
        cmd.append("--firstiter_cc")
    if flags.flatten_solvent:
        cmd.append("--flatten_solvent")
    if flags.split_random_halves:
        cmd.append("--split_random_halves")
    if flags.cpu:
        cmd.append("--cpu")
    if flags.auto_refine:
        cmd.append("--auto_refine")
    if flags.preread_images:
        cmd.append("--preread_images")
    if flags.center_classes:
        cmd.append("--center_classes")
    if flags.blush:
        cmd.append("--blush")
    if flags.auto_ignore_angles:
        cmd.append("--auto_ignore_angles")
    if flags.auto_resol_angles:
        cmd.append("--auto_resol_angles")

    # Numeric flags
    def add_num(flag: str, value: Any):
        cmd.extend([flag, str(value)])

    if flags.tau2_fudge is not None:
        add_num("--tau2_fudge", flags.tau2_fudge)
    if flags.particle_diameter is not None:
        add_num("--particle_diameter", flags.particle_diameter)
    if flags.K is not None:
        add_num("--K", flags.K)
    if flags.oversampling is not None:
        add_num("--oversampling", flags.oversampling)
    if flags.psi_step is not None:
        add_num("--psi_step", flags.psi_step)
    if flags.offset_range is not None:
        add_num("--offset_range", flags.offset_range)
    if flags.offset_step is not None:
        add_num("--offset_step", flags.offset_step)
    if flags.random_seed is not None:
        add_num("--random_seed", flags.random_seed)
    if flags.pad is not None:
        add_num("--pad", flags.pad)
    add_num("--pool", flags.pool)
    add_num("--j", flags.j)
    if flags.iter is not None:
        add_num("--iter", flags.iter)
    if flags.ini_high is not None:
        add_num("--ini_high", flags.ini_high)
    if flags.healpix_order is not None:
        add_num("--healpix_order", flags.healpix_order)
    if flags.auto_local_healpix_order is not None:
        add_num("--auto_local_healpix_order", flags.auto_local_healpix_order)
    if flags.low_resol_join_halves is not None:
        add_num("--low_resol_join_halves", flags.low_resol_join_halves)
    if flags.sym is not None:
        add_num("--sym", flags.sym)

    # Extra args (verbatim)
    if flags.extra_args:
        cmd.extend(flags.extra_args)

    return {
        "cmd": cmd,
        "workdir": project_dir,
        "output_dir": o_parent,
    }


# ---------------------------------------------------------------------------
# OPEA Component – Pattern C cleanup
# ---------------------------------------------------------------------------

@OpeaComponentRegistry.register("OPEA_OMICS_RELION")
class OpeaRelion(OpeaComponent):
    """
    OPEA component wrapping RELION 5.0 classification + auto-refine pipelines.

    Design:
      - Command-mode with full relion_refine_mpi flags from the client.
      - Input dataset is a .tar.gz on S3.
      - Each request gets its own workspace under /tmp/workspaces/user_id/job_id.
      - Workspace is ALWAYS deleted in `finally`, irrespective of success/failure.
      - Client receives metadata.json contents as base64 in the response.
      - No timeout is enforced at microservice level.
    """

    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)
        self.config = config or {}

    async def invoke(self, input: RelionInput) -> RelionOutput:
        start = time.time()
        workdir: Optional[Path] = None
        job_id: Optional[str] = None
        output_tar_s3_path: Optional[str] = None

        try:
            # 1) Per-request workspace
            workdir, job_id = create_workspace(input.user_id, input.workspace_id)

            # 2) Download dataset from S3 -> dataset.tar.gz
            tar_path = workdir / "dataset.tar.gz"
            download_from_s3(input.dataset_s3_url, tar_path)

            # 3) Extract tar.gz -> project dir
            project_dir = extract_tar(tar_path, workdir)

            # 4) Build RELION command from flags
            cmd_info = build_relion_cmd(
                project_dir=project_dir,
                num_mpi=input.num_mpi,
                flags=input.flags,
            )
            cmd = cmd_info["cmd"]
            output_dir: Path = cmd_info["output_dir"]
            logger.info(f"[RELION] Running command: {' '.join(cmd)}")

            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(input.flags.j)

            # 5) Run RELION (no timeout here)
            proc = subprocess.run(
                cmd,
                cwd=str(cmd_info["workdir"]),
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            elapsed = time.time() - start
            logger.info(
                f"[RELION] Command finished in {elapsed:.2f} seconds, rc={proc.returncode}"
            )

            status_str = "success" if proc.returncode == 0 else "error"
            msg = (
                "RELION run completed successfully."
                if proc.returncode == 0
                else "RELION run failed. See metadata for details."
            )
            if input.output_s3_prefix is not None:
                try:
                    tar_path = output_dir.parent / f"{output_dir.name}.tar.gz"
                    make_tar_gz(output_dir, tar_path)
                    dest_s3 = input.output_s3_prefix.rstrip("/") + f"/{tar_path.name}"
                    upload_file_to_s3(tar_path, dest_s3)
                    output_tar_s3_path = dest_s3
                except HTTPException:
                    # re-raise so the client sees a 500 if upload fails
                    raise

            metadata: Dict[str, Any] = {
                "status": status_str,
                "message": msg,
                "mode": input.mode,
                "elapsed_sec": elapsed,
                "returncode": proc.returncode,
                "user_id": input.user_id,
                "workspace_id": input.workspace_id,
                "job_id": job_id,
                "dataset_s3_url": input.dataset_s3_url,
                "num_mpi": input.num_mpi,
                "flags": input.flags.model_dump(),
                "command": cmd,
                "stdout_tail": proc.stdout,
                "stderr_tail": proc.stderr,
                "output_tar_s3_path": output_tar_s3_path,
            }

            # Write metadata.json for conceptual completeness (will be deleted with workspace)
            meta_path = workdir / "metadata.json"
            meta_path.write_text(json.dumps(metadata, indent=2))

            metadata_b64 = encode_base64_json(metadata)

            return RelionOutput(
                status=status_str,
                message=msg,
                metadata_b64=metadata_b64,
                job_id=job_id,
            )

        except HTTPException:
            # HTTPException is passed through unchanged
            raise

        except Exception as e:
            logger.error(f"[RELION] Unexpected error in RELION microservice: {e}")
            # We still return a 500 with a minimal metadata-like payload encoded
            fallback_meta = {
                "status": "error",
                "message": str(e),
                "user_id": getattr(input, "user_id", None),
                "workspace_id": getattr(input, "workspace_id", None),
                "job_id": job_id,
            }
            metadata_b64 = encode_base64_json(fallback_meta)
            raise HTTPException(
                status_code=500,
                detail=metadata_b64,
            )

        finally:
            # Pattern C: ALWAYS clean workspace, regardless of success or failure
            if workdir is not None:
                cleanup_workspace(workdir)

    def check_health(self) -> bool:
        """
        Simple health check: verify that the RELION executable exists.
        """
        return Path(RELION_EXE).exists()

