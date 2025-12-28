from io import BytesIO
from typing import Optional, Tuple
from pydantic import BaseModel, Field, constr, conint, confloat
import base64
import io
import tarfile
import uuid
from pathlib import Path
import re
import os
import shutil
import subprocess
import time
from fastapi import HTTPException
import json
from comps import register_statistics
from comps import (
    CustomLogger,
    OpeaComponent,
    OpeaComponentRegistry,
    statistics_dict,
)

logger = CustomLogger("autodock_vina_microservice")

MAX_RESULT_TAR_BYTES = 500 * 1024 * 1024  # 500 MB
MAX_TAR_BYTES = 500 * 1024 * 1024


def create_output_tar_b64_with_metadata(
    output_file: Path,
    metadata: dict,
    max_bytes: int = MAX_RESULT_TAR_BYTES,
) -> tuple[str, int]:
    """
    Create a tar.gz archive containing ONLY:
      - the given output_file
      - a metadata.json file with job info

    Returns:
      - base64-encoded tar.gz string
      - raw tar.gz size in bytes

    Enforces a max_bytes limit on the raw tar size.
    """
    if not output_file.is_file():
        raise ValueError(f"Output file does not exist: {output_file}")

    buf = BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        # 1) Add the output file with just its filename
        tar.add(output_file, arcname=output_file.name)

        # 2) Add metadata.json from in-memory bytes
        metadata_bytes = json.dumps(metadata, indent=2).encode("utf-8")
        info = tarfile.TarInfo(name="metadata.json")
        info.size = len(metadata_bytes)
        info.mtime = int(time.time())
        tar.addfile(info, BytesIO(metadata_bytes))

    buf.seek(0)
    tar_bytes = buf.read()

    size = len(tar_bytes)
    if size > max_bytes:
        raise ValueError(
            f"Result tar.gz is too large ({size} bytes, limit {max_bytes} bytes). "
            "Consider disabling return_tarball or reducing workspace size."
        )

    b64 = base64.b64encode(tar_bytes).decode("ascii")
    return b64, size


class VinaInput(BaseModel):
    """
    Input schema for the AutoDock Vina microservice.

    Users send JSON with:
    - identity / workspace info
    - docking dataset as base64 tar.gz
    - Vina flags as individual JSON fields

    NOTE: This microservice is restricted to a single receptor + single ligand.
    No flex/batch/dir/write_maps or advanced tuning flags are exposed.
    """

    # ---------- Multi-user / workspace ----------
    user_id: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description="Logical user identifier (no spaces). Used for workspace isolation.",
    )

    workspace_id: Optional[constr(strip_whitespace=True, min_length=1)] = Field(
        None,
        description=(
            "Optional workspace/job group identifier. "
            "If not provided, the server will generate one."
        ),
    )

    # ---------- Dataset from client ----------
    dataset_tar_b64: str = Field(
        ...,
        description=(
            "Base64-encoded tar.gz of a single docking dataset directory "
            "(e.g. containing protein.pdbqt, ligands, etc.)."
        ),
    )

    # ---------- Required Vina input files (inside the extracted dataset) ----------
    receptor: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description="Relative path/name of the receptor PDBQT inside the dataset, e.g. 'protein.pdbqt'.",
    )

    ligand: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description="Relative path/name of the ligand PDBQT inside the dataset, e.g. 'rand-1.pdbqt'.",
    )

    # ---------- Scoring function ----------
    scoring: Optional[constr(strip_whitespace=True)] = Field(
        "vina",
        description=(
            "Scoring function. For this microservice, it is restricted to 'vina'. "
            "Any other value will be rejected in validation."
        ),
    )

    # ---------- Search space (required) ----------
    center_x: float = Field(..., description="X coordinate of the search box center (Angstrom).")
    center_y: float = Field(..., description="Y coordinate of the search box center (Angstrom).")
    center_z: float = Field(..., description="Z coordinate of the search box center (Angstrom).")

    size_x: confloat(gt=0) = Field(
        ...,
        description="Size of the search box in X dimension (Angstrom, must be > 0).",
    )
    size_y: confloat(gt=0) = Field(
        ...,
        description="Size of the search box in Y dimension (Angstrom, must be > 0).",
    )
    size_z: confloat(gt=0) = Field(
        ...,
        description="Size of the search box in Z dimension (Angstrom, must be > 0).",
    )

    # ---------- Output options (optional) ----------
    out: Optional[constr(strip_whitespace=True, min_length=1)] = Field(
        None,
        description=(
            "Optional output PDBQT filename for docked poses, relative to docking directory. "
            "If not provided, Vina will choose a default based on the ligand name."
        ),
    )

    # ---------- Misc options (optional) ----------
    cpu: Optional[conint(ge=1, le=800)] = Field(
        None,
        description="Number of CPUs to use. If not set, Vina will auto-detect.",
    )

    seed: Optional[int] = Field(
        None,
        description="Explicit random seed for Vina. If not set, Vina chooses its own.",
    )

    exhaustiveness: conint(ge=1, le=1024) = Field(
        8,
        description="Exhaustiveness of the global search. Roughly proportional to time. Default: 8.",
    )

    # ---------- Tarball return option ----------
    return_tarball: bool = Field(
        False,
        description=(
            "If true, the service will return a base64-encoded tar.gz containing "
            "the Vina output PDBQT file and a metadata.json with job details."
        ),
    )


class VinaOutput(BaseModel):
    """
    Output schema for the AutoDock Vina microservice.
    """

    # High-level status
    status: str = Field(
        ...,
        description="Status of the docking request, e.g. 'success' or 'error'.",
    )
    message: Optional[str] = Field(
        None,
        description="Human-readable message summarizing the result.",
    )

    # Identity / bookkeeping
    user_id: Optional[str] = Field(
        None,
        description="User identifier echoed back from the request.",
    )
    workspace_id: Optional[str] = Field(
        None,
        description="Workspace identifier used for this job.",
    )
    job_id: Optional[str] = Field(
        None,
        description="Server-generated unique job identifier for this docking run.",
    )

    # Output files / docking result (path is only for logging/debugging)
    output_pdbqt: Optional[str] = Field(
        None,
        description=(
            "Path to the output PDBQT file with docked poses, "
            "relative to the service workspace root (may not exist anymore if cleanup is enabled)."
        ),
    )

    # Performance / debugging info
    runtime_seconds: Optional[float] = Field(
        None,
        description="Total wall-clock time spent running the docking job (seconds).",
    )
    stdout: Optional[str] = Field(
        None,
        description="Raw stdout captured from the Vina process (for debugging/logging).",
    )
    stderr: Optional[str] = Field(
        None,
        description="Raw stderr captured from the Vina process (for debugging/logging).",
    )

    # Tarball of result
    result_tar_b64: Optional[str] = Field(
        None,
        description=(
            "Base64-encoded tar.gz containing only the Vina output PDBQT file "
            "and a metadata.json with job details (present only if return_tarball=true "
            "and within size limits)."
        ),
    )
    result_tar_size_bytes: Optional[int] = Field(
        None,
        description="Size in bytes of the raw tar.gz before base64 encoding.",
    )


def sanitize_name(name: str | None) -> str | None:
    """
    Sanitize user-provided identifiers (user_id, workspace_id) to a safe form.

    - Strips whitespace
    - Replaces anything not [a-zA-Z0-9_-] with '_'
    - Returns None if name becomes empty
    """
    if not name:
        return None
    clean = re.sub(r"[^a-zA-Z0-9_-]", "_", name.strip())
    return clean or None


def build_workspace_dirs(
    base_root: Path,
    user_id_raw: str,
    workspace_id_raw: str | None,
) -> Tuple[Path, str, str]:
    """
    Build and create /<base_root>/<user_id>/<workspace_id>/<job_id>.

    Each request gets its own unique job_id directory.

    Returns:
        workspace_dir: Path to the job-specific working directory
        workspace_id:  Sanitized workspace id used
        job_id:        Newly generated job id (UUID hex)
    """
    user_id = sanitize_name(user_id_raw)
    if not user_id:
        raise ValueError("Invalid user_id after sanitization")

    workspace_id = sanitize_name(workspace_id_raw) if workspace_id_raw else None
    if not workspace_id:
        workspace_id = uuid.uuid4().hex

    job_id = uuid.uuid4().hex

    workspace_dir = base_root / user_id / workspace_id / job_id
    workspace_dir.mkdir(parents=True, exist_ok=True)

    return workspace_dir, workspace_id, job_id


def decode_tar_b64(dataset_tar_b64: str, max_bytes: int = MAX_TAR_BYTES) -> bytes:
    """
    Decode base64-encoded tar.gz data into raw bytes, enforcing a size limit.

    Raises:
        ValueError if base64 is invalid or size exceeds limit.
    """
    try:
        raw = base64.b64decode(dataset_tar_b64)
    except Exception as e:
        raise ValueError(f"Invalid base64 data: {e}") from e

    if len(raw) > max_bytes:
        raise ValueError(
            f"Decoded tar.gz is too large ({len(raw)} bytes, limit {max_bytes} bytes)"
        )

    return raw


def _is_safe_tar_member(member: tarfile.TarInfo) -> bool:
    """
    Check that a tar member's path is safe:
    - not absolute
    - does not contain '..' components
    """
    name = member.name
    # Absolute path
    if name.startswith("/"):
        return False

    # Any .. component in the path
    parts = Path(name).parts
    if ".." in parts:
        return False

    return True


def safe_extract_tar_bytes(tar_bytes: bytes, dest_dir: Path) -> None:
    """
    Safely extract a tar.gz from raw bytes into dest_dir.

    - Uses size-checked bytes (caller must ensure limit)
    - Prevents absolute and parent ('..') paths

    Raises:
        ValueError if tar is invalid or contains unsafe paths.
    """
    file_like = io.BytesIO(tar_bytes)
    try:
        with tarfile.open(fileobj=file_like, mode="r:gz") as tar:
            for member in tar.getmembers():
                if not _is_safe_tar_member(member):
                    raise ValueError(f"Unsafe path in tar member: {member.name}")
            tar.extractall(path=dest_dir)
    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Failed to extract tar.gz: {e}") from e


def find_docking_root(workspace_dir: Path) -> Path:
    """
    Heuristic to find the docking root directory after extracting the tar.

    If there's exactly one subdirectory and it's not empty, we descend into it.
    Otherwise, we treat workspace_dir as the docking root.

    This allows tarballs that either contain:
      - a single top-level folder (e.g. '5wlo/...')
      - or files directly at the root.
    """
    # List immediate children
    children = [p for p in workspace_dir.iterdir() if p.is_dir() or p.is_file()]

    # If there is exactly one child and it is a directory with content, use that
    dirs = [p for p in children if p.is_dir()]
    if len(dirs) == 1:
        # Check if this dir has some contents
        if any(True for _ in dirs[0].iterdir()):
            return dirs[0]

    # Otherwise, just treat the workspace_dir as docking root
    return workspace_dir


def _validate_rel_path(name: str, field: str) -> None:
    """
    Extra API-side validation for receptor/ligand relative paths.
    """
    p = Path(name)
    if p.is_absolute():
        raise HTTPException(
            status_code=400,
            detail=f"{field} must be a relative path inside the dataset, not an absolute path.",
        )
    if ".." in p.parts:
        raise HTTPException(
            status_code=400,
            detail=f"{field} must not contain '..' path components.",
        )


@OpeaComponentRegistry.register("OPENOMICS_AUTODOCK_VINA")
class OpeaAutodockVina(OpeaComponent):
    """
    OPEA component that wraps an AutoDock Vina docking run.

    Flow:
      - Create per-user/per-job workspace directory
      - Decode & safely extract dataset_tar_b64 into it
      - Locate receptor & ligand files
      - Build 'vina --receptor ... --ligand ... --center_x ... --size_x ...' command
      - Run Vina
      - Return stdout/stderr and output file (optionally as base64 tar.gz)
      - Delete the job directory after the response is prepared
    """

    def __init__(self, name: str, description: str, config: dict | None = None):
        config = config or {}
        # workload name: "autodock_vina" (for classification in OPEA)
        super().__init__(name, "autodock_vina", description, config)

        # Workspace root for all jobs (configurable)
        workspace_root_str = config.get("workspace_root", "/workspaces")
        self.base_workspace_dir = Path(workspace_root_str)
        self.base_workspace_dir.mkdir(parents=True, exist_ok=True)

        # Max tar size (bytes) – optional override
        self.max_tar_bytes = int(config.get("max_tar_bytes", MAX_TAR_BYTES))

        logger.info(
            f"OpeaAutodockVina initialized: workspace_root={self.base_workspace_dir}, "
            f"max_tar_bytes={self.max_tar_bytes}"
        )

    # ------------- main invoke -------------
    async def invoke(self, input: VinaInput) -> VinaOutput:
        """
        Handle a single Vina docking request.

        - Validates & extracts client dataset
        - Constructs safe Vina command using JSON flags
        - Executes Vina (no internal timeout)
        - Returns stdout / stderr / output path
        - Always deletes the per-job workspace directory after preparing the response
        """
        start_time = time.time()
        workspace_dir: Optional[Path] = None
        workspace_id: Optional[str] = None
        job_id: Optional[str] = None

        # --- Basic scoring restriction ---
        if input.scoring and input.scoring.strip().lower() != "vina":
            logger.error(f"Unsupported scoring function requested: {input.scoring}")
            raise HTTPException(
                status_code=400,
                detail="Only 'vina' scoring function is supported by this service.",
            )

        # Extra API-side validation of receptor / ligand paths
        _validate_rel_path(input.receptor, "receptor")
        _validate_rel_path(input.ligand, "ligand")

        # --- Build workspace directory: /<root>/<user>/<workspace>/<job> ---
        try:
            workspace_dir, workspace_id, job_id = build_workspace_dirs(
                self.base_workspace_dir,
                input.user_id,
                input.workspace_id,
            )
        except ValueError as e:
            logger.error(f"Invalid workspace/user id: {e}")
            raise HTTPException(status_code=400, detail=str(e))

        logger.info(
            f"Starting Vina job: user={input.user_id}, "
            f"workspace_id={workspace_id}, job_id={job_id}, dir={workspace_dir}"
        )

        svc_name = "opea_service@autodock_vina"

        try:
            # --- Decode & extract dataset tar.gz ---
            try:
                tar_bytes = decode_tar_b64(input.dataset_tar_b64, max_bytes=self.max_tar_bytes)
                safe_extract_tar_bytes(tar_bytes, workspace_dir)
            except ValueError as e:
                logger.error(f"Dataset decode/extract error: {e}")
                raise HTTPException(status_code=400, detail=f"Invalid dataset: {e}")

            # --- Determine docking root directory ---
            docking_root = find_docking_root(workspace_dir)

            # --- Resolve file paths within docking root ---
            receptor_path = (docking_root / input.receptor).resolve()
            ligand_path = (docking_root / input.ligand).resolve()

            # Ensure they are inside workspace_root for safety
            if self.base_workspace_dir not in receptor_path.parents:
                raise HTTPException(
                    status_code=400,
                    detail="Receptor path escapes workspace root.",
                )
            if self.base_workspace_dir not in ligand_path.parents:
                raise HTTPException(
                    status_code=400,
                    detail="Ligand path escapes workspace root.",
                )

            if not receptor_path.is_file():
                msg = f"Receptor file not found: {input.receptor}"
                logger.error(msg)
                raise HTTPException(status_code=400, detail=msg)

            if not ligand_path.is_file():
                msg = f"Ligand file not found: {input.ligand}"
                logger.error(msg)
                raise HTTPException(status_code=400, detail=msg)

            # --- Build Vina command from JSON flags ---
            cmd: list[str] = [
                "vina",
                "--receptor",
                str(receptor_path),
                "--ligand",
                str(ligand_path),
                "--center_x",
                str(input.center_x),
                "--center_y",
                str(input.center_y),
                "--center_z",
                str(input.center_z),
                "--size_x",
                str(input.size_x),
                "--size_y",
                str(input.size_y),
                "--size_z",
                str(input.size_z),
                "--scoring",
                "vina",  # forced
            ]

            # Output options: kept relative to docking_root
            output_pdbqt_path = None
            if input.out:
                output_pdbqt_path = (docking_root / input.out).resolve()
                cmd += ["--out", str(output_pdbqt_path)]

            # Misc options
            if input.cpu:
                cmd += ["--cpu", str(input.cpu)]
            if input.seed is not None:
                cmd += ["--seed", str(input.seed)]
            if input.exhaustiveness is not None:
                cmd += ["--exhaustiveness", str(input.exhaustiveness)]

            cmd_str = " ".join(cmd)
            logger.info(
                f"[user={input.user_id} workspace={workspace_id} job={job_id}] "
                f"Running Vina command: {cmd_str}"
            )

            # --- Run Vina (no timeout) ---
            proc = subprocess.run(
                cmd,
                cwd=str(docking_root),
                capture_output=True,
                text=True,
            )
            elapsed = time.time() - start_time
            stdout = proc.stdout or ""
            stderr = proc.stderr or ""

            if proc.returncode != 0:
                logger.error(
                    f"Vina failed with code {proc.returncode}. stderr:\n{stderr}"
                )
                if svc_name in statistics_dict:
                    statistics_dict[svc_name].append_latency(elapsed, None)

                # FULL stdout/stderr back to client
                raise HTTPException(
                    status_code=400,
                    detail={
                        "status": "error",
                        "message": f"Vina failed with code {proc.returncode}",
                        "stdout": stdout,
                        "stderr": stderr,
                        "user_id": input.user_id,
                        "workspace_id": workspace_id,
                        "job_id": job_id,
                    },
                )

            logger.info(
                f"Vina job completed successfully: user={input.user_id}, "
                f"workspace_id={workspace_id}, job_id={job_id}, elapsed={elapsed:.2f}s"
            )
            if svc_name in statistics_dict:
                statistics_dict[svc_name].append_latency(elapsed, None)

            # If user did not specify 'out', try to infer default name
            if output_pdbqt_path is None:
                ligand_stem = Path(input.ligand).stem
                guess = (docking_root / f"{ligand_stem}_out.pdbqt").resolve()
                if guess.is_file():
                    output_pdbqt_path = guess

            result_tar_b64 = None
            result_tar_size = None
            if input.return_tarball:
                try:
                    if output_pdbqt_path is None or not output_pdbqt_path.is_file():
                        logger.error(
                            "Output tarball requested but output PDBQT file does not exist."
                        )
                    else:
                        metadata = {
                            "user_id": input.user_id,
                            "workspace_id": workspace_id,
                            "job_id": job_id,
                            "runtime_seconds": elapsed,
                            "output_file": output_pdbqt_path.name,
                            "created_at": time.strftime(
                                "%Y-%m-%dT%H:%M:%SZ", time.gmtime()
                            ),
                            "vina_command": cmd_str,
                        }
                        result_tar_b64, result_tar_size = create_output_tar_b64_with_metadata(
                            output_pdbqt_path,
                            metadata,
                            max_bytes=MAX_RESULT_TAR_BYTES,
                        )
                        logger.info(
                            f"Created output-only tar.gz for user={input.user_id}, "
                            f"workspace_id={workspace_id}, job_id={job_id}, "
                            f"file={output_pdbqt_path.name}, size={result_tar_size} bytes"
                        )
                except ValueError as e:
                    # If tarball is too large or fails, we log it but do not fail the docking.
                    logger.error(f"Failed to create result tarball: {e}")

            response = VinaOutput(
                status="success",
                message="Vina docking completed successfully.",
                user_id=input.user_id,
                workspace_id=workspace_id,
                job_id=job_id,
                output_pdbqt=str(output_pdbqt_path.relative_to(self.base_workspace_dir))
                if output_pdbqt_path
                else None,
                runtime_seconds=elapsed,
                stdout=stdout,
                stderr=stderr or None,
                result_tar_b64=result_tar_b64,
                result_tar_size_bytes=result_tar_size,
            )
            return response

        except HTTPException:
            # re-raise HTTPExceptions so FastAPI returns proper error codes
            raise
        except Exception as e:
            logger.exception(
                f"Unexpected error during Vina job: user={getattr(input, 'user_id', None)}, "
                f"workspace_id={workspace_id}, job_id={job_id}: {e}"
            )
            raise HTTPException(
                status_code=500,
                detail="Internal server error while running Vina.",
            )
        finally:
            # Cleanup workspace directory on success/failure AFTER building response.
            # We also try to clean empty parent directories up to base_workspace_dir.
            if workspace_dir is not None:
                try:
                    if workspace_dir.exists():
                        shutil.rmtree(workspace_dir, ignore_errors=True)
                        logger.info(
                            f"Cleaned up job workspace_dir={workspace_dir} for user={input.user_id}, "
                            f"workspace_id={workspace_id}, job_id={job_id}"
                        )
                    # optional: clean up empty parents (workspace and user dirs)
                    parent = workspace_dir.parent
                    while parent != self.base_workspace_dir and parent.exists():
                        try:
                            next(parent.iterdir())
                            break  # not empty
                        except StopIteration:
                            # empty -> remove and go up
                            parent.rmdir()
                            parent = parent.parent
                except Exception as e:
                    logger.error(f"Failed to cleanup workspace_dir={workspace_dir}: {e}")

    # ------------- health check -------------
    async def check_health(self) -> dict:
        return {"status": "healthy"}


# Global component instance (set from server at startup)
_vina_component: OpeaAutodockVina | None = None


def init_vina_component(config: dict) -> None:
    """
    Initialize the global AutoDock Vina component with the provided config.
    Called from the server at startup.
    """
    global _vina_component
    _vina_component = OpeaAutodockVina(
        name="opea_service@autodock_vina",
        description="AutoDock Vina microservice component",
        config=config,
    )

@register_statistics(names=["opea_service@autodock_vina"])
async def autodock_vina_service(input: VinaInput) -> VinaOutput:
    """
    Service function used by register_microservice. It delegates to the
    global OpeaAutodockVina component instance.
    """
    start = time.time()

    if _vina_component is None:
        logger.error("Vina component not initialized")
        raise HTTPException(status_code=500, detail="Vina component not initialized")

    try:
        return await _vina_component.invoke(input)
    finally:
        if "opea_service@autodock_vina" in statistics_dict:
            statistics_dict["opea_service@autodock_vina"].append_latency(
                time.time() - start, None
            )

