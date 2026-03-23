import base64
import io
import json
import shutil
import subprocess
import tarfile
import time
import uuid
from pathlib import Path
from typing import Any, Dict, Optional, Literal

from fastapi import HTTPException
from pydantic import BaseModel, Field, conint, constr, field_validator

from comps import CustomLogger, OpeaComponent, OpeaComponentRegistry

logger = CustomLogger("autodock_microservice")

# ---------------------------------------------------------------------
# STANDARD API ERROR SCHEMA
# ---------------------------------------------------------------------


class ApiError(BaseModel):
    """
    Standard error payload for all HTTPException responses.
    """

    error_type: str
    message: str
    details: Optional[Dict[str, Any]] = None


def _api_error(
    status_code: int,
    message: str,
    error_type: str = "BadRequest",
    details: Optional[Dict[str, Any]] = None,
) -> HTTPException:
    payload = ApiError(error_type=error_type, message=message, details=details).model_dump()
    return HTTPException(status_code=status_code, detail=payload)


# ---------------------------------------------------------------------
# INPUT / OUTPUT MODELS
# ---------------------------------------------------------------------


class AutoDockInput(BaseModel):
    """
    Input schema for AutoDock-GPU microservice.
    """

    # Multi-user / workspace
    user_id: constr(strip_whitespace=True, min_length=1) = Field(
        ..., description="Logical user identifier (no spaces)."
    )
    workspace_id: Optional[constr(strip_whitespace=True, min_length=1)] = Field(
        None,
        description="Optional workspace group id; defaults to 'default'.",
    )

    # Dataset (.tar.gz -> base64)
    dataset_tar_b64: str = Field(
        ...,
        description="Base64-encoded .tar.gz of the docking input folder.",
    )

    # Paths relative to the dataset folder (e.g. protein.maps.fld, rand-0.pdbqt)
    ffile: str = Field(
        ...,
        description="Grid descriptor FLD file (for --ffile / -M), relative to dataset folder.",
    )
    lfile: str = Field(
        ...,
        description="Ligand PDBQT file (for --lfile / -L), relative to dataset folder.",
    )

    # AutoDock flags
    nrun: conint(gt=0) = Field(
        20, description="# LGA runs (for --nrun / -n)."
    )
    lsmet: Literal["ad", "sw"] = Field(
        "ad", description="Local search method (for --lsmet / -l)."
    )
    nev: conint(gt=0) = Field(
        2_500_000,
        description="# Score evaluations max per LGA run (for --nev / -e).",
    )

    seed: Optional[str] = Field(
        None,
        description="One to three comma-separated integers for --seed (e.g. '11,23').",
    )

    # resnam is REQUIRED (used to select output files)
    resnam: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description=(
            "Base name for docking outputs (for --resnam / -N). "
            "Service will return <resnam>.dlg and <resnam>.xml."
        ),
    )

    @field_validator("seed")
    @classmethod
    def _validate_seed(cls, v: Optional[str]) -> Optional[str]:
        if v is None:
            return v
        parts = v.split(",")
        if not 1 <= len(parts) <= 3:
            raise ValueError("seed must be 1–3 comma-separated integers, e.g. '11,23'.")
        for p in parts:
            int(p)  # raises if not int
        return v


class AutoDockOutput(BaseModel):
    """
    Response schema:

    - On success:
        status = "success"
        message = summary
        metadata_b64 = base64(JSON metadata.json)
        results_tar_b64 = base64(tar.gz of [resnam.dlg, resnam.xml, metadata.json])
        error = null

    - On failure:
        status = "failure"
        message = summary
        error = full stderr string from AutoDock
        metadata_b64 = null
        results_tar_b64 = null
    """

    status: Literal["success", "failure"]
    message: str
    metadata_b64: Optional[str] = None
    results_tar_b64: Optional[str] = None
    error: Optional[str] = None


# ---------------------------------------------------------------------
# UTILS
# ---------------------------------------------------------------------

BASE_WORKDIR = Path("/tmp/autodock_inputs")

def _safe_extract_tar_gz(tar_bytes: bytes, dest_dir: Path) -> None:
    """Safely extract a tar.gz archive into dest_dir, avoiding path traversal."""
    try:
        dest_dir.mkdir(parents=True, exist_ok=True)
        buf = io.BytesIO(tar_bytes)

        with tarfile.open(fileobj=buf, mode="r:gz") as tf:
            for member in tf.getmembers():
                dest_root = dest_dir.resolve()
                member_path = dest_dir / member.name
                if dest_root not in member_path.parents and member_path != dest_root:
                    raise _api_error(
                        400,
                        message=f"Unsafe path in archive: {member.name}",
                        error_type="BadRequest",
                    )
            tf.extractall(path=dest_dir)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Tar extract failed: {e}")
        raise _api_error(
            400,
            message="Invalid docking dataset archive",
            error_type="BadRequest",
            details={"reason": str(e)},
        )


def _decode_dataset_tar(dataset_tar_b64: str, work_dir: Path) -> None:
    """Decode base64 → tar.gz → extract."""
    try:
        raw = base64.b64decode(dataset_tar_b64)
    except Exception:
        raise _api_error(
            400,
            message="Invalid base64 in dataset_tar_b64",
            error_type="BadRequest",
        )
    _safe_extract_tar_gz(raw, work_dir)


def _get_dataset_root_dir(work_dir: Path) -> Path:
    """
    Detect dataset root directory under work_dir.

    Typical case:
      work_dir/
        1mzc/      <-- dataset_root
          protein.maps.fld
          rand-0.pdbqt
          ...

    If exactly one subdirectory exists: use that.
    Otherwise: fall back to work_dir itself.
    """
    subdirs = [p for p in work_dir.iterdir() if p.is_dir()]
    if len(subdirs) == 1:
        dataset_root = subdirs[0]
        logger.info(f"Detected dataset root directory: {dataset_root}")
        return dataset_root

    logger.warning(
        f"Could not uniquely detect dataset root under {work_dir}, "
        f"using {work_dir} as dataset root."
    )
    return work_dir


# ---------------------------------------------------------------------
# AUTODOCK EXECUTION
# ---------------------------------------------------------------------

def _build_autodock_command(input_obj: AutoDockInput, dataset_root: Path) -> Any:
    """
    Build autodock_cpu_64wi command and verify ffile/lfile exist in dataset_root.
    """
    ffile_path = dataset_root / input_obj.ffile
    lfile_path = dataset_root / input_obj.lfile

    if not ffile_path.is_file():
        raise _api_error(
            400,
            message=f"--ffile not found in dataset folder: {input_obj.ffile}",
            error_type="BadRequest",
        )
    if not lfile_path.is_file():
        raise _api_error(
            400,
            message=f"--lfile not found in dataset folder: {input_obj.lfile}",
            error_type="BadRequest",
        )

    cmd = [
        "autodock_cpu_64wi",
        "--ffile",
        input_obj.ffile,
        "--lfile",
        input_obj.lfile,
        "--nrun",
        str(input_obj.nrun),
        "--lsmet",
        input_obj.lsmet,
        "--nev",
        str(input_obj.nev),
        "--resnam",
        input_obj.resnam,
    ]

    if input_obj.seed:
        cmd.extend(["--seed", input_obj.seed])

    return cmd


def _run_autodock(cmd: Any, cwd: Path) -> Dict[str, Any]:
    """
    Run AutoDock-GPU with no internal timeout; return metadata dict (in-memory).
    cwd is the dataset_root (e.g. /tmp/.../1mzc).
    """
    logger.info(f"Running in {cwd}: {' '.join(cmd)}")
    start = time.time()

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            capture_output=True,
            text=True,
        )
    except Exception as e:
        logger.error(f"Execution failed: {e}")
        raise _api_error(
            500,
            message="Failed to execute AutoDock-GPU",
            error_type="ServerError",
            details={"reason": str(e)},
        )

    runtime = time.time() - start

    metadata: Dict[str, Any] = {
        "command": cmd,
        "cwd": str(cwd),
        "exit_code": proc.returncode,
        "runtime_sec": runtime,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "output_files": [],
    }

    for p in cwd.rglob("*"):
        if p.is_file():
            metadata["output_files"].append(str(p.relative_to(cwd)))

    return metadata


def _write_metadata_and_encode(work_dir: Path, metadata: Dict[str, Any]) -> str:
    """
    Write metadata.json in work_dir and return base64(JSON).
    Called only on success.
    """
    try:
        metadata_path = work_dir / "metadata.json"
        with metadata_path.open("w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2, default=str)

        raw = metadata_path.read_bytes()
        return base64.b64encode(raw).decode("utf-8")
    except Exception as e:
        logger.error(f"Failed to write/encode metadata.json: {e}")
        raise _api_error(
            500,
            message="Failed to write metadata.json",
            error_type="ServerError",
            details={"reason": str(e)},
        )


def _create_results_tar_b64(work_dir: Path, dataset_root: Path, resnam: str) -> str:
    """
    Create tar.gz (base64) containing ONLY:
      - <resnam>.dlg       (from dataset_root)
      - <resnam>.xml       (from dataset_root)
      - metadata.json      (from work_dir)
    Called only on success.
    """
    try:
        dlg_path = dataset_root / f"{resnam}.dlg"
        xml_path = dataset_root / f"{resnam}.xml"
        metadata_path = work_dir / "metadata.json"

        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tf:
            for p in (dlg_path, xml_path, metadata_path):
                if p.is_file():
                    tf.add(p, arcname=p.name)

        buf.seek(0)
        raw = buf.read()
        return base64.b64encode(raw).decode("utf-8")
    except Exception as e:
        logger.error(f"Failed to create results tar.gz: {e}")
        raise _api_error(
            500,
            message="Failed to create results archive",
            error_type="ServerError",
            details={"reason": str(e)},
        )


# ---------------------------------------------------------------------
# OPEA COMPONENT
# ---------------------------------------------------------------------

@OpeaComponentRegistry.register("OPEA_OMICS_AUTODOCK")
class Opea_AutoDock(OpeaComponent):
    """
    AutoDock-GPU microservice component.
    """

    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config or {})
        logger.info("Opea_AutoDock initialized")

    async def check_health(self) -> dict:
        return {"status": "healthy"}

    async def invoke(self, input: AutoDockInput) -> AutoDockOutput:
        """
        One invocation = one job in its own workspace directory.
        """
        request_id = str(uuid.uuid4())
        effective_workspace = input.workspace_id or "default"
        work_dir = BASE_WORKDIR / input.user_id / effective_workspace / request_id

        logger.info(
            f"New request: user_id={input.user_id}, workspace={effective_workspace}, "
            f"request_id={request_id}, work_dir={work_dir}"
        )

        work_dir.mkdir(parents=True, exist_ok=True)

        try:
            # 1) Extract dataset into work_dir
            _decode_dataset_tar(input.dataset_tar_b64, work_dir)

            # 2) Detect dataset root directory (e.g. work_dir/1mzc)
            dataset_root = _get_dataset_root_dir(work_dir)

            # 3) Build command (validate ffile/lfile inside dataset_root)
            cmd = _build_autodock_command(input, dataset_root)

            # 4) Run docking (no timeout), with cwd = dataset_root
            metadata = _run_autodock(cmd, dataset_root)

            # 5) Handle success vs failure
            if metadata["exit_code"] != 0:
                error_str = metadata.get("stderr") or metadata.get("stdout") or ""
                msg = (
                    "AutoDock-GPU docking failed. "
                    "See 'error' field for full stdout/stderr."
                )
                return AutoDockOutput(
                    status="failure",
                    message=msg,
                    error=error_str,
                    metadata_b64=None,
                    results_tar_b64=None,
                )

            # SUCCESS path
            metadata_b64 = _write_metadata_and_encode(work_dir, metadata)
            results_tar_b64 = _create_results_tar_b64(work_dir, dataset_root, input.resnam)

            msg = (
                "AutoDock-GPU docking completed successfully. "
                f"Runtime: {metadata['runtime_sec']:.2f} sec."
            )
            return AutoDockOutput(
                status="success",
                message=msg,
                metadata_b64=metadata_b64,
                results_tar_b64=results_tar_b64,
                error=None,
            )

        except HTTPException:
            raise
        except Exception as e:
            logger.error(f"Unexpected error in AutoDock-GPU microservice: {e}")
            raise _api_error(
                500,
                message="Unexpected server error",
                error_type="ServerError",
                details={"reason": str(e)},
            )
        finally:
            try:
                if work_dir.exists():
                    shutil.rmtree(work_dir)
                    logger.info(f"Deleted workspace: {work_dir}")
            except Exception as e:
                logger.error(f"Failed to delete workspace {work_dir}: {e}")

