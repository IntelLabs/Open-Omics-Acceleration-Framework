import base64
import json
import os
import re
import secrets
import shutil
import subprocess
import time
from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from pydantic import BaseModel, Field, constr, conint

from comps import (
    CustomLogger,
    OpeaComponent,
    OpeaComponentRegistry,
)

logger = CustomLogger("gromacs_microservice")


class GromacsInput(BaseModel):
    """
    GROMACS microservice input schema (workflow-only).

    Workflow:
      - PDB is provided by the client as base64 (pdb_b64) and a filename (pdb_filename).
      - Service copies run_commands.sh and mdtut_*.mdp from /input into a per-job workspace.
      - Service writes the PDB file into that workspace.
      - Then runs: ./run_commands.sh <pdb_filename>
    """

    user_id: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description="Logical user identifier for workspace isolation.",
    )
    job_label: Optional[constr(strip_whitespace=True, min_length=1)] = Field(
        None,
        description="Optional human-friendly label for this job.",
    )

    pdb_filename: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description=(
            "Desired PDB filename inside the workspace (e.g. 'protein.pdb'). "
            "This name will be passed as argument to run_commands.sh."
        ),
    )
    pdb_b64: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description="Base64-encoded contents of the PDB file.",
    )

    env: Dict[str, str] = Field(
        default_factory=dict,
        description="Extra environment variables for this run.",
    )


class Artifact(BaseModel):
    name: str
    b64: str
    encoding: str = "base64"
    mime: str = "application/json"


class GromacsOutput(BaseModel):
    status: str
    message: Optional[str] = None
    artifacts: Optional[List[Artifact]] = None
    mode: Optional[str] = None
    runtime_sec: Optional[float] = None
    workspace: Optional[str] = None
    commands: Optional[List[str]] = None
    files: Optional[List[str]] = None
    job_id: Optional[str] = None
    user_id: Optional[str] = None
    job_label: Optional[str] = None


def _ensure_paths() -> None:
    for p in ("/input", "/output", "/input/jobs"):
        os.makedirs(p, exist_ok=True)


def _safe_user_id(user_id: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_.-]", "_", user_id)


def _mk_user_jobdir(user_id: str, job_label: Optional[str]) -> str:
    _ensure_paths()
    safe_uid = _safe_user_id(user_id)
    user_root = os.path.join("/input/jobs", safe_uid)
    os.makedirs(user_root, exist_ok=True)

    ts = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    suffix = secrets.token_hex(4)
    base = f"{safe_uid}_{ts}_{suffix}"
    if job_label:
        safe_label = _safe_user_id(job_label)
        base = f"{safe_uid}_{safe_label}_{ts}_{suffix}"

    job_dir = os.path.join(user_root, base)
    os.makedirs(job_dir, exist_ok=True)
    return job_dir


def _write_pdb_from_b64(pdb_b64: str, dest_path: str) -> None:
    """
    Decode base64 PDB content and write to dest_path.
    """
    try:
        raw = base64.b64decode(pdb_b64)
    except Exception as e:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid pdb_b64 (base64 decode failed): {e}",
        )
    try:
        with open(dest_path, "wb") as f:
            f.write(raw)
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to write PDB file to {dest_path}: {e}",
        )


def _list_workspace_files(dirpath: str) -> List[str]:
    files: List[str] = []
    for root, dirs, filenames in os.walk(dirpath):
        for name in filenames:
            full = os.path.join(root, name)
            rel = os.path.relpath(full, dirpath)
            files.append(rel)
    files.sort()
    return files


def _cleanup_dir(path: str) -> None:
    try:
        shutil.rmtree(path, ignore_errors=True)
        logger.info(f"Cleaned job directory: {path}")
    except Exception as e:
        logger.warning(f"Failed to cleanup job directory {path}: {e}")


def _run(
    cmd: List[str],
    env: Dict[str, str],
    timeout_sec: Optional[int],
    cwd: Optional[str] = None,
    stdin_text: Optional[str] = None,
) -> tuple[int, str, str]:
    """
    Run a subprocess with merged env and optional timeout.

    If timeout_sec is None -> no timeout is enforced.
    """
    logger.info(f"Exec (cwd={cwd or os.getcwd()}): {' '.join(cmd)}")
    merged_env = os.environ.copy()
    merged_env.update(env or {})

    p = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE if stdin_text is not None else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=merged_env,
        cwd=cwd,
    )
    try:
        if timeout_sec is None:
            out, err = p.communicate(input=stdin_text)
        else:
            out, err = p.communicate(input=stdin_text, timeout=timeout_sec)
        return p.returncode, out, err
    except subprocess.TimeoutExpired:
        p.kill()
        raise HTTPException(
            status_code=504,
            detail=f"Command timed out after {timeout_sec}s: {' '.join(cmd)}",
        )

@OpeaComponentRegistry.register("OPENOMICS_GROMACS")
class OpeaGromacs(OpeaComponent):
    def __init__(
        self,
        name: str,
        description: str,
        config: Dict[str, Any] | None = None,
    ):
        super().__init__(name, "gromacs", description, config or {})
        _ensure_paths()
        logger.info("GROMACS microservice initialized")

    async def invoke(self, input: GromacsInput) -> GromacsOutput:
        if not input:
            raise HTTPException(status_code=400, detail="Missing input payload")

        start = time.time()

        job_dir = _mk_user_jobdir(input.user_id, input.job_label)
        executed_cmds: List[str] = []
        job_id = os.path.basename(job_dir)

        code = 0
        out = err = ""
        stdin_text: Optional[str] = None

        try:
            # 1) Stage workflow assets: run_commands.sh + mdtut_*.mdp
            src_run = os.path.join("/input", "run_commands.sh")
            if not os.path.exists(src_run):
                raise HTTPException(
                    status_code=500,
                    detail="run_commands.sh not found under /input for workflow mode",
                )
            dst_run = os.path.join(job_dir, "run_commands.sh")
            shutil.copy2(src_run, dst_run)
            os.chmod(dst_run, 0o755)

            for f in os.listdir("/input"):
                if f.startswith("mdtut_") and f.endswith(".mdp"):
                    src_mdp = os.path.join("/input", f)
                    shutil.copy2(src_mdp, job_dir)

            # 2) Write the PDB file from base64
            pdb_name = os.path.basename(input.pdb_filename)
            pdb_dst = os.path.join(job_dir, pdb_name)
            _write_pdb_from_b64(input.pdb_b64, pdb_dst)

            # 3) Build and run the workflow command
            cmd = ["./run_commands.sh", pdb_name]
            executed_cmds.append(" ".join(cmd))

            env = {**input.env, "GROMACS_WRKDIR": job_dir}
            code, out, err = _run(
                cmd,
                env,
                timeout_sec=None,
                cwd=job_dir,
                stdin_text=stdin_text,
            )
            logger.info(f"Exit code: {code}")

            # 4) Discover files and build metadata
            workspace_files = _list_workspace_files(job_dir)

            metadata: Dict[str, Any] = {
                "status": "success" if code == 0 else "error",
                "mode": "workflow",
                "job_id": job_id,
                "user_id": input.user_id,
                "job_label": input.job_label,
                "commands": executed_cmds,
                "workspace_files": workspace_files,
                "exit_code": code,
                "stdout_snip": (out or "")[:2000],
                "stderr_snip": (err or "")[:2000],
            }

            runtime = time.time() - start
            metadata["runtime_sec"] = runtime

            # 5) On error: raise with metadata
            if code != 0:
                raise HTTPException(
                    status_code=500,
                    detail=metadata,
                )

            # 6) Prepare metadata.json as the ONLY inline artifact
            artifacts: List[Artifact] = []
            meta_bytes = json.dumps(metadata, indent=2).encode("utf-8")
            artifacts.append(
                Artifact(
                    name=f"metadata_{job_id}.json",
                    b64=base64.b64encode(meta_bytes).decode("utf-8"),
                    encoding="base64",
                    mime="application/json",
                )
            )

            msg = (
                f"Workflow finished in {runtime:.2f}s (rc={code}). "
                f"job_id={job_id}"
            )

            return GromacsOutput(
                status="success",
                message=msg,
                artifacts=artifacts or None,
                mode="workflow",
                runtime_sec=runtime,
                workspace=job_dir,
                commands=executed_cmds,
                files=workspace_files,
                job_id=job_id,
                user_id=input.user_id,
                job_label=input.job_label,
            )

        except HTTPException:
            raise
        except Exception as e:
            logger.error(f"Unexpected error in GROMACS invoke: {e}")
            raise HTTPException(
                status_code=500,
                detail=str(e),
            )
        finally:
            _cleanup_dir(job_dir)

    def check_health(self) -> bool:
        try:
            code, out, err = _run(
                ["gmx", "-version"],
                env={},
                timeout_sec=15,
                cwd=None,
            )
            return code == 0
        except Exception as e:
            logger.error(f"Health check failed: {e}")
            return False

