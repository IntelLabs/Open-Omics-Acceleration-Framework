#!/usr/bin/env python3
import argparse
import base64
import json
from pathlib import Path
from typing import Any, Dict
import requests

def _encode_pdb(path: str) -> tuple[str, str]:
    """
    Read PDB file from disk and return (pdb_filename, pdb_b64).
    """
    p = Path(path)
    if not p.is_file():
        raise SystemExit(f"PDB file not found: {p}")

    pdb_filename = p.name
    raw = p.read_bytes()
    pdb_b64 = base64.b64encode(raw).decode("utf-8")
    return pdb_filename, pdb_b64


def _print_basic_response(res_json: Dict[str, Any]) -> None:
    print("Message      :", res_json.get("message"))
    print("Status       :", res_json.get("status"))
    print("Mode         :", res_json.get("mode"))
    print("Runtime (s)  :", res_json.get("runtime_sec"))
    print("User ID      :", res_json.get("user_id"))
    print("Job ID       :", res_json.get("job_id"))
    print("Job label    :", res_json.get("job_label"))
    print("Workspace    :", res_json.get("workspace"))

    cmds = res_json.get("commands") or []
    print("Commands:")
    for c in cmds:
        print("  -", c)

    files = res_json.get("files") or []
    print("Files in workspace:")
    for f in files:
        print("  -", f)


def _save_metadata_artifacts(res_json: Dict[str, Any]) -> None:
    """
    Save any artifacts returned by the service.

    In the current design, we expect exactly one artifact:
      - metadata_<job_id>.json
    """
    arts = res_json.get("artifacts") or []
    if not arts:
        print("No artifacts in response.")
        return

    print("Artifacts:")
    for i, a in enumerate(arts, 1):
        name = a.get("name") or f"artifact_{i}"
        encoding = a.get("encoding") or "base64"
        b64 = a.get("b64")

        if encoding != "base64" or not b64:
            # treat as text (should not happen in current design)
            text = b64 or ""
            Path(name).write_text(text)
            print(f"  - Saved text artifact: {name}")
            continue

        raw = base64.b64decode(b64)
        Path(name).write_bytes(raw)
        print(f"  - Saved binary artifact: {name} ({len(raw)} bytes)")


def main():
    parser = argparse.ArgumentParser(
        description="Client for OPENOMICS GROMACS microservice (workflow-only)."
    )
    parser.add_argument("--host", default="127.0.0.1", help="Microservice host")
    parser.add_argument("--port", default="9000", help="Microservice port")
    parser.add_argument(
        "--user-id",
        required=True,
        help="Logical user identifier (used for per-user job directories).",
    )
    parser.add_argument(
        "--job-label",
        help="Optional human-friendly label for this job.",
    )

    parser.add_argument(
        "--pdb-file",
        required=True,
        help="Path to local PDB file to send to the workflow.",
    )
    parser.add_argument(
        "--pdb-filename",
        help=(
            "Logical PDB filename inside the workspace (e.g. 'protein.pdb'). "
            "If omitted, the basename of --pdb-file is used."
        ),
    )

    args = parser.parse_args()

    url = f"http://{args.host}:{args.port}/v1/gromacs"
    default_name, pdb_b64 = _encode_pdb(args.pdb_file)
    pdb_filename = args.pdb_filename or default_name

    payload: Dict[str, Any] = {
        "user_id": args.user_id,
        "job_label": args.job_label,
        "pdb_filename": pdb_filename,
        "pdb_b64": pdb_b64,
    }

    print("POST", url)
    response = requests.post(url, json=payload)
    print("Status:", response.status_code)

    if response.status_code != 200:
        try:
            err = response.json()
            print("Error JSON:", json.dumps(err, indent=2))
        except Exception:
            print("Error text:", response.text)
        raise SystemExit(1)

    res = response.json()
    _print_basic_response(res)
    _save_metadata_artifacts(res)


if __name__ == "__main__":
    main()

