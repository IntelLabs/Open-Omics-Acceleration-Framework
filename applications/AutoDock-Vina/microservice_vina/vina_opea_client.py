#!/usr/bin/env python

import argparse
import base64
import json
import tarfile
from pathlib import Path

import requests


def create_tar_b64(dataset_dir: Path) -> str:
    """
    Create a tar.gz archive of dataset_dir and return it as base64-encoded string.
    """
    if not dataset_dir.is_dir():
        raise ValueError(f"Dataset directory does not exist or is not a directory: {dataset_dir}")

    # Create tar.gz in memory
    from io import BytesIO

    buffer = BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
        # Put the directory into the tar with its basename as top-level
        tar.add(dataset_dir, arcname=dataset_dir.name)

    buffer.seek(0)
    raw_bytes = buffer.read()
    return base64.b64encode(raw_bytes).decode("ascii")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Client for OPEA AutoDock Vina microservice."
    )

    # ---- Microservice location ----
    parser.add_argument("--host", default="127.0.0.1", help="Microservice host (default: 127.0.0.1)")
    parser.add_argument("--port", default="9000", help="Microservice port (default: 9000)")

    # ---- User / workspace ----
    parser.add_argument("--user-id", required=True, help="Logical user identifier, e.g. 'sri'")
    parser.add_argument(
        "--workspace-id",
        default=None,
        help="Optional workspace identifier; if omitted, server generates one.",
    )

    # ---- Dataset directory ----
    parser.add_argument(
        "--dataset-dir",
        required=True,
        help="Path to local docking dataset directory (e.g. '5wlo'). "
             "This folder will be tar.gz'ed and sent to the server.",
    )

    # ---- Required Vina files inside dataset ----
    parser.add_argument(
        "--receptor",
        required=True,
        help="Receptor PDBQT file name inside dataset (e.g. 'protein.pdbqt').",
    )
    parser.add_argument(
        "--ligand",
        required=True,
        help="Ligand PDBQT file name inside dataset (e.g. 'rand-1.pdbqt').",
    )

    # ---- Tarball download options ----
    parser.add_argument(
        "--download-tarball",
        action="store_true",
        help="If set, request a tar.gz of the remote job output+metadata and save it locally.",
    )
    parser.add_argument(
        "--tarball-output",
        default="vina_job_output.tar.gz",
        help="Filename to save the returned tarball (default: vina_job_output.tar.gz).",
    )

    # ---- Search space (required) ----
    parser.add_argument("--center-x", type=float, required=True, help="Center X (Angstrom).")
    parser.add_argument("--center-y", type=float, required=True, help="Center Y (Angstrom).")
    parser.add_argument("--center-z", type=float, required=True, help="Center Z (Angstrom).")
    parser.add_argument("--size-x", type=float, required=True, help="Box size X (Angstrom).")
    parser.add_argument("--size-y", type=float, required=True, help="Box size Y (Angstrom).")
    parser.add_argument("--size-z", type=float, required=True, help="Box size Z (Angstrom).")

    # ---- Output options ----
    parser.add_argument(
        "--out",
        default=None,
        help="Optional output PDBQT file name (relative to dataset), e.g. 'rand-1_out.pdbqt'.",
    )

    # ---- Misc options ----
    parser.add_argument("--cpu", type=int, default=None, help="Number of CPUs to use.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for Vina.")
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=8,
        help="Vina exhaustiveness (default: 8).",
    )
    parser.add_argument(
        "--scoring",
        default="vina",
        choices=["vina"],
        help="Scoring function (only 'vina' is supported).",
    )

    args = parser.parse_args()

    # ---- Prepare dataset ----
    dataset_dir = Path(args.dataset_dir).resolve()
    try:
        dataset_tar_b64 = create_tar_b64(dataset_dir)
    except Exception as e:
        print(f"Error creating tar.gz from dataset directory: {e}")
        raise SystemExit(1)

    # ---- Build request payload ----
    payload = {
        "user_id": args.user_id,
        "workspace_id": args.workspace_id,
        "dataset_tar_b64": dataset_tar_b64,
        "receptor": args.receptor,
        "ligand": args.ligand,
        "scoring": "vina",  # fixed
        "center_x": args.center_x,
        "center_y": args.center_y,
        "center_z": args.center_z,
        "size_x": args.size_x,
        "size_y": args.size_y,
        "size_z": args.size_z,
        "out": args.out,
        "cpu": args.cpu,
        "seed": args.seed,
        "exhaustiveness": args.exhaustiveness,
        "return_tarball": bool(args.download_tarball),
    }

    # Remove keys with value None so JSON is clean
    payload = {k: v for k, v in payload.items() if v is not None}

    url = f"http://{args.host}:{args.port}/v1/vina"

    print("Sending request to:", url)
    print("Payload fields (excluding dataset_tar_b64 content):")
    payload_preview = {
        k: ("<base64-tar.gz>" if k == "dataset_tar_b64" else v)
        for k, v in payload.items()
    }
    print(json.dumps(payload_preview, indent=2))

    # ---- Send request ----
    try:
        resp = requests.post(url, json=payload)
    except Exception as e:
        print(f"Request failed: {e}")
        raise SystemExit(1)

    print("\nHTTP status:", resp.status_code)

    # ---- Try to parse JSON response ----
    try:
        data = resp.json()
    except Exception:
        print("Response is not valid JSON:")
        print(resp.text)
        raise SystemExit(1)

    # ---- Handle non-200 HTTP errors (HTTPException, etc.) ----
    if resp.status_code != 200:
        print("\nServer returned an error response.")
        # FastAPI HTTPException format: {"detail": ...}
        if isinstance(data, dict) and "detail" in data:
            detail = data["detail"]
            if isinstance(detail, dict):
                print("\nError detail JSON:")
                print(json.dumps(detail, indent=2))

                # If the server included stdout/stderr inside detail, show them
                if "stdout" in detail and detail["stdout"]:
                    print("\n=== Vina stdout (from error) ===")
                    print(detail["stdout"])
                if "stderr" in detail and detail["stderr"]:
                    print("\n=== Vina stderr (from error) ===")
                    print(detail["stderr"])
            else:
                # detail is a simple string
                print("\nError detail:")
                print(detail)
        else:
            print("\nFull error JSON:")
            print(json.dumps(data, indent=2))

        raise SystemExit(1)

    # ---- HTTP 200: success path ----
    print("\nResponse JSON:")
    print(json.dumps(data, indent=2))

    # Convenience: highlight stdout/stderr if present
    if "stdout" in data and data["stdout"]:
        print("\n=== Vina stdout ===")
        print(data["stdout"])
    if "stderr" in data and data["stderr"]:
        print("\n=== Vina stderr ===")
        print(data["stderr"])

    # Convenience: save result tarball if requested and returned
    if args.download_tarball and "result_tar_b64" in data and data["result_tar_b64"]:
        try:
            raw = base64.b64decode(data["result_tar_b64"])
            out_path = Path(args.tarball_output).resolve()
            with out_path.open("wb") as f:
                f.write(raw)
            size = len(raw)
            print(f"\nSaved workspace tarball to: {out_path} ({size} bytes)")
        except Exception as e:
            print(f"\nFailed to decode/save result tarball: {e}")


if __name__ == "__main__":
    main()

