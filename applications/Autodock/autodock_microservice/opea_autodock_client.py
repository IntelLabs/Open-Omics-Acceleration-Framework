#!/usr/bin/env python3
import argparse
import base64
import json
import tarfile
from io import BytesIO
from pathlib import Path
from typing import Optional

import requests


def create_tar_b64(dataset_dir: Path) -> str:
    """
    Create an in-memory tar.gz from dataset_dir and return as base64-encoded string.
    """
    if not dataset_dir.is_dir():
        raise ValueError(f"Dataset directory does not exist or is not a directory: {dataset_dir}")

    buf = BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        tf.add(dataset_dir, arcname=dataset_dir.name)

    buf.seek(0)
    data = buf.read()
    return base64.b64encode(data).decode("utf-8")


def decode_metadata_b64(metadata_b64: str) -> dict:
    raw = base64.b64decode(metadata_b64.encode("utf-8"))
    return json.loads(raw.decode("utf-8"))


def save_results_tar_b64(results_tar_b64: str, output_path: Path) -> None:
    raw = base64.b64decode(results_tar_b64.encode("utf-8"))
    output_path.write_bytes(raw)


def _print_http_error(res: requests.Response) -> None:
    """
    Pretty-print HTTPException errors using the standard API error schema:
      {
        "detail": {
          "error_type": "...",
          "message": "...",
          "details": {...} | null
        }
      }
    Falls back to raw text if schema is not present.
    """
    print("Server returned HTTP error:", res.status_code)
    try:
        body = res.json()
    except Exception:
        print(res.text)
        return

    detail = body.get("detail")
    if isinstance(detail, dict):
        err_type = detail.get("error_type")
        msg = detail.get("message")
        details = detail.get("details")

        if err_type:
            print(f"  error_type: {err_type}")
        if msg:
            print(f"  message   : {msg}")
        if details:
            print("  details   :")
            print(json.dumps(details, indent=2))
    else:
        # Unexpected structure – just dump JSON
        print(json.dumps(body, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Client for OPEA AutoDock-GPU microservice."
    )

    # Microservice location (default port 9012)
    parser.add_argument("--host", default="127.0.0.1", help="Microservice host")
    parser.add_argument("--port", type=int, default=9012, help="Microservice port")

    # User / workspace
    parser.add_argument("--user-id", required=True, help="Logical user identifier")
    parser.add_argument(
        "--workspace-id",
        default=None,
        help="Optional workspace/job ID (for grouping jobs).",
    )

    # Dataset
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        required=True,
        help="Directory containing ligand PDBQT, maps FLD, etc.",
    )

    # AutoDock flags
    parser.add_argument(
        "--ffile",
        required=True,
        help="Grid descriptor FLD file name inside dataset-dir (e.g. protein.maps.fld)",
    )
    parser.add_argument(
        "--lfile",
        required=True,
        help="Ligand PDBQT file name inside dataset-dir (e.g. rand-0.pdbqt)",
    )
    parser.add_argument(
        "--nrun",
        type=int,
        default=20,
        help="Number of LGA runs (default: 20).",
    )
    parser.add_argument(
        "--lsmet",
        default="ad",
        choices=["ad", "sw"],
        help="Local search method (ad or sw).",
    )
    parser.add_argument(
        "--seed",
        default=None,
        help="Seeds for RNG (e.g. '11,23' or '11,23,37').",
    )
    parser.add_argument(
        "--nev",
        type=int,
        default=2_500_000,
        help="Max score evaluations per LGA run.",
    )
    parser.add_argument(
        "--resnam",
        required=True,
        help="Base name for docking outputs (for --resnam / -N).",
    )
    parser.add_argument(
        "--results-out",
        type=Path,
        default=Path("autodock_results.tar.gz"),
        help="Where to save results tar.gz (resnam.dlg/xml + metadata.json).",
    )

    args = parser.parse_args()

    url = f"http://{args.host}:{args.port}/v1/autodock"

    dataset_tar_b64 = create_tar_b64(args.dataset_dir)

    payload = {
        "user_id": args.user_id,
        "workspace_id": args.workspace_id,
        "dataset_tar_b64": dataset_tar_b64,
        # Relative paths from extracted root
        "ffile": args.ffile,
        "lfile": args.lfile,
        "nrun": args.nrun,
        "lsmet": args.lsmet,
        "seed": args.seed,
        "nev": args.nev,
        "resnam": args.resnam,
    }

    print(f"POST {url}")
    res = requests.post(url, json=payload, timeout=None)

    print(f"Status: {res.status_code}")
    if res.status_code != 200:
        _print_http_error(res)
        return

    res_json = res.json()
    print("Raw JSON response:")
    print(json.dumps(res_json, indent=2))

    status = res_json.get("status")
    if status == "failure":
        # Job-level failure (AutoDock non-zero exit)
        error_str = res_json.get("error", "")
        print("\n=== AutoDock FAILURE ===")
        print(error_str)
        return

    if status != "success":
        print("\nUnexpected status field:", status)
        return

    print("\n=== AutoDock SUCCESS ===")

    # Decode metadata.json
    metadata_b64 = res_json.get("metadata_b64")
    if metadata_b64:
        metadata = decode_metadata_b64(metadata_b64)
        print("\nDecoded metadata.json:")
        print(json.dumps(metadata, indent=2))
    else:
        print("\nNo metadata_b64 found in response.")

    # Decode tar.gz with resnam.dlg, resnam.xml, metadata.json
    results_tar_b64 = res_json.get("results_tar_b64")
    if results_tar_b64:
        save_results_tar_b64(results_tar_b64, args.results_out)
        print(f"\nSaved results tar.gz to: {args.results_out}")
    else:
        print("\nNo results_tar_b64 found in response.")


if __name__ == "__main__":
    main()

