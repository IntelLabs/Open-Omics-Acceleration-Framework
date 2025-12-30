#!/usr/bin/env python

import argparse
import base64
import json
from pathlib import Path

import requests


def build_flags_from_args(args) -> dict:
    """
    Build the 'flags' dict for the RELION microservice
    from CLI arguments.
    """
    flags = {
        "i": args.i,
        "o": args.o,
        "ref": args.ref,
        "K": args.K,
        "auto_refine": args.auto_refine,
        "dont_combine_weights_via_disc": args.dont_combine_weights_via_disc,
        "ctf": args.ctf,
        "zero_mask": args.zero_mask,
        "norm": args.norm,
        "scale": args.scale,
        "firstiter_cc": args.firstiter_cc,
        "flatten_solvent": args.flatten_solvent,
        "split_random_halves": args.split_random_halves,
        "cpu": True,  
        "oversampling": args.oversampling,
        "random_seed": args.random_seed,
        "preread_images": args.preread_images,
        "center_classes": args.center_classes,
        "blush": args.blush,
        "auto_ignore_angles": args.auto_ignore_angles,
        "auto_resol_angles": args.auto_resol_angles,
        "tau2_fudge": args.tau2_fudge,
        "particle_diameter": args.particle_diameter,
        "psi_step": args.psi_step,
        "offset_range": args.offset_range,
        "offset_step": args.offset_step,
        "pad": args.pad,
        "pool": args.pool,
        "j": args.j,
        "iter": args.iter,
        "ini_high": args.ini_high,
        "healpix_order": args.healpix_order,
        "auto_local_healpix_order": args.auto_local_healpix_order,
        "low_resol_join_halves": args.low_resol_join_halves,
        "sym": args.sym,
        "extra_args": args.extra_args if args.extra_args else None,
    }
    if flags.get("ref") is None:
        flags.pop("ref", None)
    if flags.get("extra_args") is None:
        flags.pop("extra_args", None)
    return flags


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Client for OPEA RELION microservice (classification + auto-refine)."
    )

    # ---- Microservice location ----
    parser.add_argument("--host", default="127.0.0.1", help="Microservice host")
    parser.add_argument("--port", type=int, default=9000, help="Microservice port")

    # ---- Identity ----
    parser.add_argument("--user-id", required=True, help="Logical user id")
    parser.add_argument(
        "--workspace-id",
        default=None,
        help="Logical workspace/job id (optional, used as label)",
    )

    # ---- Dataset S3 URL ----
    parser.add_argument(
        "--dataset-s3-url",
        required=True,
        help="S3 URL to .tar.gz RELION dataset (e.g. s3://bucket/relion_benchmark.tar.gz)",
    )
    parser.add_argument(
            "--output-s3-prefix",
            default=None,
            help=(
                "Optional S3 prefix to upload the parent folder of flags.o.\n"
                "Example: s3://my-bucket/relion_outputs/job_2d_001"
                ),
            )


    # ---- Mode ----
    parser.add_argument(
        "--mode",
        choices=["2d", "3d", "autorefine"],
        default="2d",
        help="RELION mode: 2d, 3d, or autorefine.",
    )

    # ---- MPI ----
    parser.add_argument(
        "--num-mpi",
        type=int,
        default=2,
        help="mpirun -n value (number of MPI ranks).",
    )

    # ---- relion_refine_mpi flags ----
    parser.add_argument(
        "--i",
        required=True,
        help="flags.i: input star file (relative to project root), e.g. Particles/shiny_2sets.star",
    )
    parser.add_argument(
        "--o",
        required=True,
        help="flags.o: output root (relative), e.g. 2D/2D or 3D/3D or Refine3D/job019/run",
    )
    parser.add_argument(
        "--ref",
        default=None,
        help="flags.ref: reference path (relative), e.g. emd_2660.map:mrc or Class3D/job016/run_it025_class004_box256.mrc",
    )
    parser.add_argument(
        "--K",
        type=int,
        default=None,
        help="flags.K: number of classes (required for 2D/3D classification).",
    )

    parser.add_argument(
        "--auto-refine",
        action="store_true",
        help="flags.auto_refine: enable 3D auto-refinement mode.",
    )

    # Boolean flags
    parser.add_argument(
        "--dont-combine-weights-via-disc",
        action="store_true",
        help="flags.dont_combine_weights_via_disc",
    )
    parser.add_argument(
        "--ctf",
        action="store_true",
        help="flags.ctf: perform CTF correction.",
    )
    parser.add_argument(
        "--zero-mask",
        action="store_true",
        help="flags.zero_mask",
    )
    parser.add_argument(
        "--norm",
        action="store_true",
        help="flags.norm",
    )
    parser.add_argument(
        "--scale",
        action="store_true",
        help="flags.scale",
    )
    parser.add_argument(
        "--firstiter-cc",
        action="store_true",
        help="flags.firstiter_cc",
    )
    parser.add_argument(
        "--flatten-solvent",
        action="store_true",
        help="flags.flatten_solvent",
    )
    parser.add_argument(
        "--split-random-halves",
        action="store_true",
        help="flags.split_random_halves",
    )
    parser.add_argument(
        "--tau2-fudge",
        type=float,
        default=None,
        help="flags.tau2_fudge",
    )
    parser.add_argument(
        "--particle-diameter",
        type=float,
        default=None,
        help="flags.particle_diameter",
    )
    parser.add_argument(
        "--oversampling",
        type=int,
        default=None,
        help="flags.oversampling",
    )
    parser.add_argument(
        "--psi-step",
        type=float,
        default=None,
        help="flags.psi_step",
    )
    parser.add_argument(
        "--offset-range",
        type=float,
        default=None,
        help="flags.offset_range",
    )
    parser.add_argument(
        "--offset-step",
        type=float,
        default=None,
        help="flags.offset_step",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="flags.random_seed",
    )
    parser.add_argument(
        "--pad",
        type=int,
        default=None,
        help="flags.pad",
    )
    parser.add_argument(
        "--pool",
        type=int,
        required=True,
        help="flags.pool: particle pool size per MPI rank (>=1, REQUIRED).",
    )
    parser.add_argument(
        "--j",
        type=int,
        required=True,
        help="flags.j: OMP threads per MPI rank (>=1, REQUIRED).",
    )
    parser.add_argument(
        "--iter",
        type=int,
        default=None,
        help="flags.iter (2D/3D only, auto-refine has its own schedule).",
    )
    parser.add_argument(
        "--ini-high",
        type=float,
        default=None,
        help="flags.ini_high",
    )
    parser.add_argument(
        "--healpix-order",
        type=int,
        default=None,
        help="flags.healpix_order",
    )
    parser.add_argument(
        "--auto-local-healpix-order",
        type=int,
        default=None,
        help="flags.auto_local_healpix_order",
    )
    parser.add_argument(
        "--low-resol-join-halves",
        type=float,
        default=None,
        help="flags.low_resol_join_halves",
    )
    parser.add_argument(
        "--sym",
        default=None,
        help="flags.sym (symmetry group), e.g. C1 or D2.",
    )

    # Extra args
    parser.add_argument(
        "--extra-args",
        nargs="*",
        default=None,
        help="Optional extra CLI args to append to relion_refine_mpi.",
    )

    # Save metadata
    parser.add_argument(
        "--save-metadata",
        default=None,
        help="If set, save decoded metadata JSON to this path.",
    )
    parser.add_argument(
        "--preread-images",
        action="store_true",
        help="flags.preread_images",
    )
    parser.add_argument(
        "--center-classes",
        action="store_true",
        help="flags.center_classes",
    )
    parser.add_argument(
        "--blush",
        action="store_true",
        help="flags.blush",
    )
    parser.add_argument(
        "--auto-ignore-angles",
        action="store_true",
        help="flags.auto_ignore_angles",
    )
    parser.add_argument(
        "--auto-resol-angles",
        action="store_true",
        help="flags.auto_resol_angles",
    )

    args = parser.parse_args()

    url = f"http://{args.host}:{args.port}/v1/relion"

    flags = build_flags_from_args(args)

    payload = {
        "user_id": args.user_id,
        "workspace_id": args.workspace_id,
        "dataset_s3_url": args.dataset_s3_url,
        "mode": args.mode,
        "num_mpi": args.num_mpi,
        "output_s3_prefix": args.output_s3_prefix,
        "flags": flags,
    }

    print(f"POST {url}")
    print("Request JSON:")
    print(json.dumps(payload, indent=2))

    resp = requests.post(url, json=payload)
    print(f"\nStatus: {resp.status_code}")

    if resp.status_code != 200:
        print("\nRaw JSON response:")
        try:
            print(json.dumps(resp.json(), indent=2))
        except Exception:
            print(resp.text)
        print("Server returned error, stopping.")
        return

    data = resp.json()
    print("\nRaw JSON response:")
    print(json.dumps(data, indent=2))

    metadata_b64 = data.get("metadata_b64")
    if metadata_b64:
        try:
            decoded = base64.b64decode(metadata_b64.encode("ascii")).decode("utf-8")
            print("\nDecoded metadata.json:")
            print(decoded)

            if args.save_metadata:
                out_path = Path(args.save_metadata)
                out_path.write_text(decoded)
                print(f"\nSaved metadata.json to: {out_path}")
        except Exception as e:
            print(f"\nFailed to decode metadata_b64: {e}")


if __name__ == "__main__":
    main()

