#*************************************************************************************
#                           The MIT License
#
#   Intel OpenOmics - dv2vcf: standalone DeepVariant 1.9.0 variant-calling stage
#   Copyright (C) 2023  Intel Corporation.
#
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the
#   "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to
#   permit persons to whom the Software is furnished to do so, subject to
#   the following conditions:
#
#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#*****************************************************************************************/
#
# Takes a single sorted BAM (e.g. the output of pipelines/fq2sortedbam) and a
# reference FASTA, and runs DeepVariant 1.9.0 on it. DeepVariant shards
# make_examples/call_variants internally via --num_shards, so no MPI/rank
# splitting is needed here -- this is a second, independent container stage
# (not nested inside the fq2sortedbam container).

import argparse
import os
import subprocess
import sys
import time

DEEPVARIANT = "/opt/deepvariant/bin/run_deepvariant"
SAMTOOLS = "samtools"

MODEL_TYPES = ["WGS", "WES", "PACBIO", "ONT_R104", "HYBRID_PACBIO_ILLUMINA"]


def run(cmd):
    print("[Info] Running:", cmd, flush=True)
    return subprocess.run(cmd, shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True, help="Reference FASTA (must have a .fai alongside it, or pass --dindex to create one).")
    parser.add_argument("--bam", required=True, help="Input sorted BAM (e.g. output of pipelines/fq2sortedbam).")
    parser.add_argument("--output", required=True, help="Output VCF path.")
    parser.add_argument("--model_type", default="WGS", choices=MODEL_TYPES, help="DeepVariant --model_type.")
    parser.add_argument("--num_shards", default=str(os.cpu_count() or 1), help="DeepVariant --num_shards (defaults to all visible CPUs).")
    parser.add_argument("--dindex", action="store_true", help="Create the reference .fai index if missing.")
    parser.add_argument("--intermediate_results_dir", default="", help="Optional dir for DeepVariant's intermediate results (defaults next to --output).")
    args = parser.parse_args()

    if args.dindex and not os.path.exists(args.ref + ".fai"):
        a = run(f"{SAMTOOLS} faidx {args.ref}")
        assert a.returncode == 0, "samtools faidx failed"
    assert os.path.exists(args.ref + ".fai"), f"Missing {args.ref}.fai -- rerun with --dindex or create it yourself."

    if not os.path.exists(args.bam + ".bai") and not os.path.exists(os.path.splitext(args.bam)[0] + ".bai"):
        a = run(f"{SAMTOOLS} index -@ {args.num_shards} {args.bam}")
        assert a.returncode == 0, "samtools index failed"

    outdir = os.path.dirname(os.path.abspath(args.output)) or "."
    os.makedirs(outdir, exist_ok=True)
    intermediate_dir = args.intermediate_results_dir or os.path.join(outdir, "intermediate_results_dir")
    os.makedirs(intermediate_dir, exist_ok=True)

    tic = time.time()
    cmd = (
        f"{DEEPVARIANT} --model_type={args.model_type} --ref={args.ref} --reads={args.bam} "
        f"--output_vcf={args.output} --intermediate_results_dir={intermediate_dir} "
        f"--num_shards={args.num_shards} --dry_run=false"
    )
    a = run(cmd)
    assert a.returncode == 0, "[Info] DeepVariant execution failed."
    print(f"[Info] DeepVariant runtime: {time.time() - tic:.2f} sec")


if __name__ == "__main__":
    main()
