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
# reference FASTA, and runs DeepVariant 1.9.0 on it.
#
# DeepVariant's own --num_shards only parallelizes the make_examples stage;
# call_variants always runs as a single TensorFlow process over the whole
# input, which becomes the bottleneck on a full 30x genome (confirmed by
# benchmarking: ~19 min for call_variants alone vs ~8 min for make_examples
# with 512 shards). To fix this we replicate the region-sharding strategy
# from pipelines/deepvariant-based-germline-variant-calling-fq2vcf/test_pipeline_final.py:
# split the genome into --bins regions (by cumulative sequence length, same
# algorithm as that script's calculate_bins()) and launch one independent
# run_deepvariant process per region (each with its own small --regions
# restriction, so its call_variants only has a fraction of the examples),
# all running concurrently. The per-bin VCFs are then concatenated in
# genome order with bcftools concat.
#
# This is a second, independent container stage (not nested inside the
# fq2sortedbam container) -- the DeepVariant binaries already live in this
# image, so "sub-containers" per bin are not needed, just parallel
# subprocesses.

import argparse
import os
import subprocess
import sys
import time

DEEPVARIANT = "/opt/deepvariant/bin/run_deepvariant"
SAMTOOLS = "samtools"
BCFTOOLS = "bcftools"

MODEL_TYPES = ["WGS", "WES", "PACBIO", "ONT_R104", "HYBRID_PACBIO_ILLUMINA"]


def run(cmd, **kwargs):
    print("[Info] Running:", cmd, flush=True)
    return subprocess.run(cmd, shell=True, **kwargs)


def calculate_bin_regions(fai_path, nbins, binrounding=1000):
    """Split the reference's sequences into nbins contiguous regions of
    roughly equal total length, same algorithm as test_pipeline_final.py's
    calculate_bins(), adapted from the MPI rank/bin form to a flat list.
    Returns a list of nbins region strings suitable for DeepVariant --regions
    (space-separated "seq:start-end" tokens; a bin can span multiple
    sequences when one finishes partway through it).
    """
    seq_names = []
    seq_lens = []
    with open(fai_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            seq_names.append(parts[0])
            seq_lens.append(int(parts[1]))

    seq_start = {}
    cumlen = 0
    for name, length in zip(seq_names, seq_lens):
        seq_start[name] = cumlen
        cumlen += length

    bin_regions = []
    seq_i = 0
    start = 0
    last = len(seq_names) - 1
    for bin_i in range(nbins):
        parts = []
        end = (bin_i + 1) * cumlen // nbins
        while seq_i < last and seq_start[seq_names[seq_i]] + seq_lens[seq_i] < end - binrounding:
            name = seq_names[seq_i]
            parts.append(f"{name}:{max(0, start - seq_start[name])}-{seq_lens[seq_i]}")
            seq_i += 1
        name = seq_names[seq_i]
        seq_end = seq_start[name] + seq_lens[seq_i]
        if bin_i == nbins - 1:
            end = cumlen  # last bin always runs to the end of the genome
        elif abs(end - seq_end) <= binrounding:
            end = seq_end
        else:
            end = (end - seq_start[name] + binrounding // 2) // binrounding * binrounding + seq_start[name]
        parts.append(f"{name}:{max(0, start - seq_start[name])}-{min(end - seq_start[name], seq_lens[seq_i])}")
        bin_regions.append(" ".join(parts))
        start = end
    return bin_regions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True, help="Reference FASTA (must have a .fai alongside it, or pass --dindex to create one).")
    parser.add_argument("--bam", required=True, help="Input sorted BAM (e.g. output of pipelines/fq2sortedbam).")
    parser.add_argument("--output", required=True, help="Output VCF path.")
    parser.add_argument("--model_type", default="WGS", choices=MODEL_TYPES, help="DeepVariant --model_type.")
    parser.add_argument("--num_shards", default=str(os.cpu_count() or 1), help="Total make_examples shards, split evenly across --bins (defaults to all visible CPUs).")
    parser.add_argument("--bins", type=int, default=16, help="Number of parallel region-sharded DeepVariant invocations (set to 1 to disable region-sharding).")
    parser.add_argument("--intraop_threads", default="16", help="TF_NUM_INTRAOP_THREADS per bin (DeepVariant's call_variants thread count).")
    parser.add_argument("--interop_threads", default="1", help="TF_NUM_INTEROP_THREADS per bin.")
    parser.add_argument("--openblas_threads", default="1", help="OPENBLAS_NUM_THREADS per bin.")
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

    nbins = max(1, args.bins)
    shards_per_bin = max(1, int(args.num_shards) // nbins)
    env = dict(os.environ)
    env["OPENBLAS_NUM_THREADS"] = args.openblas_threads
    env["TF_NUM_INTRAOP_THREADS"] = args.intraop_threads
    env["TF_NUM_INTEROP_THREADS"] = args.interop_threads

    tic = time.time()

    if nbins == 1:
        cmd = (
            f"{DEEPVARIANT} --model_type={args.model_type} --ref={args.ref} --reads={args.bam} "
            f"--output_vcf={args.output} --intermediate_results_dir={intermediate_dir} "
            f"--num_shards={args.num_shards} --dry_run=false"
        )
        a = run(cmd)
        assert a.returncode == 0, "[Info] DeepVariant execution failed."
        print(f"[Info] DeepVariant runtime: {time.time() - tic:.2f} sec")
        return

    print(f"[Info] Region-sharding DeepVariant into {nbins} bins ({shards_per_bin} make_examples shards/bin, "
          f"TF_NUM_INTRAOP_THREADS={args.intraop_threads})", flush=True)
    bin_regions = calculate_bin_regions(args.ref + ".fai", nbins)

    bin_vcfs = []
    procs = []
    for i, region in enumerate(bin_regions):
        binstr = f"{i:05d}"
        bin_outdir = os.path.join(outdir, binstr)
        os.makedirs(bin_outdir, exist_ok=True)
        bin_vcf = os.path.join(bin_outdir, "output.vcf.gz")
        bin_intermediate = os.path.join(intermediate_dir, binstr)
        os.makedirs(bin_intermediate, exist_ok=True)
        bin_vcfs.append(bin_vcf)
        logfile = os.path.join(bin_outdir, "log.txt")
        # Some bins (e.g. the one covering thousands of small decoy/HLA contigs) produce a
        # --regions string too long for the OS argv limit when inlined. Write a BED file
        # instead (our region tokens are already 0-based half-open, i.e. BED-compatible).
        bed_path = os.path.join(bin_outdir, "regions.bed")
        with open(bed_path, "w") as bed:
            for token in region.split():
                name, coords = token.rsplit(":", 1)
                start, end = coords.split("-")
                bed.write(f"{name}\t{start}\t{end}\n")
        cmd = (
            f"{DEEPVARIANT} --model_type={args.model_type} --ref={args.ref} --reads={args.bam} "
            f"--output_vcf={bin_vcf} --intermediate_results_dir={bin_intermediate} "
            f"--num_shards={shards_per_bin} --dry_run=false --regions {bed_path} "
            f"> {logfile} 2>&1"
        )
        print("[Info] Launching bin", binstr, "regions:", bed_path, flush=True)
        procs.append(subprocess.Popen(cmd, shell=True, env=env))

    failed = [i for i, p in enumerate(procs) if p.wait() != 0]
    assert not failed, f"[Info] DeepVariant failed for bin(s): {failed}"
    print(f"[Info] All {nbins} bins finished in {time.time() - tic:.2f} sec, merging VCFs...", flush=True)

    vcf_list = " ".join(bin_vcfs)
    a = run(f"{BCFTOOLS} concat {vcf_list} -Oz -o {args.output}")
    assert a.returncode == 0, "bcftools concat failed"
    a = run(f"{BCFTOOLS} index -t {args.output}")
    assert a.returncode == 0, "bcftools index failed"

    print(f"[Info] DeepVariant runtime: {time.time() - tic:.2f} sec")


if __name__ == "__main__":
    main()
