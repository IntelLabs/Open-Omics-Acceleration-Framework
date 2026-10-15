#!/bin/bash
#*************************************************************************************
#                           The MIT License
#
#   Intel OpenOmics - fq2vcf two-container chain: fq2sortedbam -> bams2vcf (DeepVariant 1.9.0)
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
# Runs the two fq2vcf containers back-to-back:
#   1. fq2sortedbam (aligner image, e.g. bwa-mem2/mm2-fast/STAR/bwa-meth) -> sorted BAM
#   2. bams2vcf (DeepVariant 1.9.0 image)                                -> VCF
#
# Usage:
#   ./run_fq2vcf_chain.sh --refdir DIR --ref FILE --readsdir DIR --read1 FILE --read2 FILE \
#       --outdir DIR --prefix NAME [options]
#
# Options (with defaults):
#   --refdir DIR            Directory containing the reference FASTA (and bwa-mem2 index). Required, must be an absolute path.
#   --ref FILE              Reference FASTA filename (relative to --refdir). Required.
#   --readsdir DIR          Directory containing the input FASTQ(.gz) files. Required, must be an absolute path.
#   --read1 FILE            Read 1 FASTQ(.gz) filename (relative to --readsdir). Required.
#   --read2 FILE            Read 2 FASTQ(.gz) filename (relative to --readsdir). Omit for single-end.
#   --outdir DIR            Output directory for both the BAM and the VCF. Required, must be an absolute path.
#   --prefix NAME           Output file prefix (default: fq2vcf_out).
#   --aligner-image NAME    fq2sortedbam image tag (default: fq2sortedbam:latest).
#   --bams2vcf-image NAME   bams2vcf image tag (default: bams2vcf:1.9.0).
#   --container-tool TOOL   docker or podman (default: docker).
#   --read-type TYPE        fq2sortedbam --read_type: short/long/rna/meth (default: short). short - bwa-mem2
#                            for short reads; long - mm2-fast for PacBio/ONT long reads; rna - STAR; meth -
#                            bwa-meth. Use long with --model-type PACBIO/ONT_R104 for long-read data.
#   --model-type TYPE       DeepVariant --model_type: WGS/WES/PACBIO/ONT_R104/HYBRID_PACBIO_ILLUMINA (default: WGS).
#   --num-shards N          Total DeepVariant make_examples shards, split across --bins (default: nproc).
#   --bins N                Number of parallel region-sharded DeepVariant invocations (default: the MPI rank
#                            count auto-detected by the fq2sortedbam stage, read back from
#                            <outdir>/num_ranks.txt; falls back to 16 if that file isn't present, e.g. with
#                            --skip-align; use 1 to disable).
#   --skip-align            Skip stage 1 and go straight to DeepVariant on an existing
#                            <outdir>/<prefix>.sorted.bam.
#   --skip-dv                Run only stage 1 (alignment), skip DeepVariant.
#   --rindex                Build the bwa-mem2 index for --ref on the fly if it's not
#                            already present in --refdir (passed through to
#                            run_fq2sortedbam.py's --rindex).
#   --intraop-threads N      TF_NUM_INTRAOP_THREADS per bin (default: 16).
#   --interop-threads N      TF_NUM_INTEROP_THREADS per bin (default: 1).
#   --openblas-threads N     OPENBLAS_NUM_THREADS per bin (default: 1).
#   --no-thread-caps         Don't set TF_NUM_INTRAOP_THREADS/TF_NUM_INTEROP_THREADS/OPENBLAS_NUM_THREADS
#                            at all; let TF/OpenBLAS use their own defaults (overrides the three flags above).
#
# Example:
#   ./run_fq2vcf_chain.sh --container-tool podman \
#       --refdir /data/ref --ref Homo_sapiens_assembly38.fasta \
#       --readsdir /data/reads --read1 R1.fastq.gz --read2 R2.fastq.gz \
#       --outdir /data/out --prefix HG001 --model-type WGS

set -euo pipefail

CONTAINER_TOOL="docker"
ALIGNER_IMAGE="fq2sortedbam:latest"
BAMS2VCF_IMAGE="bams2vcf:1.9.0"
PREFIX="fq2vcf_out"
READ_TYPE="short"
MODEL_TYPE="WGS"
NUM_SHARDS="$(nproc 2>/dev/null || echo 16)"
BINS=""
INTRAOP_THREADS="16"
INTEROP_THREADS="1"
OPENBLAS_THREADS="1"
READ2=""
SKIP_ALIGN=0
SKIP_DV=0
RINDEX=0
NO_THREAD_CAPS=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --refdir) REFDIR="$2"; shift 2 ;;
        --ref) REF="$2"; shift 2 ;;
        --readsdir) READSDIR="$2"; shift 2 ;;
        --read1) READ1="$2"; shift 2 ;;
        --read2) READ2="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        --aligner-image) ALIGNER_IMAGE="$2"; shift 2 ;;
        --bams2vcf-image) BAMS2VCF_IMAGE="$2"; shift 2 ;;
        --container-tool) CONTAINER_TOOL="$2"; shift 2 ;;
        --read-type) READ_TYPE="$2"; shift 2 ;;
        --model-type) MODEL_TYPE="$2"; shift 2 ;;
        --num-shards) NUM_SHARDS="$2"; shift 2 ;;
        --bins) BINS="$2"; shift 2 ;;
        --intraop-threads) INTRAOP_THREADS="$2"; shift 2 ;;
        --interop-threads) INTEROP_THREADS="$2"; shift 2 ;;
        --openblas-threads) OPENBLAS_THREADS="$2"; shift 2 ;;
        --no-thread-caps) NO_THREAD_CAPS=1; shift ;;
        --skip-align) SKIP_ALIGN=1; shift ;;
        --skip-dv) SKIP_DV=1; shift ;;
        --rindex) RINDEX=1; shift ;;
        -h|--help) grep '^#' "$0" | sed 's/^#//'; exit 0 ;;
        *) echo "[Error] Unknown argument: $1"; exit 1 ;;
    esac
done

: "${REFDIR:?--refdir is required}"
: "${REF:?--ref is required}"
: "${OUTDIR:?--outdir is required}"
if [[ "$SKIP_ALIGN" -eq 0 ]]; then
    : "${READSDIR:?--readsdir is required}"
    : "${READ1:?--read1 is required}"
fi

mkdir -p "$OUTDIR"

if [[ "$SKIP_ALIGN" -eq 0 ]]; then
    echo "[Info] Stage 1/2: fq2sortedbam ($ALIGNER_IMAGE) -> ${OUTDIR}/${PREFIX}.sorted.bam"
    READS_ARGS=(--reads "/input/${READ1}")
    if [[ -n "$READ2" ]]; then
        READS_ARGS=(--reads "/input/${READ1}" "/input/${READ2}")
    fi
    RINDEX_ARGS=()
    if [[ "$RINDEX" -eq 1 ]]; then
        RINDEX_ARGS=(--rindex)
    fi
    "$CONTAINER_TOOL" run --rm \
        -v "${REFDIR}:/refdir" \
        -v "${READSDIR}:/input" \
        -v "${OUTDIR}:/out" \
        "$ALIGNER_IMAGE" \
        python run_fq2sortedbam.py --ref "/refdir/${REF}" "${READS_ARGS[@]}" --output "/out/${PREFIX}" \
            --read_type "$READ_TYPE" "${RINDEX_ARGS[@]}"
else
    echo "[Info] Stage 1/2: skipped (--skip-align), using existing ${OUTDIR}/${PREFIX}.sorted.bam"
fi

if [[ "$SKIP_DV" -eq 0 ]]; then
    if [[ -z "$BINS" ]]; then
        if [[ -f "${OUTDIR}/num_ranks.txt" ]]; then
            BINS="$(cat "${OUTDIR}/num_ranks.txt")"
            echo "[Info] --bins not set; using ${BINS} (MPI ranks auto-detected by the fq2sortedbam stage)"
        else
            BINS=16
            echo "[Info] --bins not set and ${OUTDIR}/num_ranks.txt not found (e.g. --skip-align); defaulting to ${BINS}"
        fi
    fi
    echo "[Info] Stage 2/2: bams2vcf ($BAMS2VCF_IMAGE, DeepVariant 1.9.0) -> ${OUTDIR}/${PREFIX}.vcf.gz"
    # region-sharding fans out to --bins x --num_shards concurrent make_examples
    # processes (each opening several TF .so/library fds), which can exceed the
    # container's default nofile ulimit (independent of the host shell's ulimit).
    "$CONTAINER_TOOL" run --rm \
        --ulimit nofile=131072:131072 \
        -v "${REFDIR}:/refdir" \
        -v "${OUTDIR}:/bamdir" \
        -v "${OUTDIR}:/output" \
        "$BAMS2VCF_IMAGE" \
        python run_bams2vcf.py --ref "/refdir/${REF}" --bam "/bamdir/${PREFIX}.sorted.bam" \
            --output "/output/${PREFIX}.vcf.gz" --model_type "$MODEL_TYPE" --num_shards "$NUM_SHARDS" --bins "$BINS" \
            --intraop_threads "$INTRAOP_THREADS" --interop_threads "$INTEROP_THREADS" --openblas_threads "$OPENBLAS_THREADS" \
            $([[ "$NO_THREAD_CAPS" -eq 1 ]] && echo --no_thread_caps) --dindex
else
    echo "[Info] Stage 2/2: skipped (--skip-dv)"
fi

echo "[Info] Done. Output VCF: ${OUTDIR}/${PREFIX}.vcf.gz"
