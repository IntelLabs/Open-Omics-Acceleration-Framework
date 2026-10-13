# fq2vcf: OpenOmics Deepvariant based Variant Calling Pipeline  
### Overview:  
OpenOmics' fq2vcf is a highly optimized, distributed, deep learning-based short-read germline variant calling pipeline for x86 CPUs. It is built as a chain of two independent stages:

1. **`fq2sortedbam`** -- Sequence alignment and sorting stage. Takes raw FASTQ reads and a reference genome and produces a coordinate-sorted BAM file. Depending on the data type, it dispatches to one of the following aligners:
   - **bwa-mem2 v2.2.1** -- highly optimized DNA short-read alignment
   - **mm2-fast** (accelerated minimap2, based on minimap2 v2.24-r1122) -- long-read alignment (PacBio/ONT)
   - **samtools v1.16.1** -- SAM/BAM sorting
2. **`bams2vcf`** (containerized as `dv2vcf`) -- Variant calling stage. Takes the sorted BAM and reference genome and calls variants using an optimized build of **DeepVariant v1.9.0**, producing a compressed, indexed VCF file.

The following figure illustrates the pipeline:

<p align="center">
<img src="https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/blob/main/images/deepvariant-fq2vcf.jpg"/a></br>
</p>

## Instructions to build and run the pipeline
### 1. Download the code
```bash
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
cd Open-Omics-Acceleration-Framework
git checkout fq2vcf-dv1.9.0
git submodule update --init --recursive
```

### 2. Build both images (from the repository root)
```bash
docker build -f pipelines/fq2sortedbam/Dockerfile -t fq2sortedbam:latest .
docker build -f pipelines/deepvariant-based-germline-variant-calling-fq2vcf/Dockerfile_dv2vcf -t dv2vcf:1.9.0 .
```

### 3. Run both stages with the chain script
```bash
pipelines/deepvariant-based-germline-variant-calling-fq2vcf/run_fq2vcf_chain.sh \
  --container-tool docker \
  --refdir <refdir> --ref <reference.fasta> \
  --readsdir <readsdir> --read1 <r1.fastq.gz> --read2 <r2.fastq.gz> \
  --outdir <outdir> --prefix <prefix> --model-type WGS
```
Notes:
- Produces `<outdir>/<prefix>.sorted.bam` (fq2sortedbam stage), then `<outdir>/<prefix>.vcf.gz` + `.tbi` (dv2vcf stage).
- See `run_fq2vcf_chain.sh --help` for all options: `--model-type` for PacBio/ONT/hybrid, `--skip-align`/`--skip-dv` to run just one stage.
- `--refdir` is expected to already contain the bwa-mem2 index; add `--rindex` to build it on the fly if it isn't present yet.
- Favor local/fast disk (e.g. NVMe scratch) over networked/shared storage for input, output, and intermediate files -- both the alignment and DeepVariant stages are I/O-sensitive and slow down noticeably on network-backed filesystems.
- The first run against a given FASTQ pair will be slower because `fq2sortedbam` builds a pragzip index over the compressed reads from cold disk; subsequent runs on the same machine are faster once the input file is cached in the memory.
- `TF_NUM_INTRAOP_THREADS`/`TF_NUM_INTEROP_THREADS`/`OPENBLAS_NUM_THREADS` cap the thread pool so our concurrent MPI processes don't oversubscribe the node; for example, `MPI processes =  16` x `--intraop-threads 16` = 256 threads, matching a 256-core node. Tune these via `run_fq2vcf_chain.sh --intraop-threads`/`--interop-threads`/`--openblas-threads` (or the equivalent `run_dv2vcf.py --intraop_threads`/`--interop_threads`/`--openblas_threads` flags if invoking stage 2 directly).

#### Verified results (HG001, full 30x NovaSeq WGS, GRCh38, on a 512-CPU / 6-NUMA node, local NVMe storage)
| Stage | Time |
|---|---|
| Alignment (bwa-mem2, 16 ranks x 28 threads, auto-detected) | 534s |
| DeepVariant 1.9.0 (region-sharded, 16 bins x 32 shards) | 708s |
| **Total FASTQ -> VCF** | **~20.7 min** |

### 4. Validate accuracy with hap.py (optional)
Compare the output VCF against a GIAB truth set, e.g. for HG003/NA24149:
```bash
# Download the GIAB v4.2.1 truth set (once)
mkdir -p truth && cd truth
wget https://ftp.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
cd ..

docker pull jmcdani20/hap.py:v0.3.12
mkdir -p happy_results
docker run --rm \
  -v $(pwd)/truth:/benchmark \
  -v <outdir>:/output \
  -v <refdir>:/reference \
  -v $(pwd)/happy_results:/happy \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/<prefix>.vcf.gz \
  -f /benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r /reference/<reference.fasta> \
  -o /happy/happy.output \
  --engine=vcfeval --pass-only --threads $(nproc)
```
Swap in the matching GIAB directory/sample name (`HG001`/`HG002`/`HG003`/...) for other reference samples. Results (precision/recall/F1 per variant type) land in `happy_results/happy.output.*`.

#### Verified results (HG003, full 35x NovaSeq WGS, GRCh38 no-ALT reference, hap.py vs GIAB v4.2.1 truth)
| Type | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | Recall | Precision | Frac_NA | F1_Score |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| INDEL | ALL | 504501 | 501529 | 2972 | 937912 | 1185 | 414081 | 788 | 217 | 0.994109 | 0.997738 | 0.441492 | 0.995920 |
| INDEL | PASS | 504501 | 501529 | 2972 | 937912 | 1185 | 414081 | 788 | 217 | 0.994109 | 0.997738 | 0.441492 | 0.995920 |
| SNP | ALL | 3327496 | 3305923 | 21573 | 3816076 | 4457 | 504028 | 1437 | 471 | 0.993517 | 0.998654 | 0.132080 | 0.996079 |
| SNP | PASS | 3327496 | 3305923 | 21573 | 3816076 | 4457 | 504028 | 1437 | 471 | 0.993517 | 0.998654 | 0.132080 | 0.996079 |

These numbers match DeepVariant 1.9.0's own officially published full-genome HG003 benchmark (INDEL F1 0.995845, SNP F1 0.996133 -- see [metrics.md](https://github.com/google/deepvariant/blob/r1.9/docs/metrics.md)) to within ~0.0002 across Recall/Precision/F1 for both variant types.

# Results

Our latest results are published in this [blog](https://community.intel.com/t5/Blogs/Tech-Innovation/Artificial-Intelligence-AI/Intel-Xeon-is-all-you-need-for-AI-inference-Performance/post/1506083).


# Instructions to run the pipeline on an AWS ec2 instance (Single Node)
The above instructions run seamlessly on a standalone AWS ec2 instance. To run the steps, create an ec2 instance with Ubuntu-22.04 having at least 60GB of memory and 500GB of disk. The input reference sequence and the paired-ended read datasets must be downloaded and stored on the disk.

