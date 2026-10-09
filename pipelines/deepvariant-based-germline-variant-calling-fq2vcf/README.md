# fq2vcf: OpenOmics Deepvariant based Variant Calling Pipeline  
### Overview:  
OpenOmics' fq2vcf is a highly optimized, distributed, deep learning-based short-read germline variant calling pipeline for x86 CPUs. 
The pipeline comprises of:   
1. bwa-mem2 (a highly optimized version of bwa-mem) for sequence mapping  
2. SortSAM using samtools  
3. An optimized version of DeepVariant tool for Variant Calling   
The following figure illustrates the pipeline:

## Two-container chain (recommended, verified end-to-end): fq2sortedbam + dv2vcf
This runs the pipeline as two independent, sequential containers instead of merging everything into one image:
1. **`pipelines/fq2sortedbam`** (unchanged, existing pipeline) -- aligns reads with any of bwa-mem2/mm2-fast/STAR/bwa-meth and produces a sorted BAM.
2. **`Dockerfile_dv2vcf`** -- a small image containing only **DeepVariant 1.9.0** (binaries + models copied from `docker.io/google/deepvariant:1.9.0`, on an `ubuntu:22.04` base matching its own environment) plus `samtools`/`bcftools`. Takes any sorted BAM + reference and calls variants.

This avoids the fragility of copying DeepVariant's files into a different base image alongside the aligner build -- the two images stay independent and each matches its own dependencies exactly.

### DeepVariant region-sharding (`run_dv2vcf.py`)
`run_dv2vcf.py` splits the genome into `--bins` regions (same cumulative-length binning algorithm as the bare-metal `test_pipeline_final.py`) and launches one independent `run_deepvariant` subprocess per region, all running concurrently inside the single `dv2vcf` container -- no nested/sibling containers needed since the DeepVariant binaries already live in that image. Per-bin VCFs are merged in genome order with `bcftools concat`. Region lists are written to a BED file per bin rather than inlined on the command line, since the decoy/HLA-heavy bin's region list can exceed the OS `argv` length limit. Default `--bins 16`; pass `--bins 1` to disable region-sharding and fall back to a single whole-genome invocation.

Each bin's `run_deepvariant` subprocess is launched with `TF_NUM_INTRAOP_THREADS=16`, `TF_NUM_INTEROP_THREADS=1`, and `OPENBLAS_NUM_THREADS=1` (matching the bare-metal tuning for a 256-physical-core node: 16 bins x 16 threads = 256). These are configurable via `run_dv2vcf.py --intraop_threads`/`--interop_threads`/`--openblas_threads` if you need to tune for different hardware (e.g. increase `--intraop_threads` and decrease `--bins` on a smaller node, or vice versa on a larger one).

### 1. Build both images (from the repository root)
```bash
# pipelines/fq2sortedbam/Dockerfile currently fails to build as-is: the `python:3.10`
# Docker Hub tag now resolves to a Debian release without gcc-11. Use the local-test
# variant instead, which builds from your checkout and uses the default gcc:
docker build -f pipelines/fq2sortedbam/Dockerfile.localtest -t fq2sortedbam:latest .
docker build -f pipelines/deepvariant-based-germline-variant-calling-fq2vcf/Dockerfile_dv2vcf -t dv2vcf:1.9.0 .
```
Note: `pipelines/fq2sortedbam`'s bundled `bwa-mem2` submodule (`ext/safestringlib`) needs `<stdlib.h>`/`<ctype.h>` included for `abort()`/`toupper()` to compile under modern GCC -- if you hit `implicit declaration of function` errors there, add those two includes near the top of `applications/bwa-mem2/ext/safestringlib/safeclib/safeclib_private.h` (this is a vendored third-party submodule, so the fix doesn't persist across a fresh `git submodule update`; re-apply it if needed):
```bash
sed -i '/#include <stdio.h>/a #include <stdlib.h>\n#include <ctype.h>' \
  applications/bwa-mem2/ext/safestringlib/safeclib/safeclib_private.h
```

### 2. Run both stages with the chain script
```bash
pipelines/deepvariant-based-germline-variant-calling-fq2vcf/run_fq2vcf_chain.sh \
  --container-tool docker \
  --refdir <refdir> --ref <reference.fasta> \
  --readsdir <readsdir> --read1 <r1.fastq.gz> --read2 <r2.fastq.gz> \
  --outdir <outdir> --prefix <prefix> --model-type WGS
```
This runs `fq2sortedbam` to produce `<outdir>/<prefix>.sorted.bam`, then `dv2vcf` to produce `<outdir>/<prefix>.vcf.gz` (+ `.tbi`). See `run_fq2vcf_chain.sh --help` for all options (custom image tags, `podman` support, `--model-type` for PacBio/ONT/hybrid, `--skip-align`/`--skip-dv` to run just one stage).

### 3. Validate accuracy with hap.py (optional)
Compare the output VCF against a GIAB truth set, e.g. for HG001/NA12878:
```bash
# Download the GIAB v4.2.1 truth set (once)
mkdir -p truth && cd truth
wget https://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
cd ..

docker pull jmcdani20/hap.py:v0.3.12
mkdir -p happy_results
docker run --rm \
  -v $(pwd)/truth:/benchmark \
  -v <outdir>:/output \
  -v <refdir>:/reference \
  -v $(pwd)/happy_results:/happy \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/<prefix>.vcf.gz \
  -f /benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
  -r /reference/<reference.fasta> \
  -o /happy/happy.output \
  --engine=vcfeval --pass-only --threads $(nproc)
```
Swap in the matching GIAB directory/sample name (`HG001`/`HG002`/`HG003`/...) for other reference samples. Results (precision/recall/F1 per variant type) land in `happy_results/happy.output.*`.

#### Verified results (HG001, full 30x NovaSeq WGS, GRCh38, on a 512-CPU / 6-NUMA node, local NVMe storage)
| Stage | Time |
|---|---|
| Alignment (bwa-mem2, 16 ranks x 28 threads, auto-detected) | 534s |
| DeepVariant 1.9.0 (region-sharded, 16 bins x 32 shards) | 708s |
| **Total FASTQ -> VCF** | **~20.7 min** |

Running against NFS-backed storage instead of local disk adds significant I/O overhead to both the alignment (pragzip indexing) and DeepVariant stages; use local/NVMe scratch for input, output, and intermediate files where possible.

| Type | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | Recall | Precision | Frac_NA | F1_Score |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| INDEL | ALL | 467702 | 456097 | 11605 | 906287 | 1634 | 430822 | 1159 | 223 | 0.975187 | 0.996563 | 0.475370 | 0.985759 |
| INDEL | PASS | 467702 | 456097 | 11605 | 906287 | 1634 | 430822 | 1159 | 223 | 0.975187 | 0.996563 | 0.475370 | 0.985759 |
| SNP | ALL | 3254386 | 3165552 | 88834 | 3683913 | 5805 | 511273 | 2429 | 330 | 0.972703 | 0.998170 | 0.138785 | 0.985272 |
| SNP | PASS | 3254386 | 3165552 | 88834 | 3683913 | 5805 | 511273 | 2429 | 330 | 0.972703 | 0.998170 | 0.138785 | 0.985272 |

<p align="center">
<img src="https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/blob/main/images/deepvariant-fq2vcf.jpg"/a></br>
</p> 

# Using Dockerfile  (Single Node)  
### 1. Download the code :  

```bash
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/
```
### 2. Build the Docker Images
Part I: fq2bams
```bash
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy -t fq2bams -f Dockerfile_fq2bams .  
```
Part II: bams2vcf
```bash
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy -t bams2vcf -f Dockerfile_bams2vcf  .   
```

### 3. Run the Dockers  
Notes:  
<refdir> is expected to contain the bwa-mem2 index. You can index the reference during the run by enabling "--rindex" to fq2bams commandline.  

```bash
docker run  --volume <refdir>:/refdir <readsdir>:/readsdir <outdir_fq2bams>:/outdir fq2bams:latest python run_fq2bams.py --ref /refdir/<reference_file> --reads  /readsdir/<read1>  /readsdir/<read2>  --output /outdir/<outBAMfile>   

docker run  --volume <refdir>:/refdir <outdir_fq2bams>:/indir <output>:/outdir  bams2vcf:latest python run_bams2vcf.py --ref /refdir/<reference_file> --input /indir/  --output /outdir/<outVCFfile>   
```

# Results

Our latest results are published in this [blog](https://community.intel.com/t5/Blogs/Tech-Innovation/Artificial-Intelligence-AI/Intel-Xeon-is-all-you-need-for-AI-inference-Performance/post/1506083).


# Instructions to run the pipeline on an AWS ec2 instance (Single Node)
The following instructions run seamlessly on a standalone AWS ec2 instance. To run the following steps, create an ec2 instance with Ubuntu-22.04 having at least 60GB of memory and 500GB of disk. The input reference sequence and the paired-ended read datasets must be downloaded and stored on the disk.

### One-time setup
This step takes around ~15 mins to execute. During the installation process, whenever prompted for user input, it is recommended that the user select all default options.
```bash
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/
```

### Create the index files for the reference sequence
```bash
bash create_reference_index.sh
```

### Run the pipeline.
```bash
python run_fq2bams.py --ref refdir/<reference_file> --reads  readsdir/<read1>  readsdir/<read2>  --output outdir/<outBAMfile>     
python run_bams2vcf.py --ref refdir/<reference_file> --input outdir/  --output outdir/<outVCFfile>     

```


# Instructions to run the pipeline on an AWS ParallelCluster (Multi-node)   

The following instructions run seamlessly on AWS ParallelCluster. To run the following steps, first create an AWS parallelCluster as follows,
- Cluster setup: follow these steps to setup an [AWS ParallelCluster](https://docs.aws.amazon.com/parallelcluster/v2/ug/what-is-aws-parallelcluster.html).  Please see example [config file](scripts/aws/pcluster_example_config) to setup pcluster (_config_ file resides at ~/.parallelcluster/config/ on local machine). Please note: for best performance use shared file system with Amazon EBS _volume\_type = io2_ and _volume\_iops = 64000_ in the config file.
- Create pcluster: pcluster create <cluster_name>
- Login: login to the pcluster host/head node using the IP address of the cluster created in the previous step  
- Datasets: The input reference sequence and the paired-ended read datasets must be downloaded and stored in the _/sharedgp_ (pcluster shared directory defined in the config file) folder.


### One-time setup
This step takes around ~15 mins to execute. During the installation process, whenever prompted for user input, it is recommended that the user select all default options.
```bash
cd /sharedgp
wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz
tar -xzf Source_code_with_submodules.tar.gz
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/scripts/aws
bash deepvariant_setup.sh
```
### Modify _config_ file
We need a reference sequence and paired-ended read datasets. Open the "_config_" file and set the input and output directories as shown in config file.
The sample config contains the following lines to be updated.
```bash
export LD_PRELOAD=<absolute_path>/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/libmimalloc.so.2.0:$LD_PRELOAD
export INPUT_DIR=/path-to-read-datasets/
export OUTPUT_DIR=/path-to-output-directory/
export REF_DIR=/path-to-ref-directory/
REF=ref.fasta
R1=R1.fastq.gz
R2=R2.fastq.gz
```

### Create the index files for the reference sequence
```bash
bash pcluster_reference_index.sh
```

### Allocate compute nodes and install the prerequisites into the compute nodes.
```bash
bash pcluster_compute_node_setup.sh <num_nodes> <allocation_time>
# num_nodes: The number of compute nodes to be used for distributed multi-node execution.
# allocation_time: The maximum allocation time for the compute nodes in "hh:mm:ss" format. The default value is 2 hours, i.e., 02:00:00.
# Example command for allocating 4 nodes for 3 hours -
bash pcluster_compute_node_setup.sh 4 "03:00:00"
```

### Run the pipeline.
Note that the script uses default setting for creating multiple MPI ranks based on the system configuration.
```bash
bash run_pipeline_pcluster.sh
```

### Delete Cluster
```bash
pcluster delete <cluster_name>
```
