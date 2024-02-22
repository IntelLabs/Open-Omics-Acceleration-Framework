# Steps to use to reproduce benchmarking results of Open Omics Acceleration Framework v2.1 on AWS EC2 instances published in the AWS Intel blog

## Preparing an AWS EC2 instance for benchmarking

1.	Log in to your AWS account
2.	Launch a virtual machine with EC2. The configuration details of instances used for each pipeline are mentioned in the respective sections below.
  * Choose an Amazon Machine Image (AMI): From “Quick Start”, select "Ubuntu Server 22.04 LTS (HVM), SSD Volume Type".
  * Choose an Instance Type.
  * Configure Instance.
  * Add Storage: You can add storage based on the workload requirements
  * Configure security group
  * Review and launch the instance (ensure you have/create a key to ssh login in next step)
3.	Use SSH to login to the machine after the instance is up and running
  * $ ssh -i <key.pem> username@Public-DNS
4.	The logged in AWS instance machine is now ready to use – you can download Open Omics workloads and related datasets to be executed on this instance.


## AlphaFold2-based Protein Folding Pipeline

### Configuration Details

**BASELINE on m7i.24xlarge**: Test by Intel as of <11/30/23>. 1 instance, 1-socket, 1x Intel® Xeon® Platinum 8488C, 48 cores, HT On, Turbo On, Total Memory 384GB (1 slot/ 384GB/ DDR5 4800 MT/s), bios: Amazon EC2 v 1.0, ucode version: 0x2b0004b1, OS Version: Ubuntu 22.04.3 LTS, kernel version: 6.2.0-1017-aws, compiler version: g++ 11.4.0, PyTorch - v1.12.1, JAX - v0.4.14, OpenFold - v 1.0.1, Hmmer - v3.3.2, hh-suite - v3.3.0, Kalign2 – v2.04, model name & version: AlphaFold2

**Open Omics on m7i.24xlarge**: Test by Intel as of <11/30/23>. 1 instance, 1-socket, 1x Intel® Xeon® Platinum 8488C, 48 cores, HT On, Turbo On, Total Memory 384GB (1 slot/ 384GB/ DDR5 4800 MT/s), bios: Amazon EC2 v 1.0, ucode version: 0x2b0004b1, OS Version: Ubuntu 22.04.3 LTS, kernel version: 6.2.0-1017-aws, compiler version: g++ 11.4.0, workload version: Intel-python - 2022.1.0, JAX - v0.4.21, Open Omics Acceleration Framework v2.1, Open Omics AlphaFold2, - v1.0, IntelLabs Hmmer v1.0, IntelLabs hh-suite v1.0), Kalign2 – v2.04, framework version: PyTorch – v2.0.1, model name & version: AlphaFold2

**Open Omics on m7i.48xlarge**: Test by Intel as of <11/30/23>. 1 instance, 2-sockets, 1x Intel® Xeon® Platinum 8488C, 48 cores, HT On, Turbo On, Total Memory 384GB (1 slot/ 384GB/ DDR5 4800 MT/s), bios: Amazon EC2 v 1.0, ucode version: 0x2b0004b1, OS Version: Ubuntu 22.04.3 LTS, kernel version: 6.2.0-1017-aws, compiler version: g++ 11.4.0, workload version: Intel-python - 2022.1.0, JAX - v0.4.21, Open Omics Acceleration Framework v2.1, Open Omics AlphaFold2, - v1.0, IntelLabs Hmmer v1.0, IntelLabs hh-suite v1.0), Kalign2 – v2.04, framework version: PyTorch – v2.0.1, model name & version: AlphaFold2


### Step by step instructions to benchmark baseline and Open Omics Acceleration Framework

```sh
 cd ~
 git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
```

- Test dataset can be downloaded from https://www.uniprot.org/proteomes/UP000001940. Click on 'Download' and select options **Download only reviewed (Swiss-Prot:) canonical proteins (4,463)**, Format: **Fasta** and Compressed: **No**.

- Save the file as 'uniprotkb_proteome.fasta' inside  folder ~/Open-Omics-Acceleration-Framework/benchmarking/AWS-Intel-blog-v2.1-2024/


#### Baseline ([OpenFold](https://github.com/aqlaboratory/openfold))

EC2Instance: m7i.24xlarge

```sh
mkdir -p ~/data/models
bash ~/Open-Omics-Acceleration-Framework/applications/alphafold/alphafold/scripts/download_alphafold_params.sh ~/data/models/
cd ~/Open-Omics-Acceleration-Framework/
#Below script generate dataset folder: ~/celegans_samples contains 77 files and ~/celegans_samples_long contain 18 files
cd ~/Open-Omics-Acceleration-Framework/benchmarking/AWS-Intel-blog-v2.1-2024/
python3 proteome.py
cd ~
git clone --recursive https://github.com/aqlaboratory/openfold.git
cd openfold
#use mamba for faster depedency solve
mamba env create -n openfold_env -f environment.yml
conda activate openfold_env
bash scripts/download_alphafold_dbs.sh ~/data/
fastadir=./celegans_samples
#Run script
python3 run_pretrained_openfold.py \
    $fasta_dir \
    ~/data/pdb_mmcif/mmcif_files/ \
    --uniref90_database_path ~/data/uniref90/uniref90.fasta \
    --mgnify_database_path ~/data/mgnify/mgy_clusters_2018_12.fa \
    --pdb70_database_path ~/data/pdb70/pdb70 \
    --uniclust30_database_path ~/data/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    --bfd_database_path ~/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --jackhmmer_binary_path ~/miniconda3/envs/openfold_venv/bin/jackhmmer \
    --hhblits_binary_path ~/miniconda3/envs/openfold_venv/bin/hhblits \
    --hhsearch_binary_path ~/miniconda3/envs/openfold_venv/bin/hhsearch \
    --kalign_binary_path ~/miniconda3/envs/openfold_venv/bin/kalign \
    --config_preset "model_1" \
    --model_device "cpu" \
    --output_dir ./ \
    --jax_param_path ~/models/params/params_model_1.npz \
    --skip_relaxation  \
    --cpus 96 \
    --long_sequence_inference
```
Note 1: If you are using a different instance, modify --cpus flag to available vcpus. 
Note 2: Use --long_seququence_inference option if you are running long sequences.


#### Open Omics Acceleration Framework alphafold2-based-protein-folding pipeline

EC2Instance: m7i.24xlarge, m7i.48xlarge

Follow [steps](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/alphafold2-based-protein-folding) for more details:

```sh
cd ~/Open-Omics-Acceleration-Framework/pipelines/alphafold2-based-protein-folding
docker build -t alphafold .           # Build a docker image named alphafold
export DATA_DIR=~/data
export SAMPLES_DIR=<fasta_dir>
export OUTPUT_DIR=<path-to-output-directory>
export LOG_DIR=<path-to-log-directory>
docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:latest
```

## DeepVariant-based germline variant calling pipeline (fq2vcf)

### Configuration Details

**BASELINE on c7i.24xlarge**: Test by Intel as of <11/30/23>. 1 instance, 1-socket, 1x Intel® Xeon® Platinum 8488C, 48 cores, HT On, Turbo On, Total Memory 192GB (1 slot/ 192GB/ DDR5 4800 MT/s), bios: Amazon EC2 v 1.0, ucode version: 0x2b0004b1, OS Version: Ubuntu 22.04.3 LTS, kernel version: 6.2.0-1017-aws, compiler version: g++ 11.4.0, workload version: bwa-mem v0.7.17, Samtools v. 1.16.1, DeepVariant v1.5, framework version: Intel-tensorflow 2.11.0, model name & version: Inception V3

**Open Omics on c7i.24xlarge**: Test by Intel as of <11/30/23>. 1 instance, 1-socket, 1x Intel® Xeon® Platinum 8488C, 48 cores, HT On, Turbo On, Total Memory 192GB (1 slot/ 192GB/ DDR5 4800 MT/s), bios: Amazon EC2 v 1.0, ucode version: 0x2b0004b1, OS Version: Ubuntu 22.04.3 LTS, kernel version: 6.2.0-1017-aws, compiler version: g++ 11.4.0, workload version: bwa-mem2 v2.2.1, Samtools v. 1.16.1, Open Omics DeepVariant v1.0, Open Omics Acceleration Framework v.2.1, framework version: Intel-tensorflow 2.11.0, model name & version: Inception V3

**Open Omics on c7i.48xlarge**: Test by Intel as of <11/30/23>. Up to 8 instances, 16-sockets; Each socket has 1x Intel® Xeon® Platinum 8488C, 48 cores, HT On, Turbo On, Total Memory 192GB (1 slot/ 192GB/ DDR5 4800 MT/s), bios: Amazon EC2 v 1.0, ucode version: 0x2b0004b1, OS Version: Ubuntu 22.04.3 LTS, kernel version: 6.2.0-1017-aws, compiler version: g++ 11.4.0, workload version: bwa-mem2 v2.2.1, Samtools v. 1.16.1, Open Omics DeepVariant v1.0, Open Omics Acceleration Framework v.2.1, framework version: Intel-tensorflow 2.11.0, model name & version: Inception V3


### Step by step instructions to benchmark baseline and Open Omics Acceleration Framework

#### Dataset 
```sh
mkdir -p ~/HG001
wget https://genomics-benchmark-datasets.s3.amazonaws.com/google-brain/fastq/novaseq/wgs_pcr_free/30x/HG001.novaseq.pcr-free.30x.R1.fastq.gz -P ~/HG001/
wget https://genomics-benchmark-datasets.s3.amazonaws.com/google-brain/fastq/novaseq/wgs_pcr_free/30x/HG001.novaseq.pcr-free.30x.R2.fastq.gz -P ~/HG001/
wget https://broad-references.s3.amazonaws.com/hg38/v0/Homo_sapiens_assembly38.fasta -P ~/HG001/

```

#### Baseline

EC2Instance: c7i.24xlarge
Prerequisite : docker/podman 
```sh
cd ~
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
#pull deepvariant docker image
docker pull google/deepvariant:1.5.0

cd ~/Open-Omics-Acceleration-Framework/benchmarking/AWS-Intel-blog-v2.1-2024/
#copy baseline code
cp test_pipe_bwa.py ../../pipelines/deepvariant-based-germline-variant-calling-fq2vcf/
cp run_pipe_bwa.sh ../../pipelines/deepvariant-based-germline-variant-calling-fq2vcf/

#clone bwa repo
cd ../../applications/
wget https://github.com/lh3/bwa/archive/refs/tags/v0.7.17.tar.gz
tar -xvzf v0.7.17.tar.gz
cd bwa-0.7.17
make

# compile htslib
cd ../htslib
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
#make install   #uncomment this for installation

# compile samtools
cd ../samtools
autoheader
autoconf -Wno-syntax
chmod 775 configure
./configure           # Needed for choosing optional functionality
make
cd ..

cd ../../pipelines/deepvariant-based-germline-variant-calling-fq2vcf/

#create index for bwa
../../applications/bwa-0.7.17/bwa index ~/HG001/Homo_sapiens_assembly38.fasta

#run pipeline
bash run_pipe_bwa.sh
```
#### Open Omics Acceleration Framework deepvariant-based-germline-variant-calling-fq2vcf pipeline

EC2Instance: c7i.24xlarge, c7i.48xlarge

pcluster: 2 x c7i.48xlarge, 4 x c7i.48xlarge, 8 x c7i.48xlarge,

To run on c7i.24xlarge, c7i.48xlarge follow [link](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/deepvariant-based-germline-variant-calling-fq2vcf#instructions-to-run-the-pipeline-on-an-aws-ec2-instance).

To run on 2 x c7i.48xlarge, 4 x c7i.48xlarge, 8 x c7i.48xlarge follow [link](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/deepvariant-based-germline-variant-calling-fq2vcf#instructions-to-run-the-pipeline-on-an-aws-parallelcluster).



## Single-cell RNA-seq analysis pipeline

### Configuration Details

**BASELINE on r7i.24xlarge**: Test by Intel as of <11/30/23>. 1 instance, 1-socket, 1x Intel® Xeon® Platinum 8488C, 48 cores, HT On, Turbo On, Total Memory 768GB (1 slot/ 768GB/ DDR5 4800 MT/s), bios: Amazon EC2 v 1.0, ucode version: 0x2b0004b1, OS Version: Ubuntu 22.04.3 LTS, kernel version: 6.2.0-1017-aws, compiler version: g++ 11.4.0, workload version: scanpy v 1.9.1

**Open Omics on r7i.24xlarge**: Test by Intel as of <11/30/23>. 1 instance, 1-socket, 1x Intel® Xeon® Platinum 8488C, 48 cores, HT On, Turbo On, Total Memory 768GB (1 slot/ 768GB/ DDR5 4800 MT/s), bios: Amazon EC2 v 1.0, ucode version: 0x2b0004b1, OS Version: Ubuntu 22.04.3 LTS, kernel version: 6.2.0-1017-aws, compiler version: g++ 11.4.0, workload version: Open Omics Acceleration Framework v.2.1

**Open Omics on c7i.24xlarge**: Test by Intel as of <11/30/23>. 1 instance, 1-socket, 1x Intel® Xeon® Platinum 8488C, 48 cores, HT On, Turbo On, Total Memory 192GB (1 slot/ 192GB/ DDR5 4800 MT/s), bios: Amazon EC2 v 1.0, ucode version: 0x2b0004b1, OS Version: Ubuntu 22.04.3 LTS, kernel version: 6.2.0-1017-aws, compiler version: g++ 11.4.0, workload version: Open Omics Acceleration Framework v.2.1 

### Step by step instructions to benchmark baseline and Open Omics Acceleration Framework

#### Baseline ([rapids-single-cell-examples](https://github.com/NVIDIA-Genomics-Research/rapids-single-cell-example))

EC2Instance: r7i.24xlarge

```sh
git clone https://github.com/NVIDIA-Genomics-Research/rapids-single-cell-examples.git
cd rapids-single-cell-examples
conda env create --name rapidgenomics -f conda/cpu_notebook_env.yml
conda activate rapidgenomics
mkdir data
wget -P ./data https://rapids-single-cell-examples.s3.us-east-2.amazonaws.com/1M_brain_cells_10X.sparse.h5ad
cd notebooks
python -m ipykernel install --user --display-name "Python (rapidgenomics)"
```
Note: Open Jupyter notebook and Select **1M_brain_cpu_analysis.ipynb** file and run all cells.

#### Open Omics Acceleration Framework single-cell-RNA-seq-analysis pipeline

EC2Instance: r7i.24xlarge and c7i.24xlarge

Follow [link](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/single-cell-RNA-seq-analysis#run-with-jupyter-notebook-interactive) to run Open Omics Accelaration Framework single-cell-RNA-seq-analysis pipeline using interactive mode.

# Cleanup

Terminate all EC2 instances used to run benchmarks to avoid incurring charges.                                                                                      
