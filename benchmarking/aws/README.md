# Benchmarking of Open Omics Acceleration Framework on AWS

Step by step commands to benchmark OpenOmics framework on AWS

1.	Log in to your AWS account
2.	Launch a virtual machine with EC2
  * Choose an Amazon Machine Image (AMI): Select any 64-bit (x86) AMI  (say, Ubuntu Server 20.04 LTS) from “Quick Start”.
  * Choose an Instance Type.
  * Configure Instance.
  * Add Storage: You can add storage based on the workload requirements
  * Configure security group 
  * Review and launch the instance (ensure you have/create a key to ssh login in next step)
3.	Use SSH to login to the machine after the instance is up and running
  * $ ssh -i <key.pem> username@Public-DNS
4.	The logged in AWS instance machine is now ready to use – you can download OpenOmics workloads and related datasets to be executed on this instance.

## Machine configurations used for benchmarking
AWS c5.12xlarge: 1-instance AWS c5.12xlarge: 48 vCPUs (Cascade Lake), 96 GB total memory, ucode: 0x500320a, Ubuntu 22.04, 5.15.0-1004-aws

AWS m5.12xlarge: 1-instance AWS c5.12xlarge: 48 vCPUs (Cascade Lake), 192 GB total memory, ucode: 0x500320a, Ubuntu 22.04, 5.15.0-1004-aws

AWS c6i.16xlarge: 1-instance AWS c6i.16xlarge: 64 vCPUs (Ice Lake), 128 GB total memory, ucode: 0xd000331, Ubuntu 22.04, 5.15.0-1004-aws

AWS m6i.16xlarge: 1-instance AWS m6i.16xlarge: 64 vCPUs (Ice Lake), 256 GB total memory, ucode: 0xd000331, Ubuntu 22.04, 5.15.0-1004-aws


# Step by step instructions to benchmark baseline (bwa-mem) and OpenOmics BWA-MEM (bwa-mem2) on m5.12xlarge and m6i.16xlarge instances of AWS
## Step 1: Download datasets
Download reference genome: Homo_sapiens_assembly38.fasta
```sh
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
```

Download dataset #1: ERR194147 (Paired-End) from 
s3://sra-pub-run-odp/sra/ERR194147/ERR194147

Download dataset #2: ERR1955529 (Paired-End) from 
s3://sra-pub-run-odp/sra/ERR1955529/ERR1955529

Download dataset #3: ERR3239276 (Paired-End) from 
s3://sra-pub-run-odp/sra/ERR3239276/ERR3239276

## Step 2: Download and compile baseline (BWA v0.7.17) and create index
```sh
curl -L https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar jxf -
cd bwa-0.7.17 && make
./bwa index hs_asm38/Homo_sapiens_assembly38.fasta
cd ..
```

## Step 3: Download OpenOmics BWA-MEM (BWA-MEM2 v2.2.1) and create index
```sh
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
cd bwa-mem2-2.2.1_x86-linux && ./bwa-mem2 index hs_asm38/Homo_sapiens_assembly38.fasta
cd ..
```

## Step 4: Running baseline and OpenOmics BWA-MEM
The m5.12xlarge instance has 24 cores, while the m6i.16xlarge instance has 32 cores. The memory available on the instances allows for using both the threads on each core.

### Run baseline BWA-MEM
```sh
cd bwa-0.7.17
```
Sample commands shown below:\
For m5.12xlarge:
```sh
./bwa mem -t 48 hs_asm38/Homo_sapiens_assembly38.fasta ERR194147_1.fastq ERR194147_2.fastq > ERR194147.out.sam
```
For m6i.16xlarge:
```sh
./bwa mem -t 64 hs_asm38/Homo_sapiens_assembly38.fasta ERR194147_1.fastq ERR194147_2.fastq > ERR194147.out.sam
```

### Run OpenOmics BWA-MEM
```sh
cd ../bwa-mem2-2.2.1_x86-linux
```
Sample commands shown below:\
For m5.12xlarge:
```sh
./bwa-mem2 mem -t 48 hs_asm38/Homo_sapiens_assembly38.fasta ERR194147_1.fastq ERR194147_2.fastq > ERR194147.out.sam
```
For m6i.16xlarge:
```sh
./bwa-mem2 mem -t 64 hs_asm38/Homo_sapiens_assembly38.fasta ERR194147_1.fastq ERR194147_2.fastq > ERR194147.out.sam
```


# Step by step instructions to benchmark baseline (minimap2) and OpenOmics minimap2 (mm2-fast) on c5.12xlarge and c6i.16xlarge instances of AWS

## Step 1: Download datasets
Download reference genome
```sh
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
```
Link to download HG002 ONT Guppy 3.6.0 dataset:
https://precision.fda.gov/challenges/10/view \
File name: HG002_GM24385_1_2_3_Guppy_3.6.0_prom.fastq.gz

Link to download HG002 HiFi 14kb-15kb dataset:
https://precision.fda.gov/challenges/10/view \
File name: HG002_35x_PacBio_14kb-15kb.fastq.gz

Download HG002 CLR dataset from 
s3://giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_PacBio_MtSinai_NIST_subreads_fasta_10082018

Download hap2 assembly dataset- 
```sh
wget https://zenodo.org/record/4393631/files/NA24385.HiFi.hifiasm-0.12.hap2.fa.gz
```

## Step 2: Download and compile baseline (minimap2 v0.2.22)
```sh
git clone https://github.com/lh3/minimap2.git -b v2.22
cd minimap2 && make
```

## Step 3: Run baseline minimap2
```sh
./minimap2 -ax [preset] [ref-seq] [read-seq] -t 48 > minimap2output
```
Example command for ONT HG002 dataset:
```sh
./minimap2 -ax map-ont  GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz HG002_ONT.fastq -t 48 > minimap2output
```

## Step 3: Download and compile OpenOmics minimap2 (mm2-fast)
```sh
git clone --recursive https://github.com/bwa-mem2/mm2-fast.git -b mm2-fast-v2.22 mm2-fast-contrib
cd mm2-fast-contrib && make multi
```

## Step 4: Create index for OpenOmics minimap2
```sh
./build rmi.sh path-to-ref-seq <preset flags>
<preset flags> are as follows:
ONT: map-ont
HiFi: map-hifi
CLR: map-pb
Assembly: asm5
```
Example: Create OpenOmics minimap2 index for ONT datasets for GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
```sh
./build rmi.sh GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz map-ont
```

## Step 5: Run OpenOmics minimap2
```sh
./mm2-fast -ax [preset] [ref-seq] [read-seq] -t [num_threads] > mm2-fastoutput
```
Example command to run HG002 ONT dataset on c5.12xlarge
```sh
./mm2-fast -ax map-ont GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz HG002_ONT.fastq -t 48 > mm2-fastoutput
```
Example command to run HG002 ONT dataset on c6i.16xlarge
```sh
./mm2-fast -ax map-ont GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz HG002_ONT.fastq -t 64 > mm2-fastoutput
```

# Step by step instructions to benchmark OpenOmics for ATAC-Seq data analysis on multiple c5.24xlarge and c6i.32xlarge instances of AWS

## Update
```bash
apt-get update
```

## Install software
```bash
apt-get install -y git
apt-get install -y libcurl4-openssl-dev
apt-get install -y hdf5-tools
apt-get install -y rsync
apt-get install -y make
apt-get install -y gcc
apt-get install -y libblas-dev
apt-get install -y python3.7 python3-pip
ln -nsf /usr/bin/python3.7 /usr/bin/python
```
## Anaconda Environment
```bash
conda create --name Atac python=3.7
conda activate Atac
```

## Clone the libxsmm repository and set library path
```bash
cd /home/
git clone https://github.com/libxsmm/libxsmm.git
cd /home/libxsmm
git checkout b3da2b1bed9d27f9d6bae91a683f8cf76fe299b5
make -j                   # Use AVX=2 for AVX2 and AVX=3 for AVX512
cd /home/               
export LD_LIBRARY_PATH=/home/libxsmm/lib/
```

## Clone atacworks repo
```bash
git clone --branch v0.2.0 https://github.com/clara-parabricks/AtacWorks.git
```

## Clone the OpenOmics version
```bash
git clone https://github.com/IntelLabs/Trans-Omics-Acceleration-Library.git
```

## Apply patch
```bash
cd  /home/AtacWorks/
git apply /home/Trans-Omics-Acceleration-Library/applications/ATAC-Seq/AtacWorks_cpu_optimization_patch.patch
```

## Install python packages
```bash
python3.7 -m pip install -r requirements-base.txt
python3.7 -m pip install torch torchvision torchaudio
python3.7 -m pip install -r requirements-macs2.txt
```

## (Optional) Install torch-ccl
```bash
# Install torch-ccl
# git clone --branch v1.1.0 https://github.com/intel/torch-ccl.git && cd torch-ccl
# git submodule sync
# git submodule update --init --recursive
# python3.7 setup.py install
```

## Install 1D convolution module
```bash
cd /home/libxsmm/samples/deeplearning/conv1dopti_layer/Conv1dOpti-extension/ 
python setup.py install
```

## Install AtacWorks folder ans set path
```bash
cd /home/AtacWorks/
python3.7 -m pip install .
atacworks=/home/AtacWorks/
```

## Download data to train
```bash
wget https://atacworks-paper.s3.us-east-2.amazonaws.com/dsc_atac_blood_cell_denoising_experiments/50_cells/train_data/noisy_data/dsc.1.Mono.50.cutsites.smoothed.200.bw
wget https://atacworks-paper.s3.us-east-2.amazonaws.com/dsc_atac_blood_cell_denoising_experiments/50_cells/train_data/clean_data/dsc.Mono.2400.cutsites.smoothed.200.bw
wget https://atacworks-paper.s3.us-east-2.amazonaws.com/dsc_atac_blood_cell_denoising_experiments/50_cells/train_data/clean_data/dsc.Mono.2400.cutsites.smoothed.200.3.narrowPeak
```

## Download file conversion binaries and set path
```bash
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedGraphToBigWig /home/
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bigWigToBedGraph /home/
export PATH="$PATH:/home/" >> /home/.bashrc         # set the path for bedGraphToBigWig binaries 
```

## Data preprocessing

```python
python $atacworks/scripts/peak2bw.py \
    --input dsc.Mono.2400.cutsites.smoothed.200.3.narrowPeak \
    --sizes $atacworks/data/reference/hg19.chrom.sizes \
    --out_dir ./ \
    --skip 1

python $atacworks/scripts/get_intervals.py \
    --sizes $atacworks/data/reference/hg19.auto.sizes \
    --intervalsize 50000 \
    --out_dir ./ \
    --val chr20 \
    --holdout chr10

python $atacworks/scripts/bw2h5.py \
        --noisybw dsc.1.Mono.50.cutsites.smoothed.200.bw \
        --cleanbw dsc.Mono.2400.cutsites.smoothed.200.bw \
        --cleanpeakbw dsc.Mono.2400.cutsites.smoothed.200.3.narrowPeak.bw \
        --intervals training_intervals.bed \
        --out_dir ./ \
        --prefix Mono.50.2400.train \
        --pad 5000 \
        --nonzero

python $atacworks/scripts/bw2h5.py \
        --noisybw dsc.1.Mono.50.cutsites.smoothed.200.bw \
        --cleanbw dsc.Mono.2400.cutsites.smoothed.200.bw \
        --cleanpeakbw dsc.Mono.2400.cutsites.smoothed.200.3.narrowPeak.bw \
        --intervals val_intervals.bed \
        --out_dir ./ \
        --prefix Mono.50.2400.val \
        --pad 5000
```

## Set affinity and threads
```bash
export KMP_AFFINITY=compact,1,0,granularity=fine
export LD_PRELOAD=/home/libtcmalloc.so              # Copy these files in the /home folder first
export LD_PRELOAD=/home/libjemalloc.so
export OMP_NUM_THREADS=31                           # (Available cores (N) - 1)
```

## Training run (Single Socket)
```python
# In numactl command, "-C 1-31" is for running on cores 1 to 31. 
# General case for an N core machine is "-C 1-(N-1)".
# Keep batch size in config/train_config.yaml to a multiple of (N-1) for optimum performance

numactl --membind 0 -C 1-31 python $atacworks/scripts/main.py train \
        --config configs/train_config.yaml \
        --config_mparams configs/model_structure.yaml \
        --files_train $atacworks/Mono.50.2400.train.h5 \
        --val_files $atacworks/Mono.50.2400.val.h5                           
```

Option - Another option to use on machines without NUMA --- "taskset -c 1-31 python ..."

## Training run (Multiple Sockets/Nodes)

``` bash
export OMP_NUM_THREADS=30                           # (Available cores (N) - 2)
# 1. change line 23 in configs/train_config.yaml with the following
#        dist-backend: 'gloo' 
# 2. change line 22 in configs/train_config.yaml with the following
#        dist-backend: 'gloo'
# 3. Keep batch size (bs) in config/train_config.yaml to a multiple of (N-2) for optimum performance. 
#    Batch size gets multiplied by number of socket. Hence, if bs=30, no. of sockets = 16 than batch size = 30*16 = 480
# 4. Comment line the following line (79,80) in AtacWorks/claragenomics/dl4atac/utils.py and reinstall AtacWorks using "pip install ." command.
#       if (os.path.islink(latest_symlink)):
#               os.remove(latest_symlink)
# 5. Run the following Slurm batch script that uses MPI commands.

sbatch Batchfile_CPU.slurm

```

# Cleanup 

Terminate all EC2 instances used to run benchmarks to avoid incurring charges.
