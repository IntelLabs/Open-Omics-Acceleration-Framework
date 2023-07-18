# OpenOmics Deepvariant Pipeline
### OpenOmics Deepvariant Pipeline is a highly optimized, scalable, deep-learning-based short-read variant calling pipeline on x86 CPU clusters. The pipeline comprises of 1. bwa-mem2 (a highly optimized version of bwa-mem) for sequence mapping, 2. distributed SAM sorting using samtools, and 3. an optimized version of DeepVariant for variant calling.

### 0. Pipeline tools' location:   
The source code of bwa-mem2, samtools, and DeepVariant is residing in:
```Open-Omics-Acceleration-Framework/applications/ ```.
We build bwa-mem2 and samtools from source. For DeepVariant, we build a docker image and then use the built image while running the pipeline. Note that, the pre-built image is not available on dockerhub and the image needs to be built from source.
The following are pre-requisites for the pipeline.
   * Prerequisite :
        * Docker / Podman
        * gcc >= 8.5.0
        * MPI
        * make >= 4.3
        * autoconf >= 2.69
        * zlib1g-dev
        * libncurses5-dev
        * libbz2-dev
        * liblzma-dev
   * We are providing setup and installation scripts located at _Open-Omics-Acceleration-Framework/pipelines/deepvariant_. Note that, all the script by default are written for docker.

### 1. Clone the repo:
```bash
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
```

### 2. Setting Envionment
```bash
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant/
#Tested with Ubuntu 22.04.2 LTS
source setup_env.sh  my_env # Setting environment with name _my_env_. 
```
### 3. Cluster setup:  
#### 3.1.  Cluster using slurm job scheduler.
```bash
salloc --ntasks=1 --partition=<> --constraint=<machine type> --cpus-per-task=<cpus> --time=<node allocation time>
srun hostname > hostfile
srun --pty bash    ## login to a compute node from hostfile
```  

#### 3.2 Standalone machine
The pipeline can also be tested on a single standalone machine.
```bash
hostname > hostfile
```

#### Activate conda environment into cluster/standalone machine

```bash
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant/
conda activate my_env
```

### 4. Compilation of tools, creation and distribution of docker/podman image on the allocated nodes.
To run docker without sudo acess follow [link](https://docs.docker.com/engine/install/linux-postinstall/), or use "_sudo docker_" as an argument in the below script instead of docker.
```bash
# If you are using single node, comment out line No. 56, i.e., "bash load_deepvariant.sh" of setup.sh.
source setup.sh docker/podman/"sudo docker"  
* docker/podman/"sudo docker" : optional argument. It takes docker by default.

```
Note: It takes ~30 mins to create the docker image. Docker build might break if you are behind a proxy. Example to provide proxy for building a docker image is shown in the setup.sh file. [Follow](https://docs.docker.com/network/proxy/) instructions for more details.

### 5. Run the following script after the image is loaded on all the compute nodes listed in the hostfile.  
Usage: sh run_pipline.sh <#ranks> <#ppn> <reference_seq.fasta> <read_r1.gz> <read_r2.gz> [docker/podman/"sudo docker"]

* ranks: Number of mpi processes that we want the pipeline to run with  
* ppn: The number of mpi processes per compute node. If we are running 2 ranks and provide ppn as 2, then both the ranks will run on 1 compute node (assuming dual socket machine, it will run 1 rank per socket)
* reference_seq.fasta: The name of reference genome sequence file.
* read_r1.gz: The name of the R1 file of paired-end reads, must be in .gz format.
* read_r2.gz: The name of the R2 file of paired-end reads, must be in .gz format.   
* docker/podman/"sudo docker" : optional argument. It takes docker by default.
* Note: 
	* Before running the code, copy the read and the reference files to INPUT_DIR.
	* For the best performance, we advice to run 4 ranks per socket on spr nodes. So, assuming dual-socket compute node you can run 8 ranks on 1 compute node.
 	* Pipeline will not run if there are fewer than 3 CPUs per rank. 	  

```bash 
export INPUT_DIR=./path-to-input    # This directory contains Reference and Read files.
export OUTPUT_DIR=./path-to-log-dir   # This directory contains intermediate and log files.
ranks=8 
ppn=8
sh run_pipeline.sh $ranks $ppn reference_seq.fasta R1.gz R2.gz podman # Change the arguments according to the user specific ones.
```

# Results

For detailed information, please refer to the [blog](). 
