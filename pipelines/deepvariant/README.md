# Deepvariant Pipeline
### This repo presents OpenOmics DeepVariant Pipeline: A highly optimized scalable Deep Learning based short-read vaariant calling on x86 CPU clusters. The pipeline comprises of highly optimized: 1. bwa-mem2 for sequence mapping, 2. samtools for sorting output of bwa-mem2, and 3. DeepVariant for variant calling using sorted BAM records out of sorting.

### 0. Pipeline tools' location:   
bwa-mem2,samtools C/C++ based tools residing in:
```Open-Omics-Acceleration-Framework/applications/ ```.
DeepVariant is used as a Docker image. It is currently not available on Dockerhub. To use it, the user must build an image, convert it to a tar file, and distribute it across the cluster nodes. 
   * Prerequisite : 
	* Docker / Podman
	* gcc >= gcc 8.5.0 
	* MPI
	* make >=4.3
	* autoconf >= 2.69
	* zlib1g-dev
	* libncurses5-dev
	* libbz2-dev
	* liblzma-dev
   * All the script by default supports docker. If you are using podman use:  **alias docker=podman**

### 1. Download Code:
```bash
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
cd Open-Omics-Acceleration-Framework
git submodule update --init --recursive
```

### 2. Setting Envionment and Deepvariant Image
```bash
#Tested with Ubuntu 22.04.2 LTS
source setup_env.sh  my_env # Setting environment with name my_env. 
```
### 3. Cluster setup:  
#### 3.1.  Cluster using slurm job scheduler.
```bash
salloc --ntasks=1 --partition=<> --constraint=<machine type> --cpus-per-task=<cpus> --time=<node allocation time>
srun hostname > hostfile  
```  

#### 3.2 Standalone machine
```bash
hostname > hostfile
```

### 4. Compilation of tools , creation and distribution of docker/podman image in the allocated nodes.
```bash
# if you are using single node comment "bash load_deepvariant.sh" in the below script
source setup.sh 
```
Note: It takes ~20 mins to create docker image. 

### 5. Run the following script after the image is loaded on all the compute nodes listed in hostfile.  
Usage: sh run_pipline.sh <#ranks> <#ppn>  
* ranks: Number of mpi process that we want the pipeline to run on  
* ppn: mpi process per compute node. If we are running 2 ranks and provide ppn as 2, then both the ranks will run on 1 compute code (assuming dual socket machine, it will run 1 rank per socket)  

* Note: 
	* Before running the code, copy read and reference files to INPUT_DIR
	* for best performance, we advice to run 4 ranks per socket on spr nodes. So, assuming dual-socket compute node you can run 8 ranks on 1 compute node.  
	* It takes ~20 mins to create docker image. It might break while installation because of your proxy settings. [Follow](https://docs.docker.com/network/proxy/) instructions to build image using proxy.
 
```bash 
export INPUT_DIR=./    # This directory contains Reference and Read files.
export OUTPUT_DIR=./   # This directory contains intermediate and log files.
ranks=8 
ppn=8
sh run_pipeline.sh $ranks $ppn #change your reference read file names inside the script.
```

# Results

For detailed information please refer to [blog](). 
