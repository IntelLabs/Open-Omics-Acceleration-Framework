# fq2vcf: OpenOmics Deepvariant based Variant Calling Pipeline  
### Overview:  
OpenOmics's fq2vcf is a highly optimized, distributed, deep learning-based short-read germline variant calling pipeline for x86 CPUs. 
The pipeline comprises of:   
1. bwa-mem2 (a highly optimized version of bwa-mem) for sequence mapping  
2. SortSAM using samtools  
3. An optimized version of DeepVariant tools for variant calling.
The following figure illustrates the pipeline.

<p align="center">
<img src="https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/blob/main/images/deepvariant-fq2vcf.jpg"/a></br>
</p> 


# Using Dockerfile  (Single Node)  
### 1. Download the code :  

```bash
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/
```
### 2. Build Docker images
```bash
docker build -f Dockerfile_part1 -t deepvariant:part1 .      ## Part I: fq2bam
docker build -f Dockerfile_part2 -t deepvariant:part2 .      ## Part II: bam2vcf    
```

### 3. Setup Input Paramters through ```Config``` file<sup>1</sup>  
INPUT_DIR=/input       ## Location of input reads files (default location /reads in docker)    
OUTPUT_DIR=/out       ## Location of output files (default location /output in docker)  
REF_DIR=/ref          ## Location of reference sequence (default location /ref in docker)  
R1=HG001_1.fastq.gz         ## name of input reads file1  
R2=HG001_2.fastq.gz         ## name of input reads file2  
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna   ## name of reference sequence  

Update the above fields in ./config file and in extra_scripts/config file    
<sup>**1**</sup>  **All fields are mandatory**  

### 4. Run Docker images
```
docker run -v ./config:/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/scripts/aws/config -v <readdir>:/reads -v <refdir>:/ref -v <refdir>:/output -it deepvariant:part1 bash run_pipeline_ec2_part1.sh

docker run -v ./extra_scripts/config:/opt/deepvariant/config  -v <redir>:/ref -v <outdir>:/output -it deepvariant:part2 bash run_pipeline_ec2_part2.sh
```

# General Notes:
* The source code of bwa-mem2, samtools, and DeepVariant are residing in:
```Open-Omics-Acceleration-Framework/applications/ ```.
* We build bwa-mem2 and samtools from source; while for DeepVariant, we build a docker image and then use the built image while running the pipeline. Note that, the pre-built image is not available on dockerhub and the image needs to be built from source.
* We provide scripts for setting-up miniconda environment (setup_env.sh), compilation of the required pipeline tools, building & loading of DeepVariant docker image (setup.py). These scripts located at _Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf_. Note that, all the scripts are, by default, written for docker.  

* Prerequisite : The following are pre-requisites for the pipeline -  our scripts assumed these packages are already installed in the system.  
        * Docker / Podman  
        * gcc >= 8.5.0  
        * MPI  
        * make >= 4.3  
        * autoconf >= 2.69  
        * zlib1g-dev   
        * libncurses5-dev  
        * libbz2-dev  
        * liblzma-dev  

# Instructions to run the pipeline on on-prem (Single node & Multi-node)  
### 1. Download the latest release:  
```bash
wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/2.1/Source_code_with_submodules.tar.gz
tar -xzf Source_code_with_submodules.tar.gz
```

### 2. Setting Envionment
```bash
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/scripts/cluster/
#Tested with Ubuntu 22.04.2 LTS
source ../../setup_env.sh  dv_env # Setting environment with name dv_env.
```

### 3. Compute setup:  You can choose cluster (3.1) or standalone (3.2) mode for the run
#### 3.1.  Cluster using slurm job scheduler.
```bash
salloc --ntasks=1 --partition=<> --constraint=<machine type> --cpus-per-task=<cpus> --time=<node allocation time>
srun hostname > hostfile
```

#### 3.2 Standalone machine
The pipeline can also be tested on a single standalone machine.
```bash  
hostname > hostfile   
```

### 4. Compilation of tools, creation and distribution of docker/podman image on the allocated nodes.
To run docker without sudo access follow [link](https://docs.docker.com/engine/install/linux-postinstall/), or use "_sudo docker_" as an argument in the below script instead of docker.
```bash
# If you are using single node, comment out line No. 56, i.e., "bash load_deepvariant.sh" of setup.sh.
source setup.sh [docker | podman | "sudo docker"]
* docker/podman/"sudo docker" : optional argument. It takes docker by default.
```
Note: It takes ~30 mins to create the docker image. Docker build might break if you are behind a proxy. Example to provide proxy for building a docker image is shown in the setup.sh file. [Follow](https://docs.docker.com/network/proxy/) instructions for more details.

### 5. Create _config_ file
We need a reference sequence and paired-ended read datasets. Open the "config" file and set the input and output directories as shown in config file. The sample config contains the following lines to be updated.  

```bash
export LD_PRELOAD=<absolute_path>/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/libmimalloc.so.2.0:$LD_PRELOAD
export INPUT_DIR=/path-to-read-datasets/  
export OUTPUT_DIR=/path-to-output-directory/  
export REF_DIR=/path-to-ref-directory/
REF=ref.fasta   
R1=R1.fastq.gz  
R2=R2.fastq.gz  
```
### 6. Create the index files for the reference sequence
```bash  
bash create_reference_index.sh  
```

### 7. Run the pipeline
Note that the script uses default setting for creating multiple MPI ranks based on the system configuration information using hostfile.
```bash
bash run_pipeline_cluster.sh [docker | podman | "sudo docker"]
* docker/podman/"sudo docker" : optional argument. It takes docker by default.
```

# Results

Our latest results are published in this [blog](https://community.intel.com/t5/Blogs/Tech-Innovation/Artificial-Intelligence-AI/Intel-Xeon-is-all-you-need-for-AI-inference-Performance/post/1506083).


# Instructions to run the pipeline on an AWS ec2 instance (Single Node)
The following instructions run seamlessly on a standalone AWS ec2 instance. To run the following steps, create an ec2 instance with Ubuntu-22.04 having at least 60GB of memory and 500GB of disk. The input reference sequence and the paired-ended read datasets must be downloaded and stored on the disk.

### One-time setup
This step takes around ~15 mins to execute. During the installation process, whenever prompted for user input, it is recommended that the user select all default options.
```bash
wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/2.1/Source_code_with_submodules.tar.gz
tar -xzf Source_code_with_submodules.tar.gz
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/scripts/aws
bash deepvariant_ec2_setup.sh
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
bash create_reference_index.sh
```

### Run the pipeline.
Note that the script uses default setting for creating multiple MPI ranks based on the system configuration.
```bash
bash run_pipeline_ec2.sh
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
wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/2.1/Source_code_with_submodules.tar.gz
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
