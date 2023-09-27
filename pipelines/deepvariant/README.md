# OpenOmics Deepvariant Pipeline
### OpenOmics Deepvariant Pipeline is a highly optimized, distirbuted, deep-learning-based short-read germline variant calling pipeline on x86 CPU clusters. The pipeline comprises of: 1. bwa-mem2 (a highly optimized version of bwa-mem) for sequence mapping, 2. distributed SAM sorting using samtools, and 3. an optimized version of DeepVariant tools for variant calling.

### 0. Notes:
* The source code of bwa-mem2, samtools, and DeepVariant are residing in:
```Open-Omics-Acceleration-Framework/applications/ ```.
* We build bwa-mem2 and samtools from source; while for DeepVariant, we build a docker image and then use the built image while running the pipeline. Note that, the pre-built image is not available on dockerhub and the image needs to be built from source.
* We provide scripts for setting-up miniconda environment (setup_env.sh), compilation of the required pipeline tools, builing & loading of DeepVariant docker image (setup.py). These scripts located at _Open-Omics-Acceleration-Framework/pipelines/deepvariant_. Note that, all the script by default are written for docker.  

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

### 1. Clone the repo:  
```bash
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
```

### 2. Setting Envionment
```bash
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant/scripts/cluster/
#Tested with Ubuntu 22.04.2 LTS
source setup_env.sh  dv_env # Setting environment with name dv_env.
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
To run docker without sudo acess follow [link](https://docs.docker.com/engine/install/linux-postinstall/), or use "_sudo docker_" as an argument in the below script instead of docker.
```bash
# If you are using single node, comment out line No. 56, i.e., "bash load_deepvariant.sh" of setup.sh.
source setup.sh [docker | podman | "sudo docker"]
* docker/podman/"sudo docker" : optional argument. It takes docker by default.
```
Note: It takes ~30 mins to create the docker image. Docker build might break if you are behind a proxy. Example to provide proxy for building a docker image is shown in the setup.sh file. [Follow](https://docs.docker.com/network/proxy/) instructions for more details.

### 5. Create _config_ file
We need a reference sequence and paired-ended read datasets. Open the "config" file and set the input and output directories as shown in config file. The sample config contains the following lines to be updated.  

export INPUT_DIR=/path-to-reference-sequence-and-read-datasets/  
export OUTPUT_DIR=/path-to-output-directory/  
REF=ref.fasta   
R1=R1.fastq.gz  
R2=R2.fastq.gz  

### 6. Create the index files for the reference sequence
```bash  
bash create_reference_index.sh  
```

### 7. Run the pipeline
Note that the script uses default setting for creating multiple mpi ranks based on the system configuration information using hostfile.
```bash
bash run_pipeline_cluster.sh
```

# Results

For detailed information, please refer to the [blog](https://community.intel.com/t5/Blogs/Tech-Innovation/Artificial-Intelligence-AI/Intel-Xeon-is-all-you-need-for-AI-inference-Performance/post/1506083).


# Instructions to run the pipeline on an AWS ec2 instance
The following instructions run seamlessly on a standalone AWS ec2 instance. To run the following steps, create an ec2 instance with Ubuntu-22.04 having at least 60GB of memory. The input reference sequence and the paired-ended read datasets must be downloaded and stored on the disk.

### One-time setup
This step takes around ~15 mins to execute
```bash
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant/scripts/aws
bash deepvariant_ec2_setup.sh
```

### Modify _config_ file
We need a reference sequence and paired-ended read datasets. Open the "_config_" file and set the input and output directories as shown in config file.
The sample config contains the following lines to be updated.
```bash
export INPUT_DIR=/path-to-reference-sequence-and-read-datasets/
export OUTPUT_DIR=/path-to-output-directory/
REF=ref.fasta
R1=R1.fastq.gz
R2=R2.fastq.gz
```

### Create the index files for the reference sequence
```bash
bash create_reference_index.sh
```

### Run the pipeline.
Note that the script uses default setting for creating multiple mpi ranks based on the system configuration.
```bash
bash run_pipeline_ec2.sh
```
<!--

# Instructions to run the pipeline on an AWS ParallelCluster
The following instructions run seamlessly on AWS ParallelCluster. To run the following steps, create an [AWS ParallelCluster](https://docs.aws.amazon.com/parallelcluster/latest/ug/install-v3.html) with Ubuntu-22.04 using ParallelCluster configuration file and login into the host node. The input reference sequence and the paired-ended read datasets must be downloaded and stored on the disk in the _/shared_ folder.

### One-time setup
This step takes around ~15 mins to execute
```bash
cd /shared
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant/scripts/aws
bash deepvariant_setup.sh
```
### Modify _config_ file
We need a reference sequence and paired-ended read datasets. Open the "_config_" file and set the input and output directories as shown in config file.
The sample config contains the following lines to be updated.
```bash
export INPUT_DIR=/path-to-reference-sequence-and-read-datasets/
export OUTPUT_DIR=/path-to-output-directory/
REF=ref.fasta
R1=R1.fastq.gz
R2=R2.fastq.gz
```

### Create the index files for the reference sequence
```bash
bash create_reference_index.sh
```

### Allocate compute nodes and install the prerequisites into the compute nodes.
```bash
bash pcluster_compute_node_setup.sh <num_nodes> <allocation_time>
# num_nodes: The number of compute nodes to be used for distributed multi-node execution.
# allocation_time: The maximum allocation time for the compute nodes in "hh:mm:ss" format. The default value is 2 hours, i.e., 02:00:00.
Example command for allocating 4 nodes for 3 hours -
bash pcluster_compute_node_setup.sh 4 "03:00:00"
```

### Run the pipeline.
Note that the script uses default setting for creating multiple mpi ranks based on the system configuration.
```bash
bash run_pipeline_pcluster.sh
```
-->
