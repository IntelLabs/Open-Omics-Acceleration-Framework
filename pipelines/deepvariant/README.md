# Deepvariant Pipeline
### This repo presents OpenOmics DeepVariant Pipeline: A highly optimized scalable Deep Learning based short-read vaariant calling on x86 CPU clusters. The pipeline comprises of highly optimized: 1. bwa-mem2 for sequence mapping, 2. samtools for sorting output of bwa-mem2, and 3. DeepVariant for variant calling using sorted BAM records out of sorting.

0. Pipeline tools' location:   
bwa-mem2,samtools, deepvariant C/C++ based tools residing in:
```Open-Omics-Acceleration-Framework/applications/ ```.
While DeepVariant tools is a docker image; it is not present in the repo as of now; the user needs to create image, convert to tar file and distribute across the cluste nodes. 
   Prerequisite : Docker /Podman  
   All the script by default supports podman. If you are using docker use:  **alias podman=docker**
```
cd Open-Omics-Acceleration-Framework/applications/deepvariant
podman build .
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant/
podman save -o deepvariant.tar "IMAGE ID"
```
1. Download data and DeepVariant image:
```
export INPUT_DIR=./    ## temp
export OUTPUT_DIR=./   ## temp
#save image to tar file if you are using multiple nodes.
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant/
podman save -o deepvariant.tar "IMAGE ID"
```

2. Download code
```
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git  
cd Open-Omics-Acceleration-Framework/pipelines/deepvariant/  

source setup_env.sh  ## creates new_env conda environment
conda activate new_env   ## activate conda env
```
3. Compute nodes setup:  
3.1.  Compute cluster using slurm job scheduler, e.g.: Allocate 2 spr dnp8480 node having 112 cores per compute node for 12 hrs.
```
salloc --ntasks=2 --partition=nextgenq --constraint=dnp8480 --cpus-per-task=112 --time=12:0:0
srun hostname > hostfile  
```  

3.2 Standalone machine
```
hostname > hostfile
```
3.3 Aws cluster setup

[Follow Instructions](AWS_CLUSTER_SETUP.md) for detailed instructions.


4. Loads docker/podman images in the allocated nodes.
```
source setup.sh      ## tested with gcc 8.5.0
```
Note: It takes ~10 mins to load the image in all the hostfile nodes. 
You can check if the image in loaded in the compute node using:
```
$ podman images
```
The output should look like this if the image is loaded.  
REPOSITORY             TAG         IMAGE ID      CREATED      SIZE   
localhost/deepvariant  latest      9eaf947f7e2e  6 weeks ago  19.8 GB   

5. Run the following script after the image is loaded on all the compute nodes listed in hostfile.  
Usage: source run_piplin.sh <#ranks> <#threads> <#threads> <#shards> <#ppn>  
* ranks: Number of mpi process that we want the pipeline to run on  
* threads/shards: parameters to different tools in the pipeline, calculated as below  
* ppn: mpi process per compute node. If we are running 2 ranks and provide ppn as 2, then both the ranks will run on 1 compute code (assuming dual socket machine, it will run 1 rank per socket)  
Note: for best performance, we advice to run 4 ranks per socket on spr nodes. So, assuming dual-socket compute node you can run 8 ranks on 1 compute node.  
```
ppn=8  
CPUS=$(lscpu | grep -E '^CPU\(s\)' | awk  '{print $2}')
Cores=$(lscpu| grep -E '^Core\(s\)' | awk  '{print $4}')
Thread=$(lscpu | grep -E '^Thread' | awk  '{print $4}')

a=$(( $(( ${Cores}*${Thread}*2 / $ppn )) - 4 )) 
b=$(( $(( ${Cores}*${Thread} )) / $ppn ))

sh run_pipeline.sh 8 $a $a $b $ppn
```
Running on the test dataset on spr should take around ~10 mins.  

6. Results  
