# fq2vcf: OpenOmics Deepvariant based Variant Calling Pipeline  
### Overview:  
OpenOmics' fq2vcf is a highly optimized, distributed, deep learning-based short-read germline variant calling pipeline for x86 CPUs. 
The pipeline comprises of:   
1. Sequence alignment using: bwa-mem2 (a highly optimized version of bwa-mem) for short reads or mm2-fast  (a highly optimized version of minimap2) for long reads 
2. SortSAM using SAMTools  
3. An optimized version of DeepVariant tool for Variant Calling
   
The following figure illustrates the pipeline:

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
- Alignment tools (bwa-mem2/mm2-fast) can be selected using parameter _--read_type short_ for bwa-mem2 or _--read_type long_ for mm2-fast to run_fq2bams.py  
- All the alignment tools' command-line parameters are supported and can be provided as --params "'<param_list>'"  to run_fq2bams.py
- Inputs reads are expected in gz file format
- refdir is expected to contain the bwa-mem2 genome and index. You can also index the reference genome during the run by enabling "--rindex" to fq2bams command-line (for bwa-mem2). For mm2-fast, it creates the reference genome index during the alignment step on-the-fly.  

```bash
docker run  --volume <refdir>:/refdir <readsdir>:/readsdir <outdir_fq2bams>:/outdir fq2bams:latest python run_fq2bams.py --ref /refdir/<reference_file> --reads  /readsdir/<read1>  /readsdir/<read2>  --output /outdir/<outBAMfile>   

docker run  --volume <refdir>:/refdir <outdir_fq2bams>:/indir <output>:/outdir  bams2vcf:latest python run_bams2vcf.py --ref /refdir/<reference_file> --input /indir/  --output /outdir/<outVCFfile>   
```

Example using mm2-fast for aligment:  
```bash
docker run  --volume xpath:/refdir ypath:/readsdir outpath:/outdir fq2bams:latest python run_fq2bams.py --read_type long --ref /refdir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --reads  /readsdir/HG002_35x_PacBio_14kb-15kb.fastq.gz  --output /outdir/   

docker run  --volume xpath:/refdir outpath:/indir outpath:/outdir  bams2vcf:latest python run_bams2vcf.py --ref /refdir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --input /indir/  --output /outdir/hg002_pacbio_vcfile      
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
