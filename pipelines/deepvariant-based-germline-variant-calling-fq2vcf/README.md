# fq2vcf: OpenOmics Deepvariant based Variant Calling Pipeline  
### Overview:  
OpenOmics' fq2vcf is a highly optimized, distributed, deep learning-based short-read germline variant calling pipeline for x86 CPUs. 
The pipeline comprises of:   
1. bwa-mem2 (a highly optimized version of bwa-mem) for sequence mapping  
2. SortSAM using samtools  
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
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy -t bams2vcf -f Dockerfile_bamsvcf  .   
```

### 3. Run the Dockers  
Notes:  
<refdir> is expected to contain the bwa-mem2 index. You can index the reference during the run by enabling "--rindex" to fq2bams commandline.  

```bash
docker run  --volume <refdir>:/refdir <readsdir>:/readsdir <outdir_fq2bams>:/outdir fq2bams:latest python run_fq2bams.py --ref /refdir/<reference_file> --reads  /readsdir/<read1>  /readsdir/<read2>  --output /outdir/<outBAMfile>   

docker run  --volume <refdir>:/refdir <outdir_fq2bams>:/indir <output>:/outdir  bams2vcf:latest python run_bams2vcf.py --ref /refdir/$ref --input /indir/  --output /workdir/<outVCFfile>   
```

# Results

Our latest results are published in this [blog](https://community.intel.com/t5/Blogs/Tech-Innovation/Artificial-Intelligence-AI/Intel-Xeon-is-all-you-need-for-AI-inference-Performance/post/1506083).

