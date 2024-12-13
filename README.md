[![GitHub Downloads](https://img.shields.io/github/downloads/IntelLabs/Open-Omics-Acceleration-Framework/total?label=GitHub%20Downloads)](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases)
# Open Omics Acceleration Framework
Intel lab's open sourced data science framework for accelerating digital biology

# Introduction
We are in the epoch of digital biology, that is fueled by the convergence of three revolutions: 1) Measurement of biological systems at high resolution resulting in massive multi-modal, multi-scale, unstructured, distributed data, 2) Novel data science (AI and data management) techniques on this data, and 3) Wide-spread cloud use enabling massive compute and public data repositories, large collaborative projects and consortia. It will require computing and data management at unprecedented scale and speed. However, performance alone would not suffice if it significantly compromised the productivity of biologists and data scientists who are at the forefront of this transformation. 

With a goal to build a performant, cost effective and productive platform, we are building **Open Omics acceleration framework**: a one-click, containerized, customizable, open-sourced framework for accelerating digital biology research. It provides tools and pipelines in the field of genomics, transcriptomics, proteomics, drug molecule search and De novo drug design. The framework is being built with a modular design that keeps in mind the different ways the users would want to interact with it. As shown in the following block diagram, it consists of three layers:
* **Pipeline layer**: for users who are looking for one click solution to run standard pipelines. The pipelines can be accessed in the 'pipelines' subfolder. It provides instrcutions to build & run the docker images. Currently, we support the following pipelines:
  * [**fq2sortedbam**](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/fq2sortedbam): Given gzipped fastq files of an individual, this workflow performs sequence mapping ([BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)) and sorting ([SAMtools](https://github.com/samtools/samtools) sort) to output the sorted BAM file.
  * [**DeepVariant based germline pipeline for variant calling (fq2vcf)**](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/deepvariant-based-germline-variant-calling-fq2vcf): Given paired end gzipped fastq files of an individual, this workflow performs sequence mapping ([BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)), sorting ([SAMtools](https://github.com/samtools/samtools) sort) and variant calling ([Open Omics DeepVariant](https://github.com/IntelLabs/open-omics-deepvariant)) to call the variants in the genome of the individual.
  * [**AlphaFold2-based protein folding**](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/alphafold2-based-protein-folding): Given one or more protein sequences, this workflow performs preprocessing (database search and multiple sequence alignment using Open Omics [HMMER](https://github.com/IntelLabs/hmmer) and [HH-suite](https://github.com/IntelLabs/hh-suite)) and structure prediction ([Open Omics AlphaFold2](https://github.com/IntelLabs/open-omics-alphafold)) to output the structure(s) of the protein sequences. It has support for both AlphaFold2 monomer and AlphaFold2 multimer.
  * [**Single cell RNASeq analysis**](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/single-cell-RNA-seq-analysis): Given a cell by gene matrix, this [scanpy](https://github.com/scverse/scanpy) based workflow performs data preprocessing (filter, linear regression and normalization), dimensionality reduction (PCA), clustering (Louvain/Leiden/kmeans) to cluster the cells into different cell types and visualize those clusters (UMAP/t-SNE).
* **Toolkit layer**: for users who want to use individual tools or to create their own custom pipelines by combining various tools. The toolkit layer can be accessed in the 'applications' subfolder. For each tool, we provide instructions to build and run it. Currently, the tools supported include: genomics (BWA-MEM, minimap2, bcftools, SAMtools, DeepVariant), transcriptomics (STAR aligner), protein folding (AlphaFold2, ESMFold), protein structure and sequence design (RFDiffusion, ProteinMPNN, LM-design, ESM2-inv, ProtGPT2, ESM2 embeddings), molecular docking (AutoDock, AutoDock-Vina), De novo molecule generation (MoFlow).
* **Building blocks layer**: for tool developers, this layer consists of key building blocks -- biology specific and generic AI algorithms and data structures -- that can replace ones used in existing tools to accelerate them or can be used as ingredients to build new efficient tools. This layer can be accessed in the 'lib' subfolder.

<p align="center">
<img src="https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/blob/main/images/Open-Omics-Acceleration-Framework-v3.0.jpg" height="300"/a></br>
</p> 

With a goal of providing a one-stop platform, this framework brings our following repositories for digital biology under one umbrella:
* Architecture efficient versions of several popular applications as part of our toolkit layer (under 'applications' folder)

| Original Application | Our architecure-efficient version |
| --- | --- |
| Short read sequence mapping tool, [BWA-MEM](https://github.com/lh3/bwa) | [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) |
| long read sequence mapping tool, [minimap2](https://github.com/lh3/minimap2) | [mm2-fast](https://github.com/bwa-mem2/mm2-fast) |
| Deep learning based variant calling tool, [DeepVariant](https://github.com/google/deepvariant) | [Open-Omics-DeepVariant](https://github.com/IntelLabs/open-omics-deepvariant) |
| Deep learning based tool for protein structure prediction, [AlphaFold2](https://github.com/google-deepmind/alphafold) | [Open-Omics-AlphaFold](https://github.com/IntelLabs/open-omics-alphafold) |
| Tool for biological sequence analysis using profile HMMs, [HMMER](https://github.com/EddyRivasLab/hmmer) | [IntelLabs HMMER](https://github.com/IntelLabs/hmmer) |
| Tool for HMM based sensitive protein sequence searching, [HH-suite](https://github.com/soedinglab/hh-suite) | [IntelLabs HH-suite](https://github.com/IntelLabs/hh-suite) |

* [Trans-Omics Acceleration Library](https://github.com/IntelLabs/Trans-Omics-Acceleration-Library): As part of our building blocks layer (under 'lib' folder), this is a library containing architecture-efficient versions of key algorithms and data structures used for Omics analysis. 

In addition, we also use several existing AI libraries: oneDNN, oneDAL, oneCCL, Katana Graph, LIBXSMM.

# Getting Started
```sh
# Download release
wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz 
tar -xzf Source_code_with_submodules.tar.gz

# Clone master
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework

# Go to the pipelines directory
cd pipelines
# For running a specific pipeline, follow the instructions in the respective pipeline's README file.

# Go to the directory with toolkit
cd applications

# Go to the directory with biology building blocks to access Trans-Omics Acceleration Library
cd lib/tal

```


# Blogs & Related News
* [Intel Open Omics Acceleration Framework on AWS: fast, cost-efficient, and seamless](https://aws.amazon.com/blogs/hpc/intel-open-omics-acceleration-framework-on-aws-fast-cost-efficient-and-seamless/). AWS HPC blog, April, 2024.
* [Intel Xeon is all you need for AI inference: Performance Leadership on Real World Applications](https://community.intel.com/t5/Blogs/Tech-Innovation/Artificial-Intelligence-AI/Intel-Xeon-is-all-you-need-for-AI-inference-Performance/post/1506083). Blog under Intel Communities/Blogs/Tech Innovation/Artificial Intelligence (AI); July, 2023.
* [Intel and Mila Join Forces for Responsible AI](https://www.intel.com/content/www/us/en/newsroom/news/intel-mila-join-forces-for-responsible-ai.html#gs.ht5v6q). Intel newsroom, September, 2022.
* [Accelerating Genomics Pipelines Using Intel’s Open Omics Acceleration Framework on AWS](https://aws.amazon.com/blogs/hpc/accelerating-genomics-pipelines-using-intel-open-omics-on-aws/). AWS HPC blog, Aug, 2022.
* [Intel Labs Accelerates Single-cell RNA-Seq Analysis](https://community.intel.com/t5/Blogs/Tech-Innovation/Artificial-Intelligence-AI/Intel-Labs-Accelerates-Single-cell-RNA-Seq-Analysis/post/1390715).  Blog under Intel Communities/Blogs/Tech Innovation/Artificial Intelligence (AI); June, 2022.
* [Intel and MILA Join Forces to Put AI to Work in Medical Research](https://www.hpcwire.com/off-the-wire/intel-and-mila-join-forces-to-put-ai-to-work-in-medical-research/). HPCwire, April, 2021.

# Publications

* **GenDP: A Framework of Dynamic Programming Acceleration for Genome Sequencing Analysis**. Yufeng Gu, Arun Subramaniyan, Tim Dunn, Alireza Khadem, Kuan-Yu Chen, Somnath Paul, Md Vasimuddin, Sanchit Misra, David Blaauw, Satish Narayanasamy, Reetuparna Das. Proceedings of the 50th Annual International Symposium on Computer Architecture (ISCA); June, 2023. https://dl.acm.org/doi/abs/10.1145/3579371.3589060.
* **Accelerating Barnes-Hut t-SNE Algorithm by Efficient Parallelization on Multi-Core CPUs**. Narendra Chaudhary,  Alexander Pivovar, Pavel Yakovlev, Andrey Gorshkov and Sanchit Misra. arXiv preprint arXiv:2212.11506; Dec, 2022; 
doi: https://doi.org/10.48550/arXiv.2212.11506.
* **Accelerating Deep Learning based Identification of Chromatin Accessibility from noisy ATAC-seq Data**. Narendra Chaudhary, Sanchit Misra, Dhiraj Kalamkar, Alexander Heinecke, Evangelos Georganas, Barukh Ziv, Menachem Adelman and Bharat Kaul. 21st IEEE International Workshop on High Performance Computational Biology (HiCOMB) May 30, 2022. https://ieeexplore.ieee.org/abstract/document/9835674
* **Accelerating minimap2 for long-read sequencing applications on modern CPUs**. Saurabh Kalikar, Chirag Jain, Md Vasimuddin, Sanchit Misra. Nature Computational Science 2 (2), 78-83, Feb, 2022. https://rdcu.be/cHVAK.
* **GenomicsBench: A Benchmark Suite for Genomics**. Arun Subramaniyan, Yufeng Gu, Timothy Dunn, Somnath Paul, Md Vasimuddin, Sanchit Misra, David Blaauw, Satish Narayanasamy, Reetuparna Das. IEEE International Symposium on Performance Analysis of Systems and Software (ISPASS), 2021.https://ieeexplore.ieee.org/document/9408208.
* **LISA: Learned indexes for sequence analysis**. Darryl Ho, Saurabh Kalikar, Sanchit Misra, Jialin Ding, Vasimuddin Md, Nesime Tatbul, Heng Li, Tim Kraska. bioRxiv 2020.12.22.423964; doi: https://doi.org/10.1101/2020.12.22.423964.
* **Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems**. Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. https://ieeexplore.ieee.org/document/8820962.
* **Performance extraction and suitability analysis of multi- and many-core architectures for next generation sequencing secondary analysis**. Sanchit Misra, Tony Pan, Kanak Mahadik, George Powley, Priya N Vaidya, Md Vasimuddin, Srinivas Aluru. International Conference on Parallel Architectures and Compilation Techniques (PACT), 2018. https://dl.acm.org/doi/abs/10.1145/3243176.3243197.
* **Identification of Significant Computational Building Blocks through Comprehensive Deep Dive of NGS Secondary Analysis Methods**. Md Vasimuddin, Sanchit Misra, Srinivas Aluru.  BioRxiv 2018 301903. https://www.biorxiv.org/content/10.1101/301903v3.abstract.


# Privacy Notice

Please note that Intel does not collect or process any personal data of end users for Open Omics Acceleration Framework. All the inputs and outputs are limited to the end users's local system and Intel does not have access to the data. For more information on Intel's Global Privacy Notice, please click here: https://www.intel.com/content/www/us/en/privacy/intel-privacy-notice.html.
