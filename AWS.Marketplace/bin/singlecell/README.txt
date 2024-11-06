# Scanpy: OpenOmics' accelerated Single Cell RNASeq analysis pipeline
## Overview
Given a cell by gene matrix, this scanpy based pipeline performs data preprocessing (filter, linear regression and normalization), dimensionality reduction (PCA), clustering (Louvain/Leiden/kmeans) to cluster the cells into different cell types and visualize those clusters (UMAP/t-SNE). 

The pipeline is present in the form of Jupyter Notebook.

## Example run with complete 1.3M mouse brain cell data
    ``` 
    cd ~/bin/singlecell/
    bash scanpy.sh      
    ```
    Note: This script will run the Jupyter Notebook server, you can connect to the server (through the IP of this instance in your browser) to run the scanpy pipeline using the notebook.
    Runtime on m7i.48xlarge w/ default AMI hardware: < 600 sec

# Source code 
```https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/single-cell-RNA-seq-analysis```

## Run with your own input data
Currently, the pipeline only supports the analysis of 1.3M mouse brain single cell data

