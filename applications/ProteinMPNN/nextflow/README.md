# ProteinMPNN Workflow on AWS Batch (Nextflow)

## Overview
This workflow runs **ProteinMPNN** using **Nextflow + AWS Batch**, where each complex is executed as an **independent AWS Batch job** inside a Docker container.

Key points:
* 1 complex = 1 Nextflow task = 1 AWS Batch job
* Up to **8 complexes run in parallel** (`maxForks = 8`)
* Work directory is stored in S3 (`workDir`)
* Outputs are published to S3 per complex

**For setting up AWS Batch, follow [these](https://github.com/intel-sandbox/TransOmics.OpenOmicsInternal/tree/main/applications/common/nextflow) instructions.**

## Build DockerImage

```bash
cd ../
docker build -t pmpnn:latest .
cd ../nextflow/
```

## Run Command
```bash
nextflow run main.nf \
  -profile awsbatch_single_node \
  --input_dir "<input_pdbs>" \
  --tsv pmpnn_input.tsv \
  --out_s3 s3://<bucket-name>/outputs-pmpnn-test
```

## Inputs
### 1) Design list (TSV)
`pmpnn_input.tsv` defines the ProteinMPNN design jobs.
Each row corresponds to one ProteinMPNN run.

```bash
5TPN.pdb	
1A2K.pdb	
3GFT.pdb
2JQH.pdb	
6M0J.pdb	
4N6H.pdb
2RH1.pdb
5Y3Q.pdb
...
```

### 2) Input directory layout

`input_dir/` must contain the PDB structures referenced in the TSV.:
```bash
<input_pdbs>/
  5TPN.pdb	
  1A2K.pdb	
  3GFT.pdb
  2JQH.pdb	
  6M0J.pdb	
  4N6H.pdb
  2RH1.pdb
  5Y3Q.pdb
  ...
```

## Execution & Isolation
* Each row in the TSV becomes a separate Nextflow task and AWS Batch job
* Each task runs in its own Nextflow work directory (backed by S3 workDir).
* Each job runs inside an isolated container filesystem
* The required PDB file and parameters are copied into a per-task local workdir/

## Outputs

Each task produces a folder named after the input PDB (without extension)::
`<pdb_id>/`

Published to:
`s3://<bucket-name>/outputs-pmpnn-test/<pdb_id>/`

Example:
```bash
s3://<bucket-name>/outputs-pmpnn-test/<pdb_id>/
├── 5TPN/
    ├── seq/
        ├──5TPN.fa
```
## Docker Image
* Image: `docker.io/<user-name>/<image-name>:<tag>`
* Pulled automatically by AWS Batch
* Runs ProteinMPNN inside the container 
* Uses ${task.cpus} to control CPU usage during sequence sampling
