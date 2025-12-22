# AutoDock CPU Workflow on AWS Batch (Nextflow)

## Overview

This workflow runs AutoDock CPU docking using Nextflow + AWS Batch, where each protein (complex) is executed as an independent AWS Batch job inside a Docker container pulled from Docker Hub.

Key characteristics:
* 1 protein = 1 Nextflow task = 1 AWS Batch job
* Up to 8 tasks run in parallel
* Each task runs on a separate c7i.8xlarge EC2 instance
* Inputs and outputs are fully isolated per task
* Outputs are published to Amazon S3

## Workflow Execution

The workflow is launched using:

```bash
nextflow run main.nf   -profile awsbatch_single_node   --protein_list protein_list.txt   --input_dir inputs   --out_dir s3://srilekha-new-s3/autodock-outputs   -ansi-log false
```

## Input Structure
### 1. Protein List

`protein_list.txt` contains one protein ID per line:

```bash
10gs
1a30
1ac8
...
```

Each entry corresponds to one independent docking job.

### 2. Input Directory Layout

The inputs/ directory must contain one subdirectory per protein, named exactly as listed in `protein_list.txt`.

Example:
```bash
ls inputs/
10gs  1a30  1gpk  1h23  1hwi  1jla  1lbk  1n1m  1opk  1pmn  1r58  1sq5  1tt1  1unl 
```
Each protein directory must include:
* `protein.maps.fld`
* `rand-0.pdbqt (ligand file)`
* Any other AutoDock-required map files

## Execution Model & Isolation
### Task-Level Isolation
* Each protein is processed as one Nextflow task
* Each task is submitted as a separate AWS Batch job
* AWS Batch schedules one job per EC2 instance
* Each job runs inside its own Docker container with an isolated filesystem
## Parallelism
* Maximum 8 tasks run concurrently
* Controlled by:
   * executor.queueSize = 8
   * process.maxForks = 8
   
Remaining proteins are queued automatically by AWS Batch.

## Output Structure:

Outputs are published to S3 under a **run-specific directory**:
`s3://srilekha-new-s3/autodock-outputs/<run_id>/`

Each protein produces its own output folder:
`<protein_id>_job<index>/`

Example:

```bash
s3://srilekha-new-s3/autodock-outputs/20250118_153210_curious-hopper/
├── 10gs_job0/
├── 1a30_job1/
├── 1gpk_job2/
├── 1h23_job3/
└── 1hwi_job4/
```

Each folder contains:
* AutoDock output files (`.dlg`,`.pdbqt`, `logs`, `maps`, etc.)
* Files generated only for that specific protein

## Docker Image
* Image: `docker.io/809314/autodock_64wi:latest`
* Pulled automatically by AWS Batch
* Contains AutoDock CPU binary (`autodock_cpu_64wi`) and runtime dependencies

**To collect workflow outputs from the S3 bucket**:
```bash
aws s3 cp s3://<bucket-name>/autodock-outputs/ . --recursive
```
