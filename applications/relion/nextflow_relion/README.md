# RELION 2D Classification on AWS Batch (Nextflow)
## Overview
This setup runs **RELION 2D classification** using **Nextflow + AWS Batch (EC2)**, where **each dataset runs as an isolated Batch job on a c7i.8xlarge instance**.

The final working solution explicitly uses:
* ECS-optimized AMI (Amazon Linux 2023)
* EC2 Launch Template
* AWS Batch Managed Compute Environment

## What the workflow does

Input datasets live under:
```bash
relion_data/
  ├── data_1/
  ├── data_2/
  ├── ...
```

* Nextflow discovers each `data_*` directory
* Each dataset is submitted as a **separate AWS Batch job**
* Each job:
   * Runs on **one EC2 instance** (c7i.8xlarge)
   * Executes `mpirun relion_refine_mpi` inside a container
   * Writes output under:
   ```bash
   Class2D/run_<dataset_id>
   ```
* Outputs are published to S3

* Jobs are fully isolated: no shared working directories, no cross-dataset interference.

## Why AMI matters

AWS Batch on EC2 runs containers via **ECS**.

For ECS to start containers correctly, the EC2 instance must have:
* ECS agent installed
* A working container runtime (Docker or containerd)
* Correct system configuration
This is guaranteed **only** by an **ECS-optimized AMI**.

## What went wrong initially
* Default Batch compute environment(default AMI and 30 GB root disk):
   * ECS agent present
   * Container runtime not usable (dockerVersion = null)
* Result:
   * Batch jobs stuck in STARTING
   * No container logs
   * No visible runtime errors

When the root disk size alone was increased to 200 GB (keeping the default AMI):
* Jobs were able to move from STARTING → RUNNING
* RELION executions completed successfully

But to ensure a consistently healthy and supported runtime, the setup was switched to an ECS-optimized Amazon Linux 2023 (AL2023) AMI, after which jobs ran reliably.

## Why Launch Template was required

When no launch template is provided, AWS Batch uses defaults:
* Default ECS-optimized AMI (region-dependent)
* Default root disk size 

This is often sufficient for **small images**, but **not** for RELION.

## RELION-specific requirements
* Container image size ≈ 30 GB
* Docker/containerd needs 2–3× image size during pull & unpack

## Without a launch template:
* Disk fills up
* Container runtime fails silently
* ECS agent cannot start tasks
## Launch Template fixed this by
* Forcing a known-good ECS-optimized AL2023 AMI
* Increasing root volume to 200 GB (gp3)
* Making instance creation deterministic and reproducible

## Why Amazon Linux 2023 (AL2023)

Although Amazon Linux 2 is supported, AL2023 was chosen because:
* Uses containerd instead of legacy Docker
* More robust for:
   * Large images
   * Modern instance types (c7i)
   * HPC / MPI workloads

## Final Architecture:
```bash
Nextflow
↓
AWS Batch Job Queue
↓
AWS Batch Compute Environment (EC2)
↓
Launch Template
├── ECS-optimized AL2023 AMI
├── 200 GB root disk
└── ecsInstanceRole
↓
ECS Agent + containerd
↓
RELION container (mpirun relion_refine_mpi)
```

Command to run:
```bash
nextflow run main.nf -profile awsbatch_relion_2d --data_root relion_data --dataset_glob 'data_*' --star_file particles.star --out_dir s3://srilekha-new-s3/relion-outputs  -ansi-log false
```
