# GROMACS Nextflow + AWS Batch Workflow

This document explains **how the GROMACS Nextflow workflow works**, with a focus on **task isolation, execution flow**, and **output** organization when running on **AWS Batch**.

## Overview
This workflow runs multiple GROMACS protein simulations in parallel using:
* Nextflow DSL2 for orchestration
* AWS Batch as the executor
* Docker container with GROMACS pre-installed
* S3 for work directories and final outputs

**For setting up AWS Batch, follow [these](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/applications/common/nextflow) instructions.**

Each protein (PDB file) is executed as an **independent AWS Batch job**, fully isolated from all others.

## Input Structure
### 1. Protein list

`pdb_list.txt`
```bash
1abc.pdb
2xyz.pdb
3def.pdb
...
```

Each line corresponds to **one independent GROMACS run**.
### 2. Input directory
```bash
pdbs/
├── 1abc.pdb
├── 2xyz.pdb
├── 3def.pdb
```

### 3. MDP files
```bash
mdtut_minim.mdp
mdtut_nvt.mdp
mdtut_npt.mdp
mdtut_ions.mdp
mdtut_md.mdp
```
These files are staged into each task independently.

## High-Level Execution Flow

For **each protein** in `pdb_list.txt`:
* Nextflow creates one task
* One AWS Batch job is submitted
* One Docker container is started
* GROMACS runs inside the container
* Outputs are copied to S3 in an isolated directory

## Task Isolation Explained
### 1. Nextflow Work Directory (S3)

workDir = 's3://<bucket-name>/<folder-name>'

For each task, Nextflow creates a unique hash-based directory:
```bash
work-gromacs/
├── aa/bbccccdd1234/   # Protein A
├── ee/ff0011223344/   # Protein B
└── ...
```

No two proteins share the same work directory.

### 2. Container-Level Isolation

Each task runs in a separate AWS Batch container.

Inside the container:
* `/input` is used for inputs
* `/output` is used for GROMACS results

Even though these paths have the same names across jobs, **containers do not share filesystems**.
`/input` and `/output` are isolated per protein.

### 3. Input Staging

Inside the task script:
```bash
cp "${pdb_file}" /input/
cp mdtut_*.mdp /input/
cp custom_script.sh /input/
```

Each container receives its **own copy** of:
* PDB file
* MDP files
* Custom workflow script

### 4. GROMACS Execution

Execution happens via the container entrypoint:
```bash
/entrypoint.sh <protein.pdb>
```

The entrypoint:

* Activates the GROMACS environment
* Runs the full EM → NVT → NPT → MD workflow
* Writes results to `/output/output_<timestamp>/`

### 5. Output Collection (Per Protein)

After GROMACS finishes:
```bash
mkdir -p "${prefix}"
cp -r /output/output_*/* "${prefix}/"
```

This creates a **protein-specific output directory** inside the task workdir:

```bash
<task-workdir>/
└── 1abc/
    ├── em.gro
    ├── nvt.gro
    ├── md01.xtc
    ├── md01.tpr
    └── ...
```

### 6. Publishing Outputs to S3

Nextflow publishes outputs using:
```bash
publishDir params.out_dir, mode:'copy'
```

Resulting S3 layout:
```bash
s3://<bucket-name>/gromacs-outputs/
├── 1abc/
├── 2xyz/
├── 3def/
```

Each protein’s results are stored in a separate S3 directory.


## How to Run
```bash
nextflow run main.nf -profile awsbatch_single_node --pdb_list "$PWD/pdb_list.txt" --pdb_dir  "$PWD/pdbs" --custom_script "$PWD/run_commands.sh" --out_dir   s3://<bucket-name>/gromacs-outputs
```

