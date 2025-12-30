# AutoDock Vina Workflow on AWS Batch (Nextflow)

## Overview
This workflow runs **AutoDock Vina docking** using **Nextflow + AWS Batch**, where each complex is executed as an **independent AWS Batch job** inside a Docker container.

Key points:
* 1 complex = 1 Nextflow task = 1 AWS Batch job
* Up to **8 complexes run in parallel** (`maxForks = 8`)
* Work directory is stored in S3 (`workDir`)
* Outputs are published to S3 per complex

**For setting up AWS Batch, follow [these](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/applications/common/nextflow) instructions.**

## Build DockerImage

```bash
cd ../build_docker/
docker build -t autodock_vina_nextflow:latest --build-arg FLAVOR=nextflow .
cd ../nextflow/
```

## Run Command
```bash
nextflow run main.nf \
  -profile awsbatch_vina_test \
  --list "$PWD/complexes_16.txt" \
  --data_dir "$PWD/data" \
  --out_s3 s3://<bucket-name>/outputs-vina-test
```

## Inputs
### 1) Complex list

`complexes_16.txt` contains one complex ID per line:
```bash
1ac8
1gkc
1hnn
1hwi
...
```

### 2) Input directory layout

`data/` must contain one directory per complex:
```bash
data/
  1ac8/
  1gkc/
  1hnn/
  ...
```

Each complex directory must include:
* `protein.pdbqt`
* `rand-0.pdbqt`
* `vina_rand-0_config.txt`

## Execution & Isolation
* Each complex becomes a separate Nextflow task and AWS Batch job.
* Each task runs in a unique Nextflow work sandbox (backed by S3 workDir).
* Each job runs inside its own container filesystem.
* Inputs are copied into a per-task local `workdir/` before running Vina.

## Outputs

Each task produces a folder named after the complex:
`<complex>/`

Published to:
`s3://<bucket-name>/outputs-vina-test/<complex>/`

Example:
```bash
s3://<bucket-name>/outputs-vina-test/
├── 1ac8/
├── 1gkc/
├── 1hnn/
└── 1hwi/
```
## Docker Image
* Image: `docker.io/<user-name>/<image-name>:<tag>`
* Pulled automatically by AWS Batch
* Runs Vina inside the container (`vina --config ... --cpu ${task.cpus}`)
