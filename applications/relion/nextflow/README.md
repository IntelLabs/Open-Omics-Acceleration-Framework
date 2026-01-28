# RELION 2D Classification on AWS Batch (Nextflow)
## Overview
This setup runs **RELION 2D classification** using **Nextflow + AWS Batch (EC2)**, where **each dataset runs as an isolated Batch job on a c7i.8xlarge instance**.

The final working solution explicitly uses:
* default ECS-optimized AMI
* EC2 Launch Template
* AWS Batch Managed Compute Environment

**For setting up AWS Batch, follow [these](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/applications/common/nextflow) instructions.**


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

## Build the RELION container and push to ECR

### 1) Build the Docker image locally:

From the folder containing the Dockerfile:

```bash
cd ../build_docker/
docker build -t relion_common_scan:latest .
cd ../nextflow/
```

### 2) Create the ECR repository (one-time):

Create the ECR repo (skip if it already exists):
```bash
aws ecr create-repository --region us-east-2 --repository-name relion_common_scan
```
### 3) Authenticate Docker to ECR:
```bash
aws ecr get-login-password --region us-east-2 \
| docker login --username AWS --password-stdin <ACCOUNT_ID>.dkr.ecr.us-east-2.amazonaws.com

```

### 4) Tag the image for ECR:
```bash
docker tag relion_common_scan:latest <ACCOUNT_ID>.dkr.ecr.us-east-2.amazonaws.com/relion_common_scan:latest
```
### 5) Push to ECR:

```bash
docker push <ACCOUNT_ID>.dkr.ecr.us-east-2.amazonaws.com/relion_common_scan:latest
```
## Nextflow Execution (AWS Batch)
     
### Step 1: Create Launch Template(with 200gb disk volume)
```bash
aws ec2 create-launch-template \
  --region us-east-2 \
  --launch-template-name relion-ecs-200g \
  --launch-template-data '{
    "BlockDeviceMappings": [
      {
        "DeviceName": "/dev/xvda",
        "Ebs": {
          "VolumeSize": 200,
          "VolumeType": "gp3",
          "DeleteOnTermination": true
        }
      }
    ]
  }'
```

### Step 2: Create NEW Compute Environment (with 200 GB disk)

```bash
aws batch create-compute-environment \
  --region us-east-2 \
  --compute-environment-name relion-ce-200g \
  --type MANAGED \
  --state ENABLED \
  --service-role arn:aws:iam::<ACCOUNT_ID>:role/AWSBatchServiceRole \
  --compute-resources '{
    "type":"EC2",
    "minvCpus":0,
    "desiredvCpus":0,
    "maxvCpus":256,
    "instanceTypes":["c7i.8xlarge"],
    "subnets":[
      "subnet-<SUBNET-ID>",
      "subnet-<SUBNET-ID>"
    ],
    "securityGroupIds":["<SID>"],
    "instanceRole":"ecsInstanceRole",
    "launchTemplate":{
      "launchTemplateName":"relion-ecs-200g",
      "version":"$Latest"
    }
  }'
```

### Step 3: Create NEW Job Queue
```bash
aws batch create-job-queue \
  --region us-east-2 \
  --job-queue-name relion-queue-200g \
  --priority 1 \
  --compute-environment-order order=1,computeEnvironment=relion-ce-200g
```

### Step 4: Run Nextflow

Update Nextflow config to use:

queue = `relion-queue-200g`

```bash                            
nextflow run main.nf   -profile awsbatch_relion_2d   --relion_image $ACCOUNT_ID.dkr.ecr.us-east-2.amazonaws.com/relion_common_scan:latest   -ansi-log false
```
