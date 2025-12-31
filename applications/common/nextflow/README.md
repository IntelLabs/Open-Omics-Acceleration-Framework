# Nextflow + AWS Batch + Docker – Complete Setup & Execution Guide

This guide helps you run large-scale parallel bioinformatics workloads using:

- **Nextflow** → Workflow engine  
- **AWS Batch** → Scalable compute backend  
- **Docker** → Reproducible environment  
- **Amazon S3** → Storage for work and outputs  
- **EC2** → Compute nodes launched automatically  

> ⚠️ Region requirement: **United States (Ohio) – us-east-2**

---

## 1) AWS Environment Setup

---

### 1.1 Create an IAM User

*(If your IAM user already exists: go to 1.1.3: Generate Access Keys)*

* Go to AWS Console → IAM → Users
* Click **Create user**
* Enter username
* Click **Next**

---

### 1.2 Attach Required Permissions

Choose **Attach policies directly** and attach:

- `AmazonEC2ContainerRegistryFullAccess`
- `AmazonS3FullAccess`
- `AWSBatchFullAccess`

Click **Next** → **Create user**

---

### 1.3 Generate Access Key

1. Click your IAM user
2. Click Create access key
3. **Choose Command Line Interface (CLI)**
4. Tick confirmation
5. Skip “Set description tag “ - it’s optional.
6. Click on Create Access Key 
7. Copy and save the following:

- Access Key ID
- Secret Access Key
8. Click done

---

**Save the above access key and secret access key as it will be applied during aws configure on the  terminal.**


## 2) Install AWS CLI on your local machine

### A) apt (requires sudo)
```bash
sudo apt update
sudo apt install awscli
aws --version
```
### B) pip (no sudo required)

```bash
pip install --upgrade --user awscli
export PATH=$HOME/.local/bin:$PATH
echo 'export PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
aws --version
```
### C) Official AWS CLI v2 installer

```bash
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip -d $HOME/awscliv2
$HOME/awscliv2/aws/install --install-dir $HOME/.aws-cli --bin-dir $HOME/bin
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
aws --version
```

Type: `aws configure` on your local machine and then enter:

* AWS Access Key ID
* AWS Secret Access Key
* Default region → `us-east-2`
* Output format → `json`

Verify:

Type `aws sts get-caller-identity`:

You should see:
```bash
{
  "UserId": "...",
  "Account": "...",
  "Arn": "arn:aws:iam::<ACCOUNT>:user/<user_name>"
}
```

---
 
## 3) IAM Roles for AWS Batch: AWS Batch requires two IAM roles.

### 3.1 ecsInstanceRole

1. Go to **AWS Console → IAM → Roles**
2. Search for **ecsInstanceRole**
3. Open the role
4. Under **Permissions**, attach the following policies:
   - `AmazonEC2ContainerServiceforEC2Role`
   - `AmazonS3FullAccess`

---

## 4) AWS Batch Setup

### 4.1 Create Compute Environment

Navigate to:

AWS Console → Batch → Environments → Create environment → Compute environment

#### Set the following values:

- **Environment configuration**: EC2
- **Type**: Managed
- **Compute environment name**: `<your-compute-environment-name>`
- **Service Role**: `AWSBatchServiceRole`
- **Instance Role**: `ecsInstanceRole`
		   
Click **Next**, then configure compute resources:

- **Min vCPUs**: `0`
- **Desired vCPUs**: `0`
- **Max vCPUs**: `256`  
  *(Example: allows spawning 8 `c7i.8xlarge` instances with 32 vCPUs each)*


- **Instance types**:
  - Remove default instance types
  - Add: `c7i.8xlarge`

- **Allocation strategy**: Best fit progressive

Click **Next**, then configure networking:

- **VPC**: Default
- **Subnets**: Default
- **Security group**: Default

Click **Next**, review the configuration, and click **Create compute environment**.

---

### 4.2 Create a Job Queue

1. Go to **AWS Console → Batch → Job queues → Create**
2. Set:
   - **Orchestration type**: EC2
   - **Name**: `<your-workload-name>-queue`
   - **Priority**: `999`
3. Attach the **Compute Environment** created above
4. Click **Create job queue**

> Nextflow will submit jobs to this queue.

---

## 5) Docker Setup

### 5.1 Build the Docker Container on Local Machine / EC2 Instance

From the directory containing your `Dockerfile`, run:

```bash
docker build -t <user>/<image>:latest .
```
#### 5.1.1 Install AWS CLI Inside the Container

Ensure AWS CLI is installed inside the Docker image by adding the following line to your `Dockerfile`:

```dockerfile
RUN pip install awscli
```
### 5.2 Push Image to Docker Hub

Log in to Docker Hub and push the Docker image:

```bash
docker login -u <user>
docker push <user>/<image>:latest
```

**Once pushed, your image will be available on Docker Hub.**

### 5.3 Use Image in nextflow.config

```bash
process {
    container = 'docker.io/<user_name>/<image_name>:<tag>'
}
```

--- 

## 6. S3 Bucket Setup

Create a bucket,  go to:

aws -> s3-> Create bucket->bucketname->create bucket e.g.:  in the same region we have configured our awscli :`srilekha-nextflow-workflows`

---

## 7. Nextflow Pipeline Setup on local instance/machine.
Install nextflow:

### With sudo

```bash
sudo apt update
sudo apt install -y openjdk-21-jre-headless
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
which nextflow
```

### Without sudo:

```bash
curl -L https://github.com/nextflow-io/nextflow/releases/download/v25.04.0/nextflow -o nextflow
chmod +x nextflow
export PATH=/home/hgx/omics/srilekhx:$PATH
which nextflow
```

To make it permanent: so you don’t need export every time)
```bash
echo 'export PATH="$HOME/omics/srilekhx:$PATH"' >> ~/.bashrc
source ~/.bashrc
which nextflow
```

### Run command:
```bash
nextflow run main.nf   -profile awsbatch_single_node
```


### 7.1 main.nf (Define the channel(for tasks(inputs), process(command), workflow to call process with input)

```bash
nextflow.enable.dsl = 2
params.pdb_list = params.pdb_list ?: 'protein_list.txt'
params.pdb_dir  = params.pdb_dir  ?: 'inputs'
workflow {
    Channel
      .fromPath(params.pdb_list)
      .splitText()
      .set { jobs }
    RUN_GROMACS(jobs)
}
process RUN_GROMACS {
    input:
    val pdb_name
    output:
    path "out/${pdb_name}"
    script:
    """
    bash run_commands.sh ${pdb_name}
    """
}
```

### 7.2 nextflow.config

```bash
profiles {
  awsbatch {
    process.executor = 'awsbatch'
    workDir = 's3://srilekha-nextflow-workflows/work'
    process {
      queue     = 'gromacs-singlenode-queue'
      container = 'docker.io/809314/grms_opom:latest'
      cpus      = 32
      memory    = '60 GB'
    }
    executor {
      queueSize = 8   // max parallel Batch jobs
    }
  }
}
```

### Download data from s3 bucket to local machine:

```bash
mkdir -p backup
aws s3 cp s3://srilekha-new-s3/outputs-motif-woqueue ./backup –recursive
```

### To remove data from bucket:
```bash
aws s3 rm s3://srilekha-new-s3/ outputs-motif-woqueue / --recursive
```

### Command to run:

```bash
nextflow run main.nf   -profile awsbatch_single_node   --tsv contigs.tsv   --input_dir pdbs   --out_s3 "$OUT_S3"
```

### nextflow.config:
```bash
profiles {

  awsbatch_single_node {

    aws.region = 'us-east-2'
    workDir    = 's3://srilekha-new-s3/work-motif-demo'

    executor {
      // How many jobs Nextflow can submit to Batch at once
      queueSize = 16
    }

    process {
      executor  = 'awsbatch'

      // Your existing Batch queue name
      queue     = 'srilekha_8_11_queue_motif'

      // Docker image with RFdiffusion installed
      container = 'docker.io/809314/rfdiff-aws-rec:latest'

      cpus      = 32
      memory    = '60 GB'

      withName: runMotifScaffolding {
        // How many RFdiffusion jobs can run in parallel
        maxForks = 8
      }
    }
  }
}
```

### Main.nf:
```bash
nextflow.enable.dsl = 2

// ----------------------
// Parameters
// ----------------------
params.tsv       = params.tsv       ?: 'contigs_8.tsv'
params.input_dir = params.input_dir ?: 'pdbs'
params.out_s3    = params.out_s3    ?: 's3://srilekha-new-s3/outputs-motif-demo'

// ----------------------
// Create jobs channel from TSV
// ----------------------
Channel
  .fromPath(params.tsv)
  .ifEmpty { error "TSV not found: ${params.tsv}" }
  .splitText()
  .map { it.trim() }
  .filter { it && !it.startsWith('#') }     // skip empty & comment lines
  .map { line ->
    def parts = line.split(/\s+/)
    if (parts.size() < 2)
      throw new IllegalArgumentException("Bad TSV line: '${line}' (need: PDB CONTIG)")
    def pdb    = parts[0]
    def contig = parts[1]
    def base   = pdb.replaceAll(/\.pdb$/, '')
    def prefix = "design_${base}"
    tuple(pdb, contig, prefix)
  }
  .map { pdbName, contig, prefix ->
    // Build full path to PDB
    def full = params.input_dir.endsWith('/') ?
               "${params.input_dir}${pdbName}" :
               "${params.input_dir}/${pdbName}"
    tuple(file(full), contig, prefix)
  }
  .set { jobs }

// ----------------------
// RFdiffusion process
// ----------------------
process runMotifScaffolding {

  tag { prefix }

  publishDir params.out_s3, mode: 'copy', overwrite: false

  input:
  tuple path(pdb_file), val(contig), val(prefix)

  output:
  path "out/${prefix}",        emit: designed_dir
  path "logs/${prefix}.log",   emit: joblog

  shell:
  '''
  set -euo pipefail

  mkdir -p logs
  LOGFILE="logs/!{prefix}.log"

  {
    echo "=== RFdiffusion motif run ==="
    echo "PDB:    !{pdb_file}"
    echo "CONTIG: !{contig}"
    echo "PREFIX: !{prefix}"
    echo

    # Working directory for RFdiffusion
    mkdir -p run
    cp -f "!{pdb_file}" "run/"
    pdb_local="run/$(basename "!{pdb_file}")"

    # Run RFdiffusion
    python /app/RFdiffusion/scripts/run_inference.py \
        inference.output_prefix=./run/test/!{prefix} \
        inference.input_pdb="$pdb_local" \
        "contigmap.contigs=[\"!{contig}\"]" \
        inference.num_designs=1 \
        inference.precision=bfloat16 \
        inference.deterministic=True

    # Collect outputs in a clean folder
    mkdir -p "out/!{prefix}"
    if [ -d run/test ]; then
        cp -r run/test/* "out/!{prefix}/" || true
    fi

  } &> "$LOGFILE"
  '''
}

// ----------------------
// Workflow entrypoint
// ----------------------
workflow {
  runMotifScaffolding(jobs)
}
```
---

## Workflow Structure:
```text
PDB + scripts → S3 workDir
        ↓
AWS Batch job (Docker)
        ↓
Task outputs in workDir
        ↓
Nextflow publishes to S3 out_s3
```

---

## Additional IAM Policies (RELION):

### BatchComputeEnvAccess:
```bash
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Sid": "BatchComputeEnvironmentCRUD",
            "Effect": "Allow",
            "Action": [
                "batch:CreateComputeEnvironment",
                "batch:DescribeComputeEnvironments",
                "batch:UpdateComputeEnvironment",
                "batch:DeleteComputeEnvironment",
                "batch:DescribeJobQueues",
                "batch:UpdateJobQueue"
            ],
            "Resource": "*"
        }
    ]
}
```
### LaunchTemplateAccess:
```bash
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Sid": "LaunchTemplateCRUD",
            "Effect": "Allow",
            "Action": [
                "ec2:CreateLaunchTemplate",
                "ec2:CreateLaunchTemplateVersion",
                "ec2:DescribeLaunchTemplates",
                "ec2:DescribeLaunchTemplateVersions",
                "ec2:DeleteLaunchTemplate",
                "ec2:DeleteLaunchTemplateVersions"
            ],
            "Resource": "*"
        }
    ]
}
```

### EC2DescribeInstancesVolumesReadOnly:
```bash
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Sid": "VisualEditor0",
            "Effect": "Allow",
            "Action": [
                "ec2:DescribeInstances",
                "ec2:DescribeVolumes"
            ],
            "Resource": "*"
        }
    ]
}
```

### IAMReadOnlyAccess:
```bash
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "iam:GenerateCredentialReport",
                "iam:GenerateServiceLastAccessedDetails",
                "iam:Get*",
                "iam:List*",
                "iam:SimulateCustomPolicy",
                "iam:SimulatePrincipalPolicy"
            ],
            "Resource": "*"
        }
    ]
}
```

