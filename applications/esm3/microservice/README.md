# Using ESM3 as a Microservice

This directory provides a microservice interface for running ESM3. It allows you to expose ESM3 as a server and interact with it via HTTP-based requests, enabling easy integration into larger workflows or client applications.

## Overview

The setup includes two main components:

1. Server – Hosts the ESM3 model and exposes an API for inference.

2. Client – Sends requests to the server with required inputs and receives predictions as responses.

## Prerequisites

Before proceeding, ensure that you have already built the ESM3 Docker image following the provided instructions and exported the required directories [here](../README.md#Building-the-Docker-Image).

## Running the Microservice Experiments

### Workflow 1: ESMC

#### 1. Start the Server

Start the ESMC microservice container:

```bash
docker run -it
    -e HF_TOKEN="<your_huggingface_token>" \
    -v $MODELS:/models \
    -v $INPUT:/input \
    -v $OUTPUT:/output \
    esm3_image:latest \
    python microservice/esmc_tasks/esmc_embeding_microservice.py
```

* This will launch the server on port 9009 and expose the inference API at `http://<host>:9009/v1/esmcembedding`.

* To use a different port, specify it with the `--port` argument.
* `--bf16` flag accelerates performance by utilizing bfloat16 precision, which enhances computational efficiency without compromising accuracy.


#### Server Ready Check

The server is ready to accept requests when you see logs similar to the following:

```bash
$ docker run -it     -e HF_TOKEN="<your_huggingface_token>"     -v $MODELS:/models     -v $INPUT:/input     -v $OUTPUT:/output     esm3_image:latest     python microservice/esmc_tasks/esmc_embeding_microservice.py
WARNING: The MAMBA_ROOT_PREFIX environment variable is not set.
WARNING: This is required for mamba to work correctly as of 2.0.
WARNING:
WARNING: For now, we are setting 'MAMBA_ROOT_PREFIX' to '/home/esm3-base-service/conda'.
WARNING:
WARNING: Please make sure this is consistent with your installation or alternatively (by order of preference):
WARNING:   - rerun 'mamba shell init' to initialize mamba for your current shell
WARNING:   - manually set 'MAMBA_ROOT_PREFIX' to the root of your installation in your shell profile script.
WARNING:   - use the '-r,--root-prefix' CLI option when calling mamba.
WARNING:
WARNING: This message originates from /home/esm3-base-service/conda/etc/profile.d/mamba.sh
Using model: esmc_600m
Selected data type ESMC : torch.float32
model_name ESMC esmc_600m
Fetching 4 files: 100%|█████████████████████████████████████████████████████████████████████████████| 4/4 [00:00<00:00, 15621.24it/s]
Model load time 18.999433517456055
[2025-11-06 06:13:38,575] [    INFO] - opea_service_omics_esmcembedding_microservice - Starting esmcembedding microservice at http://172.17.0.2:9009/v1/esmcembedding
[2025-11-06 06:13:38,576] [    INFO] - Base service - CORS is enabled.
[2025-11-06 06:13:38,577] [    INFO] - Base service - Setting up HTTP server
[2025-11-06 06:13:38,578] [    INFO] - Base service - Uvicorn server setup on port 9009
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9009 (Press CTRL+C to quit)
```


#### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the esmc microservice to run various tasks like generate and optimize.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

##### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.2:9009/v1/esmcembedding
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.2" # Replace with the actual server IP
export PORT="9009" # Replace with the actual server port
```

You only need to do this once per session.

##### Step 2: Run Inference Examples

Below are examples of how to invoke the client for different esmc use-cases.

```bash
#example
export INPUT_FILE="some_proteins.fasta"
python esmc_tasks/esmc_embeding_client.py.py \
    --host $HOST \
    --port $PORT \
    --fasta_file  $INPUT_FILE
```

### Workflow 2: ESM3

#### 1. Start the Server

Start the ESM3 microservice container:

```bash
docker run -it     
  -e HF_TOKEN="<your_huggingface_token>"     
  -v $MODELS:/models     
  -v $INPUT:/input     
  -v $OUTPUT:/output     
  esm3_image:latest     
  python microservice/esm3_tasks/esm3_task_microservice.py --port 9008
```

* This will launch the server on port 9008 and expose the inference API at `http://<host>:9008/v1/esm3task`.

* To use a different port, specify it with the `--port` argument.
* `--bf16` flag accelerates performance by utilizing bfloat16 precision, which enhances computational efficiency without compromising accuracy.

#### Server Ready Check

The server is ready to accept requests when you see logs similar to the following:

```bash
$ docker run -it     -e HF_TOKEN="<your_huggingface_token>"     -v $MODELS:/models     -v $INPUT:/input     -v $OUTPUT:/output     esm3_image:latest     python microservice/esm3_tasks/esm3_task_microservice.py --port 9008
WARNING: The MAMBA_ROOT_PREFIX environment variable is not set.
WARNING: This is required for mamba to work correctly as of 2.0.
WARNING:
WARNING: For now, we are setting 'MAMBA_ROOT_PREFIX' to '/home/esm3-base-service/conda'.
WARNING:
WARNING: Please make sure this is consistent with your installation or alternatively (by order of preference):
WARNING:   - rerun 'mamba shell init' to initialize mamba for your current shell
WARNING:   - manually set 'MAMBA_ROOT_PREFIX' to the root of your installation in your shell profile script.
WARNING:   - use the '-r,--root-prefix' CLI option when calling mamba.
WARNING:
WARNING: This message originates from /home/esm3-base-service/conda/etc/profile.d/mamba.sh
No ESM_API_KEY found. Loading model locally...
Selected data type: torch.float32
Fetching 22 files: 100%|█████████████████████████████████████████████████████████████████████████| 22/22 [00:00<00:00, 258472.52it/s]
Model load time 7.015205144882202
[2025-11-05 11:42:18,899] [    INFO] - opea_service_omics_esm3embedding_microservice - Starting esm3embedding microservice at http://172.17.0.2:9008/v1/esm3task
[2025-11-05 11:42:18,900] [    INFO] - Base service - CORS is enabled.
[2025-11-05 11:42:18,901] [    INFO] - Base service - Setting up HTTP server
[2025-11-05 11:42:18,902] [    INFO] - Base service - Uvicorn server setup on port 9008
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9008 (Press CTRL+C to quit)
[2025-11-05 11:42:18,911] [    INFO] - Base service - HTTP server setup successful
[2025-11-05 11:42:18,912] [    INFO] - opea_service_omics_esm3embedding_microservice - OPEA_OMICS_esm3embedding server started.
```

#### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the esm3 microservice to run various tasks like generate and optimize.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

##### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.7:9008/v1/esm3task
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.2" # Replace with the actual server IP
export PORT="9008" # Replace with the actual server port
```

You only need to do this once per session.

##### Step 2: Run Inference Examples

Below are examples of how to invoke the client for different ESM3 use-cases.


##### ESM3 - Logits Embeddings: Generates embeddings with ESM3 for sequence representation and analysis

```bash
#example
export INPUT_FILE="some_proteins.fasta"
python esmc_tasks/esm3_task_client.py.py \
    --host $HOST \
    --port $PORT \
    --fasta_file  $INPUT_FILE
    --task logits_embedding
```

##### ESM3 - Folding: Predicts the 3D structure of proteins from amino acid sequences using ESM3

```bash
#example
export INPUT_FILE="some_proteins.fasta"
python esmc_tasks/esm3_task_client.py.py \
    --host $HOST \
    --port $PORT \
    --fasta_file  $INPUT_FILE
    --task fold
```

##### ESM3 - Inverse Folding: Designs protein sequences that fold into a given 3D structure
```bash
#example
export INPUT_FILE="5YH2.pdb"
python esmc_tasks/esm3_task_client.py.py \
    --host $HOST \
    --port $PORT \
    --pdb_file  $INPUT_FILE
    --task inversefold
```

##### ESM3 - Function Prediction: Predicts protein function from structural and sequence data
```bash
#example
export INPUT_FILE="5YH2.pdb"
python esmc_tasks/esm3_task_client.py.py \
    --host $HOST \
    --port $PORT \
    --pdb_file  $INPUT_FILE
    --task function_prediction
```

##### ESM3 - Prompt Sequence: Generates protein sequences based on user-provided prompts for design tasks
```bash
#example
export INPUT_FILE="prompt.fasta"
python esmc_tasks/esm3_task_client.py.py \
    --host $HOST \
    --port $PORT \
    --fasta_file  $INPUT_FILE
    --task prompt_task
```

##### ESM3 - Chain of Thought: Uses reasoning-based approaches to analyze and interpret protein data
```bash
#example
export INPUT_FILE="1utn.csv"
python esmc_tasks/esm3_task_client.py.py \
    --host $HOST \
    --port $PORT \
    --csv_file  $INPUT_FILE
    --task chain_of_thought
```
