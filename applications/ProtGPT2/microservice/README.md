# Using ProtGPT2 as a Microservice

This directory provides a microservice interface for running ProtGPT2. It allows you to expose ProtGPT2 as a server and interact with it via HTTP-based requests, enabling easy integration into larger workflows or client applications.

## Overview

The setup includes two main components:

1. Server – Hosts the ProtGPT2 model and exposes an API for inference.

2. Client – Sends requests to the server with required inputs and receives predictions as responses.

## Prerequisites

Before proceeding, ensure you have already created the ProtGPT2 Docker image using the instructions [here](../README.md#-using-docker)

## Running the Microservice

### 1. Start the Server

Start the ProtGPT2 microservice container:

```bash
docker run -it protgpt2:latest python microservice/protgpt2_opea_server.py
```

* This will launch the server on port 9010 and expose the inference API at `http://<host>:9010/v1/protgpt2`.

* To use a different port, specify it with the --port argument.

#### Server Ready Check

The server is ready to accept requests when you see logs similar to the following:
Once the server has started, you can find the API endpoint highlighted in bold within the logs.
```bash
/opt/conda/envs/protgpt2/lib/python3.11/site-packages/intel_extension_for_pytorch/nn/utils/_weight_prepack.py:5: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
  import pkg_resources
[2025-09-03 08:41:20,454] [    INFO] - ProtGPT2 - Loading model
[2025-09-03 08:41:20,454] [    INFO] - ProtGPT2 - Loading ProtGPT2 model: nferruz/ProtGPT2 with dtype=float32
/opt/conda/envs/protgpt2/lib/python3.11/site-packages/huggingface_hub/file_download.py:945: FutureWarning: `resume_download` is deprecated and will be removed in version 1.0.0. Downloads always resume when possible. If you want to force a new download, use `force_download=True`.
  warnings.warn(
config.json: 100%|██████████████████████████████████████████████████████████████████████████████████| 850/850 [00:00<00:00, 2.96MB/s]
pytorch_model.bin: 100%|████████████████████████████████████████████████████████████████████████| 3.13G/3.13G [04:29<00:00, 11.6MB/s]
vocab.json: 655kB [00:00, 25.0MB/s]
merges.txt: 314kB [00:00, 20.3MB/s]
tokenizer.json: 1.07MB [00:00, 22.4MB/s]
special_tokens_map.json: 100%|██████████████████████████████████████████████████████████████████████| 357/357 [00:00<00:00, 1.36MB/s]
[2025-09-03 08:45:54,296] [    INFO] - ProtGPT2 - Applying IPEX optimization
/opt/conda/envs/protgpt2/lib/python3.11/site-packages/intel_extension_for_pytorch/frontend.py:462: UserWarning: Conv BatchNorm folding failed during the optimize process.
  warnings.warn(
/opt/conda/envs/protgpt2/lib/python3.11/site-packages/intel_extension_for_pytorch/frontend.py:469: UserWarning: Linear BatchNorm folding failed during the optimize process.
  warnings.warn(
[2025-09-03 08:45:54,637] [    INFO] - ProtGPT2 - Model loaded in 274.18 seconds
[2025-09-03 08:45:54,638] [    INFO] - opea_service_omics_protgpt2_microservice - Starting ProtGPT2 microservice at **http://172.17.0.4:9010/v1/protgpt2**
[2025-09-03 08:45:54,639] [    INFO] - Base service - CORS is enabled.
[2025-09-03 08:45:54,640] [    INFO] - Base service - Setting up HTTP server
[2025-09-03 08:45:54,642] [    INFO] - Base service - Uvicorn server setup on port 9010
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9010 (Press CTRL+C to quit)
[2025-09-03 08:45:54,657] [    INFO] - Base service - HTTP server setup successful
[2025-09-03 08:45:54,658] [    INFO] - opea_service_omics_protgpt2_microservice - OPEA_OMICS_PROTGPT2 server started.
```


### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the ProtGPT2 microservice.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

#### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.4:9010/v1/protgpt2
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.4" # Replace with the actual server IP
export PORT="9010" # Replace with the actual server port
```

You only need to do this once per session.

### Step 2: Run Inference Examples

Below is example of how to invoke the client for ProtGPT2 runs.

#### Example: 

```bash
python protgpt2_opea_client.py --host $HOST --port $PORT --iterations 5 --seed 42 --max_length 100 --num_return_sequences 10
```
