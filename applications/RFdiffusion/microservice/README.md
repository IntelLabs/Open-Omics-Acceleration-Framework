# Using RFdiffusion as a Microservice

This directory provides a microservice interface for running RFdiffusion. It allows you to expose RFdiffusion as a server and interact with it via HTTP-based requests, enabling easy integration into larger workflows or client applications.

## Overview

The setup includes two main components:

1. Server – Hosts the RFdiffusion model and exposes an API for inference.

2. Client – Sends requests to the server with required inputs and receives predictions as responses.

## Prerequisites

Before proceeding, ensure you have already created the RFdiffusion Docker image using the instructions [here](../README.md#-using-docker)

## Running the Microservice

### 1. Start the Server

Start the RFdiffusion microservice container:

```bash
docker run --rm rfdiffusion:latest python ../microservice/opea_rfdiffusion_microservice.py --bfloat16
```

* This will launch the server on port 9000 and expose the inference API at `http://<host>:9000/v1/rfdiffusion`.

* To use a different port, specify it with the --port argument.

#### Server Ready Check

The server is ready to accept requests when you see logs similar to the following:

```bash
$ docker run rfdiffusion:latest python ../microservice/opea_rfdiffusion_microservice.py --bfloat16
/app/RFdiffusion/rfdiffusion/util.py:253: UserWarning: Using torch.cross without specifying the dim arg is deprecated.
Please either pass the dim explicitly or simply use torch.linalg.cross.
The default value of dim will change to agree with that of linalg.cross in a future release. (Triggered internally at ../aten/src/ATen/native/Cross.cpp:63.)
  Z = torch.cross(Xn, Yn)
DGL backend not selected or invalid.  Assuming PyTorch for now.
/opt/conda/envs/SE3nv/lib/python3.11/site-packages/intel_extension_for_pytorch/nn/utils/_weight_prepack.py:5: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
  import pkg_resources
[2025-09-16 13:07:30,738] [    INFO] - rfdiffusion_microservice - Info: Running in bfloat16 mode.
[2025-09-16 13:07:30,738] [    INFO] - opea_service_omics_rfdiffusion_microservice - Starting RFdiffusion microservice at http://172.17.0.6:9000/v1/rfdiffusion
[2025-09-16 13:07:30,739] [    INFO] - Base service - CORS is enabled.
[2025-09-16 13:07:30,740] [    INFO] - Base service - Setting up HTTP server
[2025-09-16 13:07:30,741] [    INFO] - Base service - Uvicorn server setup on port 9000
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9000 (Press CTRL+C to quit)
[2025-09-16 13:07:30,752] [    INFO] - Base service - HTTP server setup successful
[2025-09-16 13:07:30,753] [    INFO] - rfdiffusion_microservice - Work dir created: /tmp/rfdiffusion_inputs/server_9000
[2025-09-16 13:07:30,753] [    INFO] - opea_service_omics_rfdiffusion_microservice - OPEA_OMICS_RFDIFFUSION server started.
```


### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the RFdiffusion microservice to run various tasks like motif scaffolding and partial diffusion.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

Note: Before running the client program, make sure the default config file is set up at [here](./config/).

If you modify the config file, you must provide the correct path when running opea_rfdiffusion_client.py.
#### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.2:9000/v1/rfdiffusion
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.2" # Replace with the actual server IP
export PORT="9000" # Replace with the actual server port
```

You only need to do this once per session.

### Step 2: Run Inference Examples

Below are examples of how to invoke the client for different RFdiffusion use-cases.

#### Example: Motif Scaffolding

```bash

cd applications/RFdiffusion/microservice

export INPUT_FILE=../examples/input_pdbs/5TPN.pdb

python opea_rfdiffusion_client.py \
  --host $HOST \
  --port $PORT \
  inference.input_pdb=$INPUT_FILE \
  'contigmap.contigs=[10-40/A163-181/10-40]' \
  inference.num_designs=1

```

#### Example: Partial Diffusion

```bash

cd applications/RFdiffusion/microservice

export INPUT_FILE=../examples/input_pdbs/2KL8.pdb

python opea_rfdiffusion_client.py \
  --host $HOST \
  --port $PORT \
  inference.input_pdb=$INPUT_FILE \
  'contigmap.contigs=[79-79]' \
  inference.num_designs=1 \
  diffuser.partial_T=10

```

