# Using ProteinMPNN as a Microservice

This directory provides a microservice interface for running ProteinMPNN. It allows you to expose ProteinMPNN as a server and interact with it via HTTP-based requests, enabling easy integration into larger workflows or client applications.

## Overview

The setup includes two main components:

1. Server – Hosts the ProteinMPNN model and exposes an API for inference.

2. Client – Sends requests to the server with required inputs and receives predictions as responses.

## Prerequisites

Before proceeding, ensure you have already created the ProteinMPNN Docker image using the instructions [here](../README.md#-using-docker)

## Running the Microservice

### 1. Start the Server

Start the ProteinMPNN microservice container:

```bash
docker run -it pmpnn:latest python microservice/opea_pmpnn_server.py
```

* This will launch the server on port 9100 and expose the inference API at `http://<host>:9100/v1/protein_mpnn`.

* To use a different port, specify it with the --port argument.

#### Server Ready Check

The server is ready to accept requests when you see logs similar to the following:
Once the server has started, you can find the API endpoint highlighted in bold within the logs.
```bash
[2025-11-04 11:27:18,870] [    INFO] - opea_service_protein_mpnn_microservice - Using component name: OPEA_OMICS_PROTEINMPNN
✅ Model loaded in 0.08 seconds
[2025-11-04 11:27:18,953] [    INFO] - opea_service_protein_mpnn_microservice - Starting ProteinMPNN microservice at http://172.31.47.52:9100/v1/protein_mpnn
[2025-11-04 11:27:18,957] [    INFO] - Base service - CORS is enabled.
[2025-11-04 11:27:18,958] [    INFO] - Base service - Setting up HTTP server
[2025-11-04 11:27:18,958] [    INFO] - Base service - Uvicorn server setup on port 9100
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9100 (Press CTRL+C to quit)
[2025-11-04 11:27:18,987] [    INFO] - Base service - HTTP server setup successful
[2025-11-04 11:27:18,988] [    INFO] - opea_service_protein_mpnn_microservice - ProteinMPNN microservice server started.
```


### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the ProteinMPNN microservice.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

#### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.31.47.52:9100/v1/protein_mpnn
```

On the client machine, export the following environment variables:

```bash
export HOST="172.31.47.52" # Replace with the actual server IP
export PORT="9100" # Replace with the actual server port
```

You only need to do this once per session.

### Step 2: Run Inference Examples

Below is example of how to invoke the client for ProteinMPNN runs.

Note: Various input parameters to protein_mpnn_run.py are described in the original ProteinMPNN README. Please refer to: [here](../README.md#-using-docker)
### 🧭 Command-line Help

This script supports two modes:

🔹 standalone – run ProteinMPNN locally.

🔹 microservice – run ProteinMPNN through an API server.

### 🔹 Show help for microservice mode
```bash
python ../examples/script_example_1.py microservice --help
```

#### Example: 
#### Simple monomer example
```bash
python ../examples/script_example_1.py  microservice --input <pdbs_folder>  --output <outputs/example_1_outputs> --num_seq_per_target <int> --sampling_temp <float> --seed <int> --batch_size <int>  --host $HOST  --port $PORT
```
### Simple multi-chain example
```bash
python ../examples/script_example_2.py microservice --input <pdbs_folder>  --output <outputs/example_2_outputs>  --host $HOST  --port $PORT
```
### Directly from the .pdb path
```bash
python ../examples/script_example_3.py microservice --pdb_path <5TPN.pdb>  --output <outputs/example_3_outputs>   --host $HOST  --port $PORT
```
### Return score only (model's uncertainty)
```bash
python ../examples/script_example_3_score_only.py microservice --pdb_path <5TPN.pdb>  --output <outputs/example_3_score_only_outputs>   --host $HOST  --port $PORT
```
### Return score only (model's uncertainty) loading sequence from fasta files
```bash
python ../examples/script_example_3_score_only_from_fasta.py microservice --pdb_path <5TPN.pdb> --path_to_fasta <outputs/example_3_outputs/seqs/5TPN.fa>  --output <outputs/example_3_score_only_from_fasta_output>  --host $HOST  --port $PORT
```
### Fix some residue positions
```bash
python ../examples/script_example_4.py microservice --input <pdb_folder>  --output <outputs/example_4_outputs>  --host $HOST  --port $PORT
```
### Specify which positions to design
```bash
python ../examples/script_example_4_non_fixed.py microservice --input <pdb_folder>  --output <outputs/example_4_non_fixed_outputs>  --host $HOST  --port $PORT
```
### Tie some positions together (symmetry)
```bash
python ../examples/script_example_5.py microservice --input <pdb_folder>  --output <outputs/example_5_outputs>  --host $HOST  --port $PORT
```
### Homooligomer example
```bash
python ../examples/script_example_6.py microservice --input <pdb_folder>  --output <outputs/example_6_outputs>  --host $HOST  --port $PORT
```
### Return sequence unconditional probabilities (PSSM like)
```bash
python ../examples/script_example_7.py microservice --input <pdb_folder>  --output <outputs/example_7_outputs>  --host $HOST  --port $PORT
```
### Add amino acid bias
```bash
python ../examples/script_example_8.py microservice --input <pdb_folder>  --output <outputs/example_8_outputs>  --host $HOST  --port $PORT
```
### Use PSSM bias when designing sequences
```bash
python ../examples/script_example_pssm.py microservice --pssm_input <PSSM_inputs> --input <pdb_folder>  --output <outputs/example_pssm_outputs>  --host $HOST  --port $PORT
```




