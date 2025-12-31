# Using ESM2 as a Microservice

This directory provides a microservice interface for running ESM2. It allows you to expose ESM2 as a server and interact with it via HTTP-based requests, enabling easy integration into larger workflows or client applications.

## Overview

The setup includes two main components:

1. Server – Hosts the ESM model and exposes an API for inference.

2. Client – Sends requests to the server with required inputs and receives predictions as responses.

## Prerequisites
Before proceeding, ensure that you have already built the ESM2 Docker image using the provided instructions and exported the necessary directories as described [here](../README.md#Installation).

## Running the Microservice ESM2

### Embedding

#### 1. Start the Server

Start the esm embedding microservice container:


```bash
docker run -it \
  -v $MODELS:/checkpoints \
  -v $INPUT:/input \
  -v $OUTPUT:/output \
  esm_image:latest \
  python omics_setup/microservice/embedding/opea_embedding_microservice.py --model_name esm2_t33_650M_UR50D --port 9001
```

* This will launch the server on port 9001 and expose the inference API at `http://<host>:9001/v1/esmembedding`.

* To use a different port, specify it with the `--port` argument.
* `--bf16` flag accelerates performance by utilizing bfloat16 precision, which enhances computational efficiency without compromising accuracy.

Server Ready Check:

The server is ready to accept requests when you see logs similar to the following:

```bash
$ docker run esm:latest python ../micro
service/opea_esm_microservice.py
/app/esm/esm/util.py:253: UserWarning: Using torch.cross without specifying the dim arg is deprecated.
Please either pass the dim explicitly or simply use torch.linalg.cross.
The default value of dim will change to agree with that of linalg.cross in a future release. (Triggered internally at ../aten/src/ATen/native/Cross.cpp:63.)
  Z = torch.cross(Xn, Yn)
DGL backend not selected or invalid.  Assuming PyTorch for now.
/opt/conda/envs/SE3nv/lib/python3.11/site-packages/intel_extension_for_pytorch/nn/utils/_weight_prepack.py:5: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
  import pkg_resources
[2025-07-28 06:20:19,240] [    INFO] - esm_microservice - Starting esm microservice at http://172.17.0.2:9001/v1/esmembedding
[2025-07-28 06:20:19,245] [    INFO] - Base service - CORS is enabled.
[2025-07-28 06:20:19,245] [    INFO] - Base service - Setting up HTTP server
[2025-07-28 06:20:19,246] [    INFO] - Base service - Uvicorn server setup on port 9001
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9001 (Press CTRL+C to quit)
[2025-07-28 06:20:19,278] [    INFO] - Base service - HTTP server setup successful
```


#### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the esm microservice to run various tasks like embedding and esmfold.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

##### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.2:9001/v1/esmembedding
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.2" # Replace with the actual server IP
export PORT="9001" # Replace with the actual server port
```

You only need to do this once per session.

##### Step 2: Run Inference Examples

Below are example of how to invoke the client for esm embedig.

```bash
#example
export INPUT=../examples/data/few_proteins.fasta

python embedding/opea_embedding_client.py \
  --host $HOST \
  --port $PORT \
  --fasta_file=$INPUT \
  --repr_layers -1 \
  --include mean
```

### ESMFold

#### 1. Start the Server

Start the esmfold microservice container:

```bash
docker run -it \
  -v $MODELS:/checkpoints \
  -v $INPUT:/input \
  -v $OUTPUT:/output \
  esmfold_image:latest \
  python omics_setup/microservice/esmfold/opea_esmfold_microservice.py --port 9002
```

* This will launch the server on port 9002 and expose the inference API at `http://<host>:9002/v1/esmfold`.

* To use a different port, specify it with the `--port` argument.
* `--bf16` flag accelerates performance by utilizing bfloat16 precision, which enhances computational efficiency without compromising accuracy.

Server Ready Check

The server is ready to accept requests when you see logs similar to the following:

```bash
$ docker run -it   -v $MODELS:/checkpoints   -v $INPUT:/input   -v $OUTPUT:/output   esmfold_image:latest   python omics_setup/microservice/esmfold/opea_esmfold_microservice.py --port 9002
[2025-10-08 11:57:22,785] [    INFO] - esm_fold_microservice - Loading model
Model load time 49.38147163391113
[2025-10-08 11:58:12,167] [    INFO] - opea_service_omics_esmfold_microservice - Starting esmfold microservice at http://172.17.0.2:9002/v1/esmfold
[2025-10-08 11:58:12,171] [    INFO] - Base service - CORS is enabled.
[2025-10-08 11:58:12,173] [    INFO] - Base service - Setting up HTTP server
[2025-10-08 11:58:12,174] [    INFO] - Base service - Uvicorn server setup on port 9002
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9002 (Press CTRL+C to quit)
[2025-10-08 11:58:12,182] [    INFO] - Base service - HTTP server setup successful
[2025-10-08 11:58:12,184] [    INFO] - opea_service_omics_esmfold_microservice - OPEA_OMICS_esmfold server started.
```


#### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the esmfold microservice.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

##### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.2:9002/v1/esmfold
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.2" # Replace with the actual server IP
export PORT="9002" # Replace with the actual server port
```

You only need to do this once per session.

##### Step 2: Run Inference Examples

Below are example of how to invoke the client for esmfold.

```bash
#example
export INPUT=../examples/data/few_proteins.fasta

python esmfold/opea_esmfold_client.py \
  --host $HOST \
  --port $PORT \
  --fasta_file=$INPUT
```


### InverseFold

#### Sample sequence designs for a given structure

###### 1. Start the Server

Start the esm inversefold sample sequence microservice container:

```bash
 docker run -it \
  -v $MODELS:/checkpoints  \
  -v $INPUT:/input  \
  -v $OUTPUT:/output  \
  esm_image:latest   \
  python omics_setup/microservice/inversefold/sample_sequence/opea_sample_sequence_microservice.py
```

* This will launch the server on port 9003 and expose the inference API at `http://<host>:9003/v1/esminversefold`.

* To use a different port, specify it with the `--port` argument.
* `--bf16` flag accelerates performance by utilizing bfloat16 precision, which enhances computational efficiency without compromising accuracy.

Server Ready Check

The server is ready to accept requests when you see logs similar to the following:

```bash
$ docker run -it   -v $MODELS:/checkpoints   -v $INPUT:/input   -v $OUTPUT:/output   esm_image:latest   python omics_setup/microservice/inversefold/sample_sequence/opea_sample_sequence_microservice.py --port 9003
/home/esm-base-service/conda/envs/esm_py11/lib/python3.11/site-packages/esm/pretrained.py:215: UserWarning: Regression weights not found, predicting contacts will not produce correct results.
  warnings.warn(
Model load time 2.733304738998413
[2025-10-08 05:41:15,507] [    INFO] - opea_service_omics_esminversefold_microservice - Starting esminversefolding sample sequence microservice at http://172.17.0.2:9003/v1/esminversefold
[2025-10-08 05:41:15,509] [    INFO] - Base service - CORS is enabled.
[2025-10-08 05:41:15,510] [    INFO] - Base service - Setting up HTTP server
[2025-10-08 05:41:15,511] [    INFO] - Base service - Uvicorn server setup on port 9003
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9003 (Press CTRL+C to quit)
[2025-10-08 05:41:15,523] [    INFO] - Base service - HTTP server setup successful
```


##### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the esm inversefold sample sequence microservice.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

##### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.2:9003/v1/esminversefold
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.2" # Replace with the actual server IP
export PORT="9003" # Replace with the actual server port
```

You only need to do this once per session.

##### Step 2: Run Inference Examples

Below are example of how to invoke the client for esm inversefold sample sequence.

```bash
#example
export INPUT_PDB=../../../examples/inverse_folding/data/5YH2.pdb 

python opea_sample_sequence_client.py \
  --host $HOST \
  --port $PORT \
  --pdb_file $INPUT_PDB \
  --chain C \
  --temperature 1 \
  --num_samples 3 
```

#### Scoring sequences


###### 1. Start the Server

Start the esm inversefolding score microservice container:

```bash
docker run -it  \
  -v $MODELS:/checkpoints  \
  -v $INPUT:/input  \
  -v $OUTPUT:/output  \
  esm_image:latest  \
  python omics_setup/microservice/inversefold/scoring_sequences/opea_scoring_sequence_microservice.py
```

* This will launch the server on port 9004 and expose the inference API at `http://<host>:9004/v1/esminversefoldscore`.

* To use a different port, specify it with the `--port` argument.
* `--bf16` flag accelerates performance by utilizing bfloat16 precision, which enhances computational efficiency without compromising accuracy.

Server Ready Check

The server is ready to accept requests when you see logs similar to the following:

```bash
$ docker run -it   -v $MODELS:/checkpoints   -v $INPUT:/input   -v $OUTPUT:/output   esm_image:latest   python omics_setup/microservice/inversefold/scoring_sequences/opea_scoring_sequence_microservice.py --port 9004
/home/esm-base-service/conda/envs/esm_py11/lib/python3.11/site-packages/esm/pretrained.py:215: UserWarning: Regression weights not found, predicting contacts will not produce correct results.
  warnings.warn(
Model load time 5.408844232559204
[2025-10-08 12:17:43,020] [    INFO] - opea_service_omics_esminversefold_score_microservice - Starting esminversefolding score sequence microservice microservice microservice at http://172.17.0.2:9004/v1/esminversefoldscore
[2025-10-08 12:17:43,022] [    INFO] - Base service - CORS is enabled.
[2025-10-08 12:17:43,023] [    INFO] - Base service - Setting up HTTP server
[2025-10-08 12:17:43,024] [    INFO] - Base service - Uvicorn server setup on port 9004
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9004 (Press CTRL+C to quit)
[2025-10-08 12:17:43,037] [    INFO] - Base service - HTTP server setup successful
[2025-10-08 12:17:43,038] [    INFO] - opea_service_omics_esminversefold_score_microservice - OPEA_OMICS_esminversefold_score server started.
```


###### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the esm inversefold score microservice.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

###### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.2:9004/v1/esminversefoldscore
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.2" # Replace with the actual server IP
export PORT="9004" # Replace with the actual server port
```

You only need to do this once per session.

###### Step 2: Run Inference Examples

Below are example of how to invoke the client for esm inversefolding score.

```bash
#example
export INPUT_PDB=../../../examples/inverse_folding/data/5YH2.pdb
export INPUT_FASTA=../../../examples/inverse_folding/data/5YH2_mutated_seqs.fasta 

python opea_scoring_sequence_client.py \
  --host $HOST \
  --port $PORT \
  --pdbfile $INPUT_PDB \
  --seqfile $INPUT_FASTA \
  --chain C
```
### LM Design

#### 1. Start the Server

Start the lm design microservice container:

```bash
docker run -it  \
  -v $MODELS:/checkpoints  \
  -v $INPUT:/input  \
  -v $OUTPUT:/output  \
  esm_image:latest  \
  python omics_setup/microservice/lm_design/opea_lmdesign_microservice.py --bf16
```

* This will launch the server on port 9005 and expose the inference API at `http://<host>:9005/v1/esmlmdesign`.

* To use a different port, specify it with the `--port` argument.
* `--bf16` flag accelerates performance by utilizing bfloat16 precision, which enhances computational efficiency without compromising accuracy.

Server Ready Check

The server is ready to accept requests when you see logs similar to the following:

```bash
$ docker run -it   -v $MODELS:/checkpoints   -v $INPUT:/input   -v $OUTPUT:/output   esm_image:latest   python omics_setup/microservice/lm_design/opea_lmdesign_microservice.py --bf16
/app/esm/omics_setup/examples/lm-design/lm_design.py:486: UserWarning:
The version_base parameter is not specified.
Please specify a compatability version level, or None.
Will assume defaults for version 1.1
  @hydra.main(config_path="conf/", config_name="config")
[2025-10-14 04:44:43,181] [    INFO] - esm_lmdesign_microservice - Info: Running in bfloat16 mode.
/app/esm/omics_setup/examples/lm-design/utils/struct_models.py:76: FutureWarning: You are using `torch.load` with `weights_only=False` (the current default value), which uses the default pickle module implicitly. It is possible to construct malicious pickle data which will execute arbitrary code during unpickling (See https://github.com/pytorch/pytorch/blob/main/SECURITY.md#untrusted-models for more details). In a future release, the default value for `weights_only` will be flipped to `True`. This limits the functions that could be executed during unpickling. Arbitrary objects will no longer be allowed to be loaded via this mode unless they are explicitly allowlisted by the user via `torch.serialization.add_safe_globals`. We recommend you start setting `weights_only=True` for any use case where you don't have full control of the loaded file. Please open an issue on GitHub for any issues related to this experimental feature.
  state = torch.load(local_model_path, map_location='cpu')
ESM model loaded
Model load time 8.506359100341797
[2025-10-14 04:44:51,687] [    INFO] - opea_service_omics_esmlmdesign_microservice - Starting esmlmdesign microservice at http://172.17.0.2:9005/v1/esmlmdesign
[2025-10-14 04:44:51,690] [    INFO] - Base service - CORS is enabled.
[2025-10-14 04:44:51,692] [    INFO] - Base service - Setting up HTTP server
[2025-10-14 04:44:51,693] [    INFO] - Base service - Uvicorn server setup on port 9005
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9005 (Press CTRL+C to quit)
[2025-10-14 04:44:51,706] [    INFO] - Base service - HTTP server setup successful
[2025-10-14 04:44:51,708] [    INFO] - opea_service_omics_esmlmdesign_microservice - OPEA_OMICS_esmlmdesign server started.
```


#### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the lm-design microservice.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

##### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.2:9005/v1/esmlmdesign
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.2" # Replace with the actual server IP
export PORT="9005" # Replace with the actual server port
```

You only need to do this once per session.

##### Step 2: Run Inference Examples

Below are example of how to invoke the client for lm-design.

Fixed backbone design
```bash
#example
export INPUT=../../examples/lm-design/2N2U.pdb

 python opea_lmdesign_client.py \
  --host $HOST \
  --port $PORT \
  task=fixedbb \
  pdb_fn=$INPUT
```

Free generation design
```bash
#example
export INPUT=../../examples/lm-design/2N2U.pdb

 python opea_lmdesign_client.py \
  --host $HOST \
  --port $PORT \
  task=free_generation 
```

