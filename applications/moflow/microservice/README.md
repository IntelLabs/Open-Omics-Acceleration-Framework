# Using Moflow as a Microservice

This directory provides a microservice interface for running Moflow. It allows you to expose Moflow as a server and interact with it via HTTP-based requests, enabling easy integration into larger workflows or client applications.

## Overview

The setup includes two main components:

1. Server – Hosts the Moflow model and exposes an API for inference.

2. Client – Sends requests to the server with required inputs and receives predictions as responses.

## Prerequisites

Before proceeding, ensure that you have already built the Moflow Docker image following the provided instructions and exported the required directories [here](../README.md#Building-the-Docker-Image).

## Running the Microservice Experiments

### Workflow 1: Experiments 1, 2, and 3

#### 1. Start the Server

Start the moflow microservice container:
##### Generate zinc250k

```bash
docker run -it --rm \
  -v $MODELS:/results \
  -v $DATA_PREPROCESSING:/data_preprocessing \
  -v $OUTPUT:/output \
  moflow:latest bash -c \
  "python microservice/generate/opea_moflow_microservice.py --model_dir /results/zinc250k_512t2cnn_256gnn_512-64lin_10flow_19fold_convlu2_38af-1-1mask -snapshot model_snapshot_epoch_200 --data_name zinc250k --data_dir /data_preprocessing  --hyperparams_path moflow-params.json"
```
🔒 Note: If you use `--data_name` qm9, you must also update the model directory in the `--model_dir` argument to `/results/qm9_64gnn_128-64lin_1-1mask_0d6noise_convlu1`.
* This will launch the server on port 9006 and expose the inference API at `http://<host>:9006/v1/moflow`.

* To use a different port, specify it with the `--port` argument.
* `--bf16` flag accelerates performance by utilizing bfloat16 precision, which enhances computational efficiency without compromising accuracy.


#### Server Ready Check

The server is ready to accept requests when you see logs similar to the following:

```bash
$ docker run -it -v $MODELS:/results -v $DATA_PREPROCESSING:/data_preprocessing -v $OUTPUT:/output moflow:latest bash -c "python microservice/generate/opea_moflow_microservice.py --model_dir /results/zinc250k_512t2cnn_256gnn_512-64lin_10flow_19fold_convlu2_38af-1-1mask -snapshot model_snapshot_epoch_200 --data_name zinc250k --data_dir /data_preprocessing --hyperparams_path moflow-params.json"
[2025-09-17 12:31:45,684] [    INFO] - esm_fold_microservice - Loading model
loading hyperparamaters from /results/zinc250k_512t2cnn_256gnn_512-64lin_10flow_19fold_convlu2_38af-1-1mask/moflow-params.json
loading snapshot: /results/zinc250k_512t2cnn_256gnn_512-64lin_10flow_19fold_convlu2_38af-1-1mask/model_snapshot_epoch_200
Hyper-parameters:
--------------------  ------------------------------------------------------------------------------------------------
b_n_type              4
b_n_flow              10
b_n_block             1
b_n_squeeze           19
b_hidden_ch           [512, 512]
b_affine              True
b_conv_lu             2
a_n_node              38
a_n_type              10
a_hidden_gnn          [256]
a_hidden_lin          [512, 64]
a_n_flow              38
a_n_block             1
mask_row_size_list    [1]
mask_row_stride_list  [1]
a_affine              True
path                  results/zinc250k_512t2cnn_256gnn_512-64lin_10flow_19fold_convlu2_38af-1-1mask/moflow-params.json
learn_dist            1
seed                  1
noise_scale           0.6
--------------------  ------------------------------------------------------------------------------------------------
model.ln_var: 1.12
loading train/valid split information from: data/valid_idx_zinc.json
Loading file /data_preprocessing/zinc250k_relgcn_kekulized_ggnp.npz
Model load time 48.51513934135437
[2025-09-17 12:32:34,200] [    INFO] - opea_service_omics_moflow_microservice - Starting moflow microservice at http://172.17.0.8:9006/v1/moflow
[2025-09-17 12:32:34,203] [    INFO] - Base service - CORS is enabled.
[2025-09-17 12:32:34,205] [    INFO] - Base service - Setting up HTTP server
[2025-09-17 12:32:34,206] [    INFO] - Base service - Uvicorn server setup on port 9006
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9006 (Press CTRL+C to quit)
[2025-09-17 12:32:34,223] [    INFO] - Base service - HTTP server setup successful
[2025-09-17 12:32:34,225] [    INFO] - opea_service_omics_moflow_microservice - OPEA_OMICS_moflow server started.
```


#### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the moflow microservice to run various tasks like generate and optimize.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

##### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.8:9006/v1/moflow
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.8" # Replace with the actual server IP
export PORT="9006" # Replace with the actual server port
```

You only need to do this once per session.

##### Step 2: Run Inference Examples

Below are examples of how to invoke the client for different moflow use-cases.

```bash
cd applications/moflow/microservice/
```

###### 1. Experiment: reconstruction

```bash
#example
python generate/opea_moflow_client.py \
    --host $HOST \
    --port $PORT \
    --batch_size 256 --reconstruct
```

###### 2. Experiment: Random generation
```bash
#example
python generate/opea_moflow_client.py \
  --host $HOST \
  --port $PORT \
  --n_experiments 5 --correct_validity true --batch_size 10000 --temperature 0.85 --delta 0.05 --save_fig false --correct_validity true --random_generation
```

###### 3. Experiment: Interpolation generation & visualization

interpolation between 2 molecules (molecular graphs)
```bash
#example
python generate/opea_moflow_client.py \
  --host $HOST \
  --port $PORT \
  --batch_size 1000  --temperature 0.65   --int2point --inter_times 50  --correct_validity true
```

interpolation in a grid of molecules (molecular graphs)
```bash
#example
python generate/opea_moflow_client.py \
  --host $HOST \
  --port $PORT \
  --batch_size 1000  --temperature 0.65 --delta 5  --intgrid  --inter_times 40  --correct_validity true
```


### Workflow 2: Experiment 4

#### 1. Start the Server

Start the moflow microservice container:
##### Optimize zinc250k

```bash
docker run -it --rm \
  -v $MODELS:/results \
  -v $DATA_PREPROCESSING:/data_preprocessing \
  -v $OUTPUT:/output \
  moflow:latest bash -c \
  "python microservice/optimize/opea_moflow_optimize_microservice.py -snapshot model_snapshot_epoch_200  --hyperparams_path moflow-params.json --batch_size 256 --model_dir /results/zinc250k_512t2cnn_256gnn_512-64lin_10flow_19fold_convlu2_38af-1-1mask  --data_name zinc250k   --property_name qed --property_model_path qed_model.pt --debug false --data_dir /data_preprocessing"
```

🔒 Note: If you use `--data_name` qm9, you must also update the model directory in the `--model_dir` argument to `/results/qm9_64gnn_128-64lin_1-1mask_0d6noise_convlu1`.
* This will launch the server on port 9007 and expose the inference API at `http://<host>:9007/v1/moflow`.

* To use a different port, specify it with the `--port` argument.
* `--bf16` flag accelerates performance by utilizing bfloat16 precision, which enhances computational efficiency without compromising accuracy.

#### Server Ready Check

The server is ready to accept requests when you see logs similar to the following:

```bash
$ docker run -it   -v $MOD   -v $DATA_PREPROCESSING:/data_preprocessing   -v $OUTPUT:/output   moflow:latest bash -c   "cd .. && python microservice/optimize/opea_moflow_optimize_microservice.py -snapshot model_snapshot_epoch_200  --hyperparams_path moflow-params.json --batch_size 256 --model_dir /results/zinc250k_512t2cnn_256gnn_512-64lin_10flow_19fold_convlu2_38af-1-1mask  --data_name zinc250k   --property_name qed --property_model_path qed_model.pt --debug false --data_dir /data_preprocessing"
all import is done
25/09/17 12:52:53 | INFO | root | Loading model
Start at Time: Wed Sep 17 12:52:53 2025
loading snapshot: /results/zinc250k_512t2cnn_256gnn_512-64lin_10flow_19fold_convlu2_38af-1-1mask/model_snapshot_epoch_200
Hyper-parameters:
--------------------  ------------------------------------------------------------------------------------------------
b_n_type              4
b_n_flow              10
b_n_block             1
b_n_squeeze           19
b_hidden_ch           [512, 512]
b_affine              True
b_conv_lu             2
a_n_node              38
a_n_type              10
a_hidden_gnn          [256]
a_hidden_lin          [512, 64]
a_n_flow              38
a_n_block             1
mask_row_size_list    [1]
mask_row_stride_list  [1]
a_affine              True
path                  results/zinc250k_512t2cnn_256gnn_512-64lin_10flow_19fold_convlu2_38af-1-1mask/moflow-params.json
learn_dist            1
seed                  1
noise_scale           0.6
--------------------  ------------------------------------------------------------------------------------------------
loading train/valid split information from: data/valid_idx_zinc.json
Loading file /data_preprocessing/zinc250k_relgcn_kekulized_ggnp.npz
Load /data_preprocessing/zinc250k_relgcn_kekulized_ggnp.npz done, length: 249455
Loading trained regression model for optimization
Load data/zinc250k_property.csv done, length: 249455
Prepare data done! Time 49.60 seconds
loading qed regression model from: /results/zinc250k_512t2cnn_256gnn_512-64lin_10flow_19fold_convlu2_38af-1-1mask/qed_model.pt
Load model done! Time 49.97 seconds
Model load time 50.01725721359253
[2025-09-17 12:53:43,204] [    INFO] - opea_service_omics_moflow_optimize_microservice - Starting moflow microservice at http://172.17.0.7:9007/v1/moflow_optimize
[2025-09-17 12:53:43,207] [    INFO] - Base service - CORS is enabled.
[2025-09-17 12:53:43,209] [    INFO] - Base service - Setting up HTTP server
[2025-09-17 12:53:43,210] [    INFO] - Base service - Uvicorn server setup on port 9007
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9007 (Press CTRL+C to quit)
[2025-09-17 12:53:43,224] [    INFO] - Base service - HTTP server setup successful
[2025-09-17 12:53:43,225] [    INFO] - opea_service_omics_moflow_optimize_microservice - OPEA_OMICS_moflow server started.
```

#### 2. Query the Server from a Client Application

We provide a sample client program that interacts with the moflow microservice to run various tasks like generate and optimize.

The client can be run:

* On the same machine as the server, or

* On a different machine within the same network.

##### Step 1: Set the Server Address

After starting the server, note the API endpoint from the logs:

```bash
🔗 API Endpoint: http://172.17.0.7:9007/v1/moflow_optimize
```

On the client machine, export the following environment variables:

```bash
export HOST="172.17.0.7" # Replace with the actual server IP
export PORT="9007" # Replace with the actual server port
```

You only need to do this once per session.

#### Step 2: Run Inference Examples

Below are examples of how to invoke the client for different moflow use-cases.


##### 4. Experiment: Molecular optimization & constrained optimization

###### To optimize existing molecules to get novel molecules with optimized QED scores

```bash
#example
python optimize/opea_moflow_optimize_client.py \
  --host $HOST \
  --port $PORT \
  --topk 5 --topscore
```
###### To optimize existing molecules to get novel molecules with optimized plogp scores and constrained similarity

```bash
#example
python optimize/opea_moflow_optimize_client.py \
  --host $HOST \
  --port $PORT \
  --topk 20  --consopt  --sim_cutoff 0
```
