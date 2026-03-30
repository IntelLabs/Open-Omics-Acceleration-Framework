<div align="center">  

[Paper](https://doi.org/10.1101/2024.11.19.624167) |
[Slack](https://join.slack.com/t/boltz-community/shared_invite/zt-2zj7e077b-D1R9S3JVOolhv_NaMELgjQ) <br> <br>
</div>

## Introduction

Boltz-1 is the state-of-the-art open-source model to predict biomolecular structures containing combinations of proteins, RNA, DNA, and other molecules. It also supports modified residues, covalent ligands and glycans, as well as conditioning the prediction on specified interaction pockets or contacts. 

All the code and weights are provided under MIT license, making them freely available for both academic and commercial uses. For more information about the model, see our [technical report](https://doi.org/10.1101/2024.11.19.624167). To discuss updates, tools and applications join our [Slack channel](https://join.slack.com/t/boltz-community/shared_invite/zt-2zj7e077b-D1R9S3JVOolhv_NaMELgjQ).

---

## 🔍 Running Boltz with Docker (CPU Only)

This repository provides a unified Docker image that can run in two modes:
1.  **Batch/CLI Mode:** Process local folders containing FASTA/YAML files automatically.
2.  **Microservice Mode:** Start an OPEA-compliant REST API server for on-demand predictions.

---

### 📁 1. Preparation (Important)

Before running in either mode, create your directories and **ensure write permissions**. Since the container runs as a non-root user, it needs explicit permission to write to your host folders.

```bash
# Create folders
mkdir -p inputs outputs models

# Grant write permissions (Required for Docker user)
chmod 777 inputs outputs models
          
# Set convenience variables
export INPUT=$PWD/inputs
export OUTPUT=$PWD/outputs
export MODELS=$PWD/models
```

> ⚠️ **Note:** Place your `.fasta` or `.yaml` input files inside the `inputs` folder before running Batch Mode.

---

### 🐳 2. Build the Docker Image

To build the boltz docker image go to `applications` folder. Copy `common` folder into `boltz` folder. Then go inside `boltz`. This will become our context path. Then run:

```bash
docker build -f Dockerfile --network=host -t boltz:latest .
```

---

### 🚀 Mode A: File Processing (CLI)

Use this mode to automatically process all files in your input directory. 
Please use boltz repo for instructions on how to format the input files. Addtionally, you can use --use_msa_server as the argument for MSA generation.

Run the container mounting the volumes

```bash
docker run --rm \
  --ipc=host \
  --shm-size=100g \
  -v $INPUT:/inputs \
  -v $OUTPUT:/outputs \
  -v $MODELS:/app/.boltz_cache \
  boltz:latest \
  boltz predict /inputs/test.fasta \
  --out_dir /outputs \
  --accelerator cpu
```

**What happens:**
1. The container scans `$INPUT` for `.fasta` or `.yaml` files.
2. It runs Boltz inference on CPU.
3. Results are saved to `$OUTPUT`.
4. The container exits automatically when finished.

### multiprocess
In order to run using multiprocess

```bash
docker run -it --rm \
   --user root \
   --ipc=host \
   --privileged \
   --shm-size=100g \
   -v $INPUT:/inputs \
   -v $OUTPUT:/outputs \
   -v $MODELS:/app/.boltz_cache \
   -v $PWD/multiprocess_config.json:/app/boltz/multiprocess_config.json \
   boltz:latest \
   python common/multiprocess/multiprocess.py --json_file=multiprocess_config.json --case=2
```

Sample config looks like this [multiprocess_config.json](multiprocess_config.json).
For every file add `"/inputs/<file.name> --accelerator cpu --override --out_dir /outputs"` under `unique_args`

---

### 🌐 Mode B: Microservice (API Server)

Use this mode to keep the server running and send requests programmatically (via Python/Curl).

Start the container with the `microservice` argument:

```bash
docker run -it --rm --network=host boltz:latest python /app/boltz/microservice/boltz_microservice_server.py
```

#### Checking Status
*   **Logs:** `docker logs -f boltz-server`
*   **API Docs:** Open `http://localhost:8000/docs` in your browser.

#### Sending a Request
The server accepts a JSON payload and returns a Base64 encoded ZIP file containing the results.

**Python Client Example:**

```python
import requests
import base64

url = "http://localhost:8000/v1/boltz"

# 1. Read your YAML input
with open("my_input.yaml", "r") as f:
    payload = {"yaml_content": f.read()}

# 2. Send POST request
response = requests.post(url, json=payload)
data = response.json()

if data.get("status") == "success":
    # 3. Decode Base64 result to ZIP
    zip_content = base64.b64decode(data["result"])
    with open("results.zip", "wb") as f:
        f.write(zip_content)
    print("✅ Saved results.zip")
else:
    print(f"❌ Error: {data.get('message')}")
```

There is one complete python [Click here to see the Client Code](./microservice/boltz_microservice_client.py)
---

## License

Our model and code are released under MIT License, and can be freely used for both academic and commercial purposes.

## Cite

If you use this code or the models in your research, please cite the following paper:

```bibtex
@article{wohlwend2024boltz1,
  author = {Wohlwend, Jeremy and Corso, Gabriele and Passaro, Saro and Reveiz, Mateo and Leidal, Ken and Swiderski, Wojtek and Portnoi, Tally and Chinn, Itamar and Silterra, Jacob and Jaakkola, Tommi and Barzilay, Regina},
  title = {Boltz-1: Democratizing Biomolecular Interaction Modeling},
  year = {2024},
  doi = {10.1101/2024.11.19.624167},
  journal = {bioRxiv}
}
```

In addition if you use the automatic MSA generation, please cite:

```bibtex
@article{mirdita2022colabfold,
  title={ColabFold: making protein folding accessible to all},
  author={Mirdita, Milot and Sch{\"u}tze, Konstantin and Moriwaki, Yoshitaka and Heo, Lim and Ovchinnikov, Sergey and Steinegger, Martin},
  journal={Nature methods},
  year={2022},
}
```
