# OpenOmics GROMACS Microservice (Workflow Mode Only)
This microservice allows clients to run a **full GROMACS workflow (EM → NVT → NPT → MD01)** using a single API call.

The client sends a **PDB file as base64**, and the service:
* Executes a workflow script (`run_commands.sh`) inside an isolated job workspace
* Returns a `metadata JSON` artifact only (no heavy `.xtc / .trr` files)
* Automatically deletes the workspace after each run
## Workflow-only execution using:
```bash
./run_commands.sh <pdb_filename>
```
* Per-user, per-job isolated directory under `/input/jobs/`
* Sanitized `user_id` / `job_label` (safe against traversal)
* Base64 PDB upload
* No large outputs returned (.xtc/.trr). Only:
* metadata_<job_id>.json (base64 artifact)
* Auto cleanup of workspace after every run
* stdout/stderr snippets included in metadata
Dockerfile includes `/input` whose contents are:
```bash
run_commands.sh
mdtut_minim.mdp
mdtut_nvt.mdp
mdtut_npt.mdp
mdtut_md.mdp
```

# Features
Build image:
```bash
cd build_docker
sudo docker build -t gromacs_microservice_scan:latest .
cd ..
```

# Running the Server:                                                      
```bash
docker run --rm -it -p 9012:9012   -v $(pwd):/microservice  gromacs_microservice_scan:latest   /opt/conda/envs/grms_env/bin/python /microservice/opea_gromacs_microservice.py --port 9012
```
* -p : Maps container port 9012 → host port 9012 so the API is reachable
* `gromacs_opea:latest` : Name of the Docker image to run
* --port : Microservice listens on port 9012 inside the container
* `$(pwd)/gromacs_microservice:/microservice` : Mounts your local microservice code into the container

# Client Example
```bash
python opea_gromacs_client.py     --host 127.0.0.1     --port 9012     --user-id user     --job-label test1     --pdb-file ./1AKI.pdb
```
This will:
* Encode the PDB
* Submit a workflow job
* Save the returned `metadata_<jobid>.json`

