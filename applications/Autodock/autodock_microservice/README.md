# AutoDock-GPU OPEA Microservice
Production-ready microservice for running **AutoDock-GPU (CPU backend)** through a **FastAPI + OPEA framework**, with:
* Secure multi-user isolated workspaces
* Base64 + .tar.gz dataset upload
* Full server-side validation
* Clean job success / failure separation
* Minimal output return (.dlg, .xml, metadata.json)
* Automatic workspace cleanup
* Standardized HTTP error schema
# Docker Setup Instructions
## 1. Build the Docker Image
```zsh
cd build_docker/
docker build -t autodock_microservice_scan .
cd ..
```
**Note:Change in dockerfile `NUMWI=64` to `NUMWI=128` if you want 128 bit size of autodock executable**
## 2. Run the Microservice (Dynamic Port)
You can run on any free port:
```bash
docker run --rm -it -p 9012:9012 -v $(pwd):/opt/AutoDock/microservice autodock_microservice_scan:latest python /opt/AutoDock/microservice/server.py --port 9012
```
Service will be available at:
```bash
http://127.0.0.1:9012/v1/autodock
```
# Client Usage
## 1. Dataset Folder Example
```bash
1mzc/
├── protein.maps.fld
├── protein.maps.e.map
├── protein.maps.d.map
├── rand-0.pdbqt
├── rand-5.pdbqt
└── ...
```
## 2. Run a Docking Job
```bash
python opea_autodock_client.py   --host 127.0.0.1   --port 9012   --user-id testuser   --dataset-dir ./5wlo/   --ffile protein.maps.fld   --lfile rand-0.pdbqt   --nrun 10   --lsmet sw   --resnam rand-5 --seed 11,89,23 --nev 10000
```
What Happens Internally
### 1. Client creates:
```bash
1mzc.tar.gz → base64 → JSON request
```
### 2. Server extracts into:
```bash
/tmp/autodock_inputs/testuser/default/<uuid>/1mzc
```
### 3. AutoDock is run inside 1mzc/
### 4. Server creates:
* metadata.json
* rand-5.dlg
* rand-5.xml
### 5. Server returns:
* metadata_b64
* results_tar_b64 (only those 3 files)
### 6.Workspace is deleted automatically
