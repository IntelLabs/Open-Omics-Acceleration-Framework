# RELION OPEA Microservice

This repository provides a production-grade FastAPI microservice for running
RELION 2D classification, 3D classification, and 3D auto-refinement jobs in
an isolated per-request workspace with S3-based input/output.

## Features
- REST API with FastAPI + Pydantic
- Per-user isolated workspaces
- Input via S3 (.tar.gz) (no large file transfer through API)
- Output directory upload to S3
- Strict MPI / OMP / pool validation
- CPU-only execution
- Automatic workspace cleanup
- Structured metadata.json results


## Environment Requirements

-  Python 3.9+
-  RELION 5.0 (CPU build)
-  OpenMPI / Intel MPI
-  AWS CLI
-  Valid IAM Role or AWS credentials with:
  * `s3:GetObject`
  * `s3:PutObject`

## Build Dockerfile
```bash
cd build_docker/
docker build --no-cache -t relion_microservice_scan:latest .
cd .. 
##  Configurable RELION Binary

RELION path is configurable through environment:

```bash
export RELION_EXE=/opt/relion_5.0_cpu_benchmark_prefix/bin/relion_refine_mpi
```

## Workflow
S3 (.tar.gz dataset) -> Per-user workspace (/tmp/workspaces/<user>/<job_id>) -> relion_refine_mpi (MPI + OpenMP) -> Output directory detected from --o path -> Upload output(.tar.gz) to S3 -> Workspace automatically deleted -> Client receives metadata.json 
## How to Run
1. Start the server:
```bash
sudo docker run --rm --ipc=host -it -p 9012:9012 -v "$(pwd)":/microservice relion_microservice_scan:latest service --port 9012
```
Server will be available at:
```bash
http://<HOST_IP>:9012/v1/relion
```

2. Run the client(Autorefine Example):
```bash
python relion_opea_client.py   --host 127.0.0.1   --port 9012   --user-id srilekha   --workspace-id job_refine   --dataset-s3-url s3://srilekha-new-s3/data_3d_refine.tar.gz   --mode autorefine   --num-mpi 7 --i particles.star --ref run_it025_class002_box256.mrc --o Refine/rf  --j 24 --pool 48 --save-metadata relion_refine_ok.json --healpix-order 2 --ini-high 50 --auto-refine --split-random-halves  --auto-local-healpix-order 4 --low-resol-join-halves 40 --sym D2 --offset-range 5 --offset-step 2 --norm --scale --zero-mask --oversampling 1  --auto-ignore-angles --auto-resol-angles --ctf --ctf --particle-diameter 200 --flatten-solvent --dont-combine-weights-via-disc --preread-images --firstiter-cc --output-s3-prefix s3://srilekha-new-s3/relion_outputs
```

## API Contract
Input (RelionInput)

* `user_id`
* `workspace_id (optional)`
* `dataset_s3_url`
* `output_s3_prefix (optional)`
* `mode`: 2d, 3d, autorefine
* `num_mpi`
* `flags (structured RELION flags)`
Output (RelionOutput)
```bash
{
  "status": "success",
  "message": "...",
  "metadata_b64": "...",
  "job_id": "job_2d_001_abcd1234"
}
```
**Note: This microservice has been tested on the RELION benchmark dataset (`ftp://ftp.mrc-lmb.cam.ac.uk/pub/scheres/relion_benchmark.tar.gz`) and the EMPIAR-10204 tutorial dataset (`ftp://ftp.mrc-lmb.cam.ac.uk/pub/scheres/relion50_tutorial_precalculated_results.tar.gz`).**
