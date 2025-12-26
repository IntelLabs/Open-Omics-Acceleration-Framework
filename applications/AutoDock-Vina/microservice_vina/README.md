# AutoDock Vina OPEA Microservice

This microservice exposes **AutoDock Vina docking** as a REST API using the OPEA framework.  
It supports **multi-user workspaces**, **secure dataset uploads**, **validated inputs**, and **automatic output tarball download with metadata**.

---

## 1. Features

- **Multi-user & multi-workspace isolation**
  - Each request runs under `/workspaces/<user_id>/<workspace_id>/<job_id>/`.
- **Secure `tar.gz` dataset upload**
  - Client sends a base64-encoded tar.gz.
  - Server safely extracts it (no absolute paths, no `..` path traversal).
- **Validated inputs via Pydantic**
  - All key Vina flags (`--receptor`, `--ligand`, `--center_x`, `--size_x`, etc.) are validated on the server.
  - `scoring` is restricted to `"vina"` only.
- **Automatic output tarball with metadata**
  - If requested, the server returns a base64-encoded tar.gz containing:
    - The Vina output PDBQT file (e.g. `rand-1_out.pdbqt`)
    - A `metadata.json` with `user_id`, `workspace_id`, `job_id`, runtime, and the exact Vina command.
- **Workspace cleanup**
  - Job directories under `/workspaces/<user_id>/<workspace_id>/<job_id>/` are removed after each request when `--cleanup-on-success` is enabled.
- **Structured error handling**
  - Uses FastAPI + `HTTPException` for proper HTTP errors (e.g. missing ligand file, invalid flags).
  - Vina runtime failures (e.g. invalid `--seed`) return HTTP 4xx with detailed JSON (`stdout`/`stderr`).

---

## 2. Docker Image

This microservice runs fully inside a Docker container.  
The Dockerfile installs AutoDock Vina, `opea-comps`, and the Python runtime.

From the project’s `build_docker/` directory:

```bash
cd build_docker/
docker build -t autodock_vina_microservice:latest --build-arg FLAVOR=microservice
cd ..
``` 
After this step, the image `autodock_vina_microservice:latest` will be available locally.
You can verify with:
```bash
docker images | grep autodock_vina_microservice
```
## 3. Running the Server (Docker):
### 3.1 Fix permissions (because container user is non-root)
From the host machine, run:
```bash
chown -R 1001:1001 .
```
This ensures the non-root user inside the container (UID 1001) can read/write the code and workspace directories.
### 3.2 Start the server
Then run command:
```bash
docker run --rm -it   -v "$(pwd):/microservice"   -v "$(pwd)/workspaces:/workspaces"   -p 9012:9012   autodock_vina_microservice:latest   python /microservice/vina_opea_server.py     --port 9012     --workspace-root /workspaces     --max-tar-size-mb 500
```
Arguments:
* --port: Port for the HTTP microservice (default: 9000).
* -p : `HOST_PORT` : `CONTAINER_PORT`
* --workspace-root: Root directory for workspaces inside the container (default: /workspaces).
* --max-tar-size-mb: Max allowed size (MB) of decoded dataset tar.gz (default: 500 MB).
* --cleanup-on-success: Every job directory is removed after the request completes (whether Vina succeeded, failed, or raised HTTPException).(This deletes `/workspaces/<user>/<workspace>/<job>/`)
startup the server logs something like:
```bash
[2025-12-04 09:39:52,264] [    INFO] - autodock_vina_microservice_server - Starting AutoDock Vina microservice at http://172.17.0.2:9012/v1/vina
[2025-12-04 09:39:52,265] [    INFO] - autodock_vina_microservice - OpeaAutodockVina initialized: workspace_root=/workspaces, max_tar_bytes=524288000, cleanup_on_success=True
[2025-12-04 09:39:52,266] [    INFO] - Base service - CORS is enabled.
[2025-12-04 09:39:52,267] [    INFO] - Base service - Setting up HTTP server
[2025-12-04 09:39:52,267] [    INFO] - Base service - Uvicorn server setup on port 9012
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:9012 (Press CTRL+C to quit)
[2025-12-04 09:39:52,282] [    INFO] - Base service - HTTP server setup successful
```
## 4. Run the Client
```bash
python vina_opea_client.py   --host 127.0.0.1   --port 9012   --user-id sri   --dataset-dir 5wlo   --receptor protein.pdbqt   --ligand rand-0.pdbqt   --center-x 16.459   --center-y -19.946   --center-z -5.850   --size-x 18.375   --size-y 18.375   --size-z 18.375 --cpu 4 --exhaustiveness 4 --download-tarball --tarball-output vina_result.tar.gz
```
Where:
* `--dataset-dir 5wlo` points to a local folder like:
```bash
5wlo/
  protein.pdbqt
  flex-xray.pdbqt
  rand-1.pdbqt
  rand-0.pdbqt
  rand-2.pdbqt
  protein.A.map
  protein.C.map
  ...
```
You can run the server and client on **different machines** as long as:
* The client machine can reach the server machine’s public IP and port (no firewall block).
## Notes
* Only `"vina"` scoring is supported (`--scoring` is fixed to `vina`).
* This microservice currently supports a **single receptor – single ligand** use case per request.

