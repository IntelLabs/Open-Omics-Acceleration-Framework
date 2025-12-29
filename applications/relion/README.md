# Open-Omics-Relion
**Open-Omics-Relion** is a Dockerized RELION 5.0 setup for running benchmark workloads using Intel oneAPI and Intel MPI. It supports **2D classification**, **3D classification**, and **auto-refinement** modes using official test data, designed for **reproducibility**, **performance testing**, and **ease of deployment**.

---

## Step 1: Download Benchmark Dataset
```zsh
wget ftp://ftp.mrc-lmb.cam.ac.uk/pub/scheres/relion_benchmark.tar.gz
tar -xzvf relion_benchmark.tar.gz
```
This will extract a folder named `relion_benchmark`.
## Step 2: Build the Docker Image
Build the docker image with `Dockerfile`, run:
```zsh
cd build_docker/
docker build -t relion_nru .
cd ..
```
Verify the image was built:
```zsh
docker images | grep -i relion_nru
```
## Step 3: Change Ownership of Benchmark Data
To avoid permission issues when mounting the directory (RELION runs as a non-root user UID 1001):
```zsh
cd relion_benchmark/
sudo chown 1001:1001 $(pwd)
```

## Step 4: Run RELION Benchmark
Launch the container with required options for shared memory and MPI support:
```zsh
sudo docker run --rm --net=host --ipc=host --pid=host --ulimit stack=67108864 --shm-size=2g --cap-add=SYS_PTRACE -e I_MPI_DEBUG=5 -e I_MPI_SHM_LMT=shm -e I_MPI_FABRICS=shm:tcp -it -v $(pwd):/opt/relion_5.0/relion_benchmark relion_nru:latest
```
**Notes**
- Modify entrypoint.sh if you wish to customize CPU/thread usage.
- The container is designed to work on systems with Intel MPI and oneAPI properly set up inside.
- Ensure Docker has access to sufficient system shared memory (e.g. via /dev/shm or --shm-size if needed).
- The container uses a non-root user (UID 1001). Make sure mounted volumes are writable by this user.
### Available Run Modes
You can specify different modes as an argument to the Docker command. Each mode will automatically create an output folder inside the `relion_benchmark` directory.

| Mode         | Command Example (append to `docker run`)                    | Description                  | Output Folder |
|--------------|--------------------------------------------------------------|------------------------------|--------------|
| *(default)*   | `relion_nru:latest`               | Run 3D classification        |     `3D/`      |
| `3d`          | `relion_nru:latest 3d`               | Run 3D classification        |     `3D/`    |
| `2d`         | `relion_nru:latest 2d`                                      | Run 2D classification        |   `2D/`   |
| `autorefine` | `relion_nru:latest autorefine`                              | Run 3D auto-refinement       |  `3D_AUTO/`   |
| `custom`     | `relion_nru:latest relion_refine_mpi [your-flags]`          | Run any custom RELION command |  User-defined (`--o`) |

---

**The Original README for Relion starts here:**


RELION 5.0-beta
===============

RELION (for REgularised LIkelihood OptimisatioN) is a stand-alone computer
program for Maximum A Posteriori refinement of (multiple) 3D reconstructions
or 2D class averages in cryo-electron microscopy. It is developed in the
research group of Sjors Scheres at the MRC Laboratory of Molecular Biology.

If RELION is useful in your work, please cite our papers.

Comprehensive documentation of RELION and tutorials are stored [here](https://relion.readthedocs.io/).

## Installation

See our [installation instructions](https://relion.readthedocs.io/en/release-5.0/Installation.html).

You will have to set up a Python environment to use Python modules (e.g. Blush, ModelAngelo and DynaMight).
Thus, please read the above instructions carefully even if you are familiar with earlier versions.

## Class Ranker

The default model for the class ranker has been trained and tested in Python 3.9.12 with Pytorch 1.10.0 and Numpy 1.20.0.
If you wish to retrain the class ranker model with your own data, please refer to [this repo](https://github.com/3dem/relion-classranker).
