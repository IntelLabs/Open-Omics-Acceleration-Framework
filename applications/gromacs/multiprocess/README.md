# GROMACS Multiprocess Runner 
This setup runs multiple full GROMACS pipelines in parallel on a single node using:
* CPU pinning (`numactl`)
* per-protein isolated working directories
* automatic MPI sizing (`ntmpi = pinned CPUs`)
* `OMP_NUM_THREADS=1` to avoid oversubscription
# Prerequisites
Docker image: `gromacs_opea_multi_omp1`
* Files mounted into the container:
* [multiprocess.py](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/applications/common/multiprocess/multiprocess.py)
* `new_script.sh` (wrapper)
* `run_commands.sh`
* `mdtut_*.mdp`
* protein PDBs
## Case 1: Auto-discover proteins
* All `.pdb` files in a folder are treated as independent GROMACS jobs
* One protein → one isolated working directory
* Best when all proteins share the same workflow
Run:

```bash
cp ../../common/multiprocess/multiprocess.py .
sudo docker run --rm --privileged -u 0 -v $(pwd)/run_gromacs_multi_one.sh:/input/run_gromacs_multi_one.sh -v $(pwd)/multiprocess.py:/input/multiprocess.py -v $(pwd)/gromacs_case1.json:/input/gromacs_case1.json -v $(pwd)/proteins:/input/proteins -v $(pwd)/output:/tmp/work_gromacs -it gromacs_opea_multi_omp1:latest bash
```
Inside container:
```bash
/opt/conda/envs/grms_env/bin/python multiprocess.py --json_file=gromacs_case1.json --case=1
```
## Outputs
Logs:
```bash
gmx_logs/<timestamp>/<protein>.txt
```

Per-protein results:
```bash
/tmp/work_gromacs/<protein>_<timestamp>/
```
## case 2: Explicit protein list (manual control)

Still:
* one protein = one isolated directory
* same `run_commands.sh`
* same wrapper (`run_gromacs_multi_one.sh`)

Useful when:
* proteins are in different folders
* you want a specific subset
* you want ordering control
Run:
```bash
/opt/conda/envs/grms_env/bin/python multiprocess.py --json_file=gromacs_case2.json --case=2
```

