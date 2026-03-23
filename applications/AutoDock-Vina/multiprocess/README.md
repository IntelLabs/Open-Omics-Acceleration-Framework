# Multiprocess Task Runner – AutoDock Vina Integration
## Overview
This repo provides a **generic multiprocess launcher** [multiprocess.py](https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/applications/common/multiprocess/multiprocess.py) plus **AutoDock Vina wrapper scripts** that let you run many docking jobs in parallel with:

* automatic process count selection (based on cores/NUMA),
* NUMA memory binding (numactl -m) and CPU pinning (numactl -C),
* per-task logs under a timestamped log folder.

The runner reads a JSON config and supports three common execution patterns:
* Case 1: one complex, many ligands (auto-discovered from `input_dir`)
* Case 2: multiple complexes / arbitrary ligand list (explicit `unique_args`)
* Case 3: one complex, many ligands, but each ligand can have different parameters (e.g., exhaustiveness/seed)
## Prerequisites
* Linux node with `numactl`, `lscpu`
* Python packages:
   * `psutil`
   * `absl-py`
## Build Docker
```bash
cd ../build_docker/
docker build -t autodock_vina_multiprocess:latest --build-arg FLAVOR=multiprocess .
cd ../multiprocess/
```

## Output and Logging
* Each run creates a timestamped log directory under `output_log_dir` (e.g., `vina_logs/2025_12_15_074501_/`)
* Each task writes a separate log: <unique_id>.txt (`stdout+stderr`)
* Each wrapper script writes docking output files under the complex directory:
```bash
results/<ligand>_out.pdbqt
```

## Expected Data Layout
* Each complex directory contains:
* ligand files: `ligands/*.pdbqt`
* per-ligand Vina config files named like: `vina_<ligand_basename>_config.txt`
Example:
```bash
data/
  1ac8/
    vina_rand-0_config.txt
    vina_rand-1_config.txt
    ligands/
      rand-0.pdbqt
      rand-1.pdbqt
    results/   (created automatically)
```
The wrapper scripts assume the config file name is derived from the ligand name:
```bash
vina_${lig_base}_config.txt
```

## How CPU threads are chosen
* `multiprocess.py` sets:
```bash
OMP_NUM_THREADS = total_cores / max_processes
```
  * Each wrapper uses:
     ```bash
       threads="${OMP_NUM_THREADS:-1}"
     ```
  * passes this to Vina:
     ```bash
	 vina --cpu "${threads}"
	 ```
This ensures each Vina process uses only the CPU cores assigned to it.
## Case 1: Single complex, run all ligands in `input_dir`

**Wrapper**: `run_vina_case1_out.sh`

**Inputs**:
    
	*  `complex_dir=<complex_path>`
	
    *  `input_file=<ligand_path>`

**Behaviour**:

* cd <complex_dir>
* finds vina_<lig_base>_config.txt
* writes output to results/<lig_base>_out.pdbqt
**Run**:
  
```bash
cp ../../common/multiprocess/multiprocess.py .
sudo docker run --rm --privileged -u 0  -v $(pwd)/run_vina_case1_out.sh:/work/run_vina_case1_out.sh -v $(pwd)/multiprocess.py:/work/multiprocess.py -v $(pwd)/vina_case1.json:/work/vina_case1.json -v $(pwd)/data:/work/data -it autodock_vina_multiprocess:latest bash
```
Inside docker run
```bash
chmod +x run_vina_case1_out.sh
python multiprocess.py --json_file=vina_case1.json --case=1
```
## Case 2: Multi-protein / arbitrary ligand list (explicit tasks)

**Wrapper**:`run_vina_case2_out.sh`

**Run**:

```bash
python multiprocess.py --json_file=vina_multi_protein_case2.json --case=2
```
## Case 3: Single complex + per-task unique parameters (seed/exhaustiveness)

**Wrapper**:`new_script_case3.sh`

Supports optional parameters:

 * exhaustiveness=<int>
 * seed=<int>
 
**Run**:

```bash
python multiprocess.py --json_file=vina_case3_new.json --case=3
```

