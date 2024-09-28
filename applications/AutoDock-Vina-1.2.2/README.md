#AutoDock Vina: Docking and virtual screening program

**AutoDock Vina** is one of the **fastest** and **most widely used** **open-source** docking engines. It is a turnkey computational docking program that is based on a simple scoring function and rapid gradient-optimization conformational search. It was originally designed and implemented by Dr. Oleg Trott in the Molecular Graphics Lab, and it is now being maintained and develop by the Forli Lab at The Scripps Research Institute.

* AutoDock4.2 and Vina scoring functions
* Support of simultaneous docking of multiple ligands and batch mode for virtual screening
* Support of macrocycle molecules
* Hydrated docking protocol
* Can write and load external AutoDock maps
* Python bindings for Python 3

## AutoDock Vina 1.2.2 Docker Setup

This repository contains a Docker setup for running AutoDock Vina 1.2.2. The following instructions guide you through cloning therepository, building the Docker image, setting up input/output directories, and running AutoDock Vina inside the Docker container.
### 1. Clone the Repository

First, clone the repository containing AutoDock Vina 1.2.2:

```bash
git clone https://github.com/intel-sandbox/TransOmics.OpenOmicsInternal/tree/main/applications/AutoDock-Vina-1.2.2
cd AutoDock-Vina-1.2.2
```bash

### 2. Build the Docker Image
Build the Docker image using the following command:
```bash
docker build -t docker_vina .
```bash

### 3. Setup Input and Output Directories
Create the input and output directories:
```bash
mkdir input_local
mkdir output_local
```bash
### 4. Prepare Input Files
Add your receptor (.pdbqt), ligand (.pdbqt), and dependent map files to the input_local directory

### 5. Set Environment Variables
Set environment variables to reference the input and output directories:
```bash
export INPUT_VINA=$PWD/input_local
export OUTPUT_VINA=$PWD/output_local
```bash
### 6. Run the Docker Container
Run the Docker container and execute AutoDock Vina:
```bash
docker run -it -v $INPUT_VINA:/input -v $OUTPUT_VINA:/output docker_vina sh -c "cd /input && vina --receptor protein.pdbqt --ligand rand-1.pdbqt --out /output/rand-1_out.pdbqt --center_x 16.459 --center_y -19.946 --center_z -5.850 --size_x 18 --size_y 18 --size_z 18 --seed 1234 --exhaustiveness 64"
```bash
This command will run AutoDock Vina on your receptor and ligand files, placing the result in the output_local directory.

### 7. Expected Output
After running the above command, you should find the output file (rand-1_out.pdbqt) in the output_local directory.

