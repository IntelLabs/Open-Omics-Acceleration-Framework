## AutoDock Vina: Docking and virtual screening program

**AutoDock Vina** is one of the **fastest** and **most widely used** **open-source** docking engines. It is a turnkey computational docking program that is based on a simple scoring function and rapid gradient-optimization conformational search. It was originally designed and implemented by Dr. Oleg Trott in the Molecular Graphics Lab, and it is now being maintained and develop by the Forli Lab at The Scripps Research Institute.

* AutoDock4.2 and Vina scoring functions
* Support of simultaneous docking of multiple ligands and batch mode for virtual screening
* Support of macrocycle molecules
* Hydrated docking protocol
* Can write and load external AutoDock maps
* Python bindings for Python 3

## Documentation

The installation instructions, documentation and tutorials can be found on [readthedocs.org](https://autodock-vina.readthedocs.io/en/latest/).

## Docker Setup

This repository contains a Docker setup for running AutoDock Vina version 1.2.2. The following instructions will guide you through cloning the repository, building the Docker image, setting up input/output directories, and running AutoDock Vina inside the Docker container.
### 1. Clone the Repository

First, clone the repository and navigate to the AutoDock Vina folder:

```bash
git clone https://github.com/intel-sandbox/TransOmics.OpenOmicsInternal.git
cd TransOmics.OpenOmicsInternal/applications/AutoDock-Vina-1.2.2
```

### 2. Build the Docker Image
Build the Docker image using the following command:
```bash
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy="127.0.0.1,localhost,apt.repo.inel.com" -t docker_vina .
```

### 3. Setup Input and Output Directories
Create the input and output directories on your local machine:
```bash
mkdir input_local
mkdir output_local
```
### 4. Prepare Input Files
Add your receptor (`.pdbqt`), ligand (`.pdbqt`), and dependent map files to the `input_local` directory

### 5. Set Environment Variables
Set environment variables to reference the input and output directories:
```bash
export INPUT_VINA=$PWD/input_local
export OUTPUT_VINA=$PWD/output_local
```
### 6. EXAMPLE: Running with Provided Data
We have provided a folder for the protein `5wlo`, which contains all the necessary grid maps, protein PDBQT, and ligand files. Create an output folder for `5wlo` to store the results, and then export these directories for use with Docker:
```bash
mkdir 5wlo_output
export INPUT_VINA=$PWD/5wlo
export OUTPUT_VINA=$PWD/5wlo_output
```
### 7. Add the necessary permission to output directory
```bash
sudo chmod -R 777 $OUTPUT_VINA
```

### 8. Run the Docker Container
Before running the Docker container, verify that the Docker image named `docker_vina` has been successfully generated:
```bash
docker images | grep docker_vina
```
If the image is listed, proceed to run the Docker container and execute  Autodock Vina:

```bash
docker run -it -v $INPUT_VINA:/input -v $OUTPUT_VINA:/output <docker_image_name> sh -c "cd /input && vina --receptor protein.pdbqt --ligand rand-1.pdbqt --out /output/rand-1_out.pdbqt --center_x 16.459 --center_y -19.946 --center_z -5.850 --size_x 18 --size_y 18 --size_z 18 --seed 1234 --exhaustiveness 64"
```
This command will run AutoDock Vina on your receptor and ligand files, placing the result in the `output_local` directory. For `5wlo` protein output folder will be `5wlo_output`.

### 9. Expected Output

After running the above command, you should find the output file (`rand-1_out.pdbqt`) in the `output_local` directory.

## Citations
* [J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. Journal of Chemical Information and Modeling.](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00203)
* [O. Trott and A. J. Olson. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. Journal of computational chemistry, 31(2), 455-461.](https://onlinelibrary.wiley.com/doi/10.1002/jcc.21334)
