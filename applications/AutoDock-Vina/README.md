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

This guide will help you clone the repository, build the Docker image, set up input/output directories, and run AutoDock Vina (v1.2.2) inside a Docker container.

### 1. Clone the Repository
First, clone the repository and navigate to the AutoDock Vina folder:

```bash
git clone https://github.com/intel-sandbox/TransOmics.OpenOmicsInternal.git
cd TransOmics.OpenOmicsInternal/applications/AutoDock-Vina
```

### 2. Build the Docker Image
Build the Docker image using the following command:
```bash
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy="127.0.0.1,localhost,apt.repo.inel.com" -t docker_vina .
```
Make sure Docker is installed and properly set up on your system.

### 3. Setup Input and Output Directories
Create directories for input and output on your local machine that will be mounted to the Docker container during execution:
```bash
mkdir -p <input_dir>
mkdir -p <output_dir>
```
Add your receptor (`.pdbqt`), ligand (`.pdbqt`) and grid map files to the input directory of your protein.

We have provided a sample protein `5wlo` with all the necessary files (receptor, ligand, grid maps). Create an output directory for storing results specific to `5wlo`:
```bash
mkdir 5wlo_output
```
Set the environment variables for the `5wlo` protein:
```bash
export INPUT_VINA=$PWD/5wlo
export OUTPUT_VINA=$PWD/5wlo_output
```
Add the necessary permissions to the output folder so Docker can write to it:
```bash
sudo chmod -R a+w $OUTPUT_VINA
```
### 4. Run the Docker Container
Check if the Docker image was built successfully by listing Docker images:
```bash
docker images | grep docker_vina
```
If the image is listed, run AutoDock Vina with the following command:

```bash
docker run -it -v $INPUT_VINA:/input -v $OUTPUT_VINA:/output docker_vina vina --receptor protein.pdbqt --ligand rand-1.pdbqt --out /output/rand-1_out.pdbqt --center_x 16.459 --center_y -19.946 --center_z -5.850 --size_x 18 --size_y 18 --size_z 18 --seed 1234 --exhaustiveness 64
```
This command will process your receptor and ligand files and place the results in the specified output directory.

### 5. Expected Output

After running the above command, you should find the output file (`rand-1_out.pdbqt`) in the output directory (e.g `5wlo_output` for the 5wlo example).

## Citations
* [J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. Journal of Chemical Information and Modeling.](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00203)
* [O. Trott and A. J. Olson. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. Journal of computational chemistry, 31(2), 455-461.](https://onlinelibrary.wiley.com/doi/10.1002/jcc.21334)
