## Openomics-Autodock-Vina
This guide will help you clone the repository, build the Docker image, set up input/output directories, and run AutoDock Vina inside a Docker container.
## Docker Setup Instructions
### 1. Build the Docker Image                                                                                                   
Build the Docker image by running the following command:
```bash
docker build -t docker_vina .
```
This command creates an image with the tag `docker_vina`. Confirm the image is built by listing Docker images:
```bash
docker images | grep docker_vina
```

### 2. Choose and Download Protein Complex Data
Select any protein complex from the available dataset of **140 protein** complexes which you can download from (https://zenodo.org/records/4031961/files/data.zip?download=1). This guide uses the **5wlo** protein as an example.

1) Run the below commnads to make data download script executable, download the complete dataset and extract the data for `5wlo`:

```bash
chmod +x data_download_script.sh
bash data_download_script.sh 5wlo
```

2) Create an output directory to store results specific to `5wlo`:
```bash
mkdir 5wlo_output                                                                                                               
```

3) Set the environment variables for the `5wlo` protein as follows:
```bash                                                                                                                         
export INPUT_VINA=$PWD/5wlo
export OUTPUT_VINA=$PWD/5wlo_output
```

4) Add the necessary permissions to output folder for Docker to write to it:
```bash
sudo chmod -R a+w $OUTPUT_VINA
```

### 3. Run the Docker Container
Verify that the Docker image was built successfully by listing Docker images:
```bash
docker images | grep docker_vina                                                                                                
```
If the image is listed, run AutoDock Vina with the following command:
```bash                                                                                                                         
docker run -it -v $INPUT_VINA:/input -v $OUTPUT_VINA:/output docker_vina:latest vina --receptor protein.pdbqt --ligand rand-1.pdbqt --out /output/rand-1_out.pdbqt --center_x 16.459 --center_y -19.946 --center_z -5.850 --size_x 18 --size_y 18 --size_z 18 --seed 1234 --exhaustiveness 64
```
This command will process your receptor and ligand files and place the results in the specified output directory.
### 4. Expected Output                                                                                                           
After running the above command, you should find the output file (`rand-1_out.pdbqt`) in the output directory, such as `5wlo_output` for this example.

---
The original README content of AutoDock-Vina follows:

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

## Citations
* [J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. Journal of Chemical Information and Modeling.](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00203)
* [O. Trott and A. J. Olson. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. Journal of computational chemistry, 31(2), 455-461.](https://onlinelibrary.wiley.com/doi/10.1002/jcc.21334)
