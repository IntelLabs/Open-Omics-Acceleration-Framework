## Open-Omics-Autodock-Vina
Open-Omics-Autodock-Vina is a fast, efficient molecular docking software used to predict ligand-protein binding poses and affinities. It features a refined scoring function, parallel execution on multicore CPUs and user-friendly configuration.

## Docker Setup Instructions


### 1. Build the Docker Image 
To build the Docker image with the tag `docker_vina`, use the following commands based on your machine's proxy requirements:
* For machine without a proxy:
```bash
cd build_docker/
docker build -t docker_vina:latest --build-arg FLAVOR=microservice .
cd ..
```
* For machine with a proxy:
```bash
cd build_docker/
docker build --build-arg http_proxy=<http_proxy> --build-arg https_proxy=<https_proxy> --build-arg no_proxy=<no_proxy_ip> FLAVOR=microservice -t docker_vina .
cd ..
```


### 2. Choose and Download Protein Complex Data
Select any protein complex from the available dataset of **140** protein-ligand complexes(https://zenodo.org/records/4031961) which you can download from (https://zenodo.org/records/4031961/files/data.zip?download=1). This guide uses the **5wlo** protein as an example.

1) Run the below commands to make data download script executable, download the complete dataset and extract the data for `5wlo`:

```bash
chmod +x data_download_script.sh
bash data_download_script.sh 5wlo
```
**Note: You can replace 5wlo with any other complex name from the complete dataset available in `data_original/data` directory.**

2) Create an output directory to store results specific to `5wlo`:
```bash
mkdir -p 5wlo_output                                                                                                               
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

## Additional Ways to Run Autodock-vina

OpenOmics supports multiple ways to run Autodock-vina depending on your workflow and scale. Choose the mode that best fits your use case:

### 1. Run as a Microservice
If you want to expose Autodock-vina as a service that can be queried over an API, you can deploy it as a microservice.
Refer to [here](https://github.com/sri480673/Open-Omics-Acceleration-Framework-srilekha/tree/autodock_microservice/applications/AutoDock-Vina/microservice/README.md) for setup instructions and API usage details.

### 2. Run on Cloud Instances (Nextflow)
To run many Autodock-vina tasks over multiple cloud instances, we provide a Nextflow-based option.
This allows you to scale Autodock-vina tasks easily across cloud infrastructure. Refer to [here](https://github.com/sri480673/Open-Omics-Acceleration-Framework-srilekha/tree/autodock_microservice/applications/AutoDock-Vina/nextflow/README.md)

### 3. Run Multiple Processes (Parallel Local Execution)
To run many Autodock-vina tasks on a single machine, you can use the multiprocess tool, which batches your tasks and executes them in parallel using all the available cores.
It automatically determines and configures the optimal level of parallelism.
Read more about the multiprocess tool [here](https://github.com/sri480673/Open-Omics-Acceleration-Framework-srilekha/tree/autodock_microservice/applications/AutoDock-Vina/multiprocess/README.md)


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
