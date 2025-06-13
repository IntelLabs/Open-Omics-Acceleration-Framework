## 🔍 Running Inference with Boltz Docker

Follow the steps below to run inference using the Boltz Docker container:

---

### 🐳 1. Build the Docker Image

From the root of the project directory, build the Docker image:

```bash
docker build -t boltz1 .
```

---

### 📁 2. Create and Set Output Directory Permissions

Create an output folder and give it proper write permissions:

```bash
mkdir -p <output_folder_location> <model_folder_location>
chmod a+w <output_folder_location> <model_folder_location>
          
export OUTPUT=$PWD/<output_folder_location>
export MODELS=$PWD/<model_folder_location>
export INPUT=$PWD/<input_folder_location>
```

> ⚠️ Docker needs write permissions in the `<output_folder_location>` and `<model_folder_location>`  folder. `<input_folder_location>` is the folder contaning the input `.yaml` or `.fasta` file

Example

```bash
mkdir -p ./output ./model
chmod a+w ./output ./model
          
export OUTPUT=$PWD/output
export MODELS=$PWD/model
export INPUT=$PWD/examples/
```

---

### 🚀 3. Run Inference

In order to do inferencing few things needs to be done
Mount the volumes for input folder and output folder. Pass the mounted volumes to boltz as arguments. So the docker run command looks like

```bash
docker run -it \
  --shm-size=100g \
  -v $INPUT:/app/boltz/input \
  -v $MODELS:/home/boltz-service/.boltz/ \
  -v $OUTPUT:/app/boltz/output \
  boltz1
```

> 📝 The `--shm-size=100g` flag avoids shared memory issues during data loading with PyTorch.

---

### ✅ Output

Results will be written to the <output_folder_location> folder.

Boltz currently accepts three input formats:

1. Fasta file, for most use cases

2. A comprehensive YAML schema, for more complex use cases

3. A directory containing files of the above formats, for batched processing

## For more information checkout [boltz](https://github.com/jwohlwend/boltz)

## License

Our model and code are released under MIT License, and can be freely used for both academic and commercial purposes.


## Cite

If you use this code or the models in your research, please cite the following paper:

```bibtex
@article{wohlwend2024boltz1,
  author = {Wohlwend, Jeremy and Corso, Gabriele and Passaro, Saro and Reveiz, Mateo and Leidal, Ken and Swiderski, Wojtek and Portnoi, Tally and Chinn, Itamar and Silterra, Jacob and Jaakkola, Tommi and Barzilay, Regina},
  title = {Boltz-1: Democratizing Biomolecular Interaction Modeling},
  year = {2024},
  doi = {10.1101/2024.11.19.624167},
  journal = {bioRxiv}
}
```

In addition if you use the automatic MSA generation, please cite:

```bibtex
@article{mirdita2022colabfold,
  title={ColabFold: making protein folding accessible to all},
  author={Mirdita, Milot and Sch{\"u}tze, Konstantin and Moriwaki, Yoshitaka and Heo, Lim and Ovchinnikov, Sergey and Steinegger, Martin},
  journal={Nature methods},
  year={2022},
}
```
