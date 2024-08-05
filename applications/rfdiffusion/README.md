## Description

RFdiffusion is an open source method for structure generation, with or without conditional information (a motif, target etc). It can perform a whole range of protein design challenges as we have outlined in the RFdiffusion paper.

## Things Diffusion can do

Motif Scaffolding
Unconditional protein generation
Symmetric unconditional generation (cyclic, dihedral and tetrahedral symmetries currently implemented, more coming!)
Symmetric motif scaffolding
Binder design
Design diversification ("partial diffusion", sampling around a design)

# Table of contents

RFdiffusion
Description
Table of contents
Getting started / installation
Conda Install SE3-Transformer
Get PPI Scaffold Examples
Usage
Running the diffusion script
Basic execution - an unconditional monomer
Motif Scaffolding
The "active site" model holds very small motifs in place
The inpaint_seq flag
A note on diffuser.T
Partial diffusion
Binder Design
Practical Considerations for Binder Design
Fold Conditioning
Generation of Symmetric Oligomers
Using Auxiliary Potentials
Symmetric Motif Scaffolding.
A Note on Model Weights
Things you might want to play with at inference time
Understanding the output files
Docker
Conclusion


# Getting started / installation

Thanks to Sergey Ovchinnikov, RFdiffusion is available as a Google Colab Notebook if you would like to run it there!

We strongly recommend reading this README carefully before getting started with RFdiffusion, and working through some of the examples in the Colab Notebook.

If you want to set up RFdiffusion locally, follow the steps below:

To get started using RFdiffusion, clone the repo:



```bash
source setup_conda.sh
```

You'll then need to download the model weights into the RFDiffusion directory.


```bash
#Usage
mkdir models 
bash download_models.sh /path/to/download/directory/models
```

# Run a docker container
```bash
mkdir output
export MODEL_DIR=<path-to-database-directory>
export OUTPUT_DIR=<path-to-output-directory>

docker run -v $MODEL_DIR:/app/RFdiffusion/models 
           -v $OUTPUT_FOLDER:/app/RFdiffusion/example_outputs rfdiffusion:latest 
           ./scripts/run_inference.py inference.output_prefix=example_outputs/design_motifscaffolding                                                                            inference.input_pdb=examples/input_pdbs/5TPN.pdb 
           'contigmap.contigs=[10-40/A163-181/10-40]' inference.num_designs=1
```


# Run a RFdiffusion Standalone

```bash
source setup_conda.sh

cp -r models RFdiffusion/
```

# Optional you can Install manually
To get started using RFdiffusion, clone the repo:

```bash
git clone https://github.com/RosettaCommons/RFdiffusion.git
```

```bash
cd RFdiffusion
mkdir models && cd models
wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt

Optional:
wget http://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt

# original structure prediction weights
wget http://files.ipd.uw.edu/pub/RFdiffusion/1befcb9b28e2f778f53d47f18b7597fa/RF_structure_prediction_weights.pt
```

# Conda Install SE3-Transformer

Ensure that you have either Anaconda or Miniconda installed.

You also need to install NVIDIA's implementation of SE(3)-Transformers Here is how to install the NVIDIA SE(3)-Transformer code:

```bash
conda env create -f env/SE3nv.yml

conda activate SE3nv
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../.. # change into the root directory of the repository
pip install -e . # install the rfdiffusion module from the root of the repository
```

Anytime you run diffusion you should be sure to activate this conda environment by running the following command:

```bash
conda activate SE3nv
```

# Get PPI Scaffold Examples

To run the scaffolded protein binder design (PPI) examples, we have provided some example scaffold files (examples/ppi_scaffolds_subset.tar.gz). You'll need to untar this:

```bash
tar -xvf examples/ppi_scaffolds_subset.tar.gz -C examples/
```
We will explain what these files are and how to use them in the Fold Conditioning section.

# Running the diffusion script


