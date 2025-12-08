## Overview: Protein Design Pipeline  
RFdiffusion based protein design is a cutting-edge approach applying diffusion generative models to create novel protein structures with specific shapes and functions. RFdiffusion learns how natural proteins fold and evolve by iteratively “denoising” random structures into realistic 3D backbones that obey biophysical constraints. This allows researchers to design proteins that can bind to targets, form nanostructures, or catalyze reactions, even if such proteins have never existed in nature. By guiding the diffusion process with structural or functional constraints—like binding sites, symmetry, or motifs — RFdiffusion can generate candidate proteins that are both computationally stable and experimentally viable, dramatically accelerating the discovery of new therapeutics, enzymes, and biomaterials.   

This pipeline performs de novo protein design — generating, sequencing, and validating entirely new protein structures that may not exist in nature.
It combines three state-of-the-art AI tools:   
1. **RFdiffusion (RFD)** — Generative diffusion model that designs novel 3D protein backbones   
2. **ProteinMPNN (PMPNN)** — Neural sequence design model that fits amino acids to a given backbone    
3. **AlphaFold2 (AF2)**  — Structure predictor that validates whether the designed sequence actually folds into the desired structure    

The pipeline is computationally demanding. Here we present a highly accelerated version of the pipeline with exactly same features (and command-line options) and ouptut.  
OpenOmics provides the following two variations of the pipeline:
1. `with_pdb_db`: it uses pdb database during the pre-processing stage of AF2
2. `without_pdb_db`: it uses generated pdb from RFD during pre-processing stage of AF2, thus avoiding the use of full blown pdb database    

Option 2 above incurs lightweight AF2 preprocessing as compared to option 1.  

### Build (Docker images)
```bash
git clone https://github.com/intel-sandbox/TransOmics.OpenOmicsInternal.git -b super_pipeline
cd TransOmics.OpenOmicsInternal/pipelines/protein-design
```

```bash
./build_dockers.sh
```
The above command build four docker images corresponding to the three pipeline tools (2 docker images for AF2: Pre-processing and Model Inference)    

### Run (Docker container): with_pdb_db

1. **Download PDB Database for AF2 pre-processing**:
  Follow these instructions from https://github.com/google-deepmind/alphafold/tree/main/scripts to download the AF2 pdb database.
  Download size zipped: ~556 GB, and unzipped size: 2.62 TB.  Store the database at location <af2_pdb_db>

2. **Generate proteins**:    
   a. Motif Scaffolding
    ```bash
    $ protein_design.sh --prep with_pdb_db  \
     --input_pdb <input.pdb> \
     --contig <contigs>  \
     --output_dir <output_folder> \
     --num_designs <num_designs>   \
     --precision <float32/bfloat16>  \
     --db_af2_path <af2_pdb_db>
    ```
    Example:
    ```bash
    $ protein_design.sh --prep with_pdb_db  \
     --input_pdb ./input_pdbs/5TPN.pdb  \
     --contig '[10-40/A163-181/10-40]'  \
     --output_dir ./motif_output  \
     --num_designs 1  \
     --precision bfloat16  \
     --db_af2_path /home/alphafold_database
    ```
    
    b. Binder Design
    ```bash
    $ protein_design.sh --prep with_pdb_db \
     --input_pdb ./input_pdbs/insulin_target.pdb \
     --contig "['A1-150/0 70-100']" \
     --hotspot_res [A59,A83,A91] \
     --output_dir ./binder_output  \
     --num_designs 1  \
     --precision bfloat16  \
     --num_designs 1 \
     --noise_scale_ca 0   \
     --noise_scale_frame 0   \
     --db_af2_path /home/alphafold_database
    ```

    c. Unconditional monomer
    ```bash
    $ protein_design.sh --prep with_pdb_db  \
     --contig "[100-100]"  \
     --output_dir ./uncond_output  \
     --num_designs 1  \
     --precision bfloat16  \
     --db_af2_path /home/alphafold_database
    ```
    
### Run (Docker container): without_pdb_db
**Generate Proteins:**
1. Motif Scaffolding
```bash
$ protein_design.sh --prep without_pdb_db \
 --mode_name motif  \
 --input_pdb ./input_pdbs/5TPN.pdb  \
 --contigs "10-40/A163-181/10-40"  \
 --iterations 50   --num_designs 1   \
 --output_dir ./motifscaffolding_output
```

2. Binder Design
```bash
$ protein_design.sh --prep without_pdb_db \
 --mode_name binder  \
 --input_pdb ./input_pdbs/insulin_target.pdb  \
 --contigs "['/A1-150/0 70-100/']" \
 --hotspot "A59,A83,A91"  --iterations 50  \
 --num_designs 1   \
 --output_dir ./binder_output
```

3. Unconditional-monomer
```bash
$ protein_design.sh --prep without_pdb_db \
 --mode_name unconditional  \
 --contigs "100-100" \
 --hotspot ""  --iterations 50  \
 --num_designs 1   \
 --output_dir ./unconditional_monomer_output
```
