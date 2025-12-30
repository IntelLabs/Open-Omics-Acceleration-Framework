# Open-Omics-Gromacs
Open-Omics-Gromacs is a dockerized GROMACS setup for high-performance molecular dynamics simulations. It uses multi-threading and Intel MKL FFT for enhanced computational efficiency on modern Intel CPUs.

# Docker Setup Instructions

## 1. Build the Docker Image

### If Your Machine Uses a Proxy
In some cases, systems within a corporate or institutional network require a proxy to connect to the internet. To enable internet connectivity during the Docker build process, specify the proxy environment variables (`http_proxy`, `https_proxy`, and `no_proxy`) during the docker build process.
```zsh
docker build --build-arg http_proxy=<http_proxy> --build-arg https_proxy=<https_proxy> --build-arg no_proxy=<no_proxy_ip> -t grms_dock .
```
**Example Proxy Variables**:
```zsh
http_proxy=http://proxy.example.com:port
https_proxy=http://proxy.example.com:port
no_proxy=localhost,127.0.0.1,.example.com
```

### If Your Machine Doesnot Use Proxy:
```zsh
docker build -t grms_dock .
```

After successfully building the Docker image, verify its presence:
```zsh
docker images |  grep grms_dock
```
## 2. Prepare Input and Output Directories
After Building docker image, set up the necessary **input and output directories**:
### Steps:
1) **If you already have a protein PDB file, you can use it for GROMACS simulation. Otherwise, you can download a sample protein PDB file using the following command:**
```zsh
wget https://files.rcsb.org/download/1AKI.pdb
```

2) **Create a directory for storing simulation output:**
```zsh
mkdir -p grms_output
```

3) **Set environment variables for easy access to input and output directories:**
```zsh
export INPUT_GMX_CPU=$PWD/grms_input
export OUTPUT_GMX_CPU=$PWD/grms_output
```
4) **Prepare Input Directory (`grms_input`) with Necessary Files:**<br>
Copy the following files to the `grms_input` directory:
- `.mdp` files: These are already provided to define simulation parameters.
- `<protein_name>.pdb`: Replace `protein_name` with the name of the downloaded PDB file or your own PDB file.
- `run_commands.sh`: The script used to run the complete workflow of the simulation process.
5) **Set appropriate permissions:**
```zsh
chown -R 1001:1001 $INPUT_GMX_CPU
chown -R 1001:1001 $OUTPUT_GMX_CPU
chmod -R a+w $OUTPUT_GMX_CPU
chmod +x grms_input/run_commands.sh
```

## 3. Run the Docker Container
										   
### Run Complete Workflow with a Selected PDB File
```zsh
docker run -v $INPUT_GMX_CPU:/input -v $OUTPUT_GMX_CPU:/output -it grms_dock <protein_name>.pdb
```
**Explanation**:
* `-v $INPUT_GMX_CPU:/input` mounts the local input directory to `/input` inside the container.
* `-v $OUTPUT_GMX_CPU:/output` mounts the local output directory to `/output` inside the container.
* `<protein_name>.pdb` is the selected PDB file for simulation.

**Note**: If you encounter errors related to missing atoms, non-standard residues, incomplete terminal groups, residue numbering, or naming mismatches during the pdb2gmx stage in GROMACS, PDBFixer(https://github.com/openmm/pdbfixer) is a helpful tool that can clean and fix PDB files before passing them to GROMACS. Additionally, pay close attention to warnings and notes from GROMACS, as they often provide valuable insights and indicate necessary modifications for a successful run.

### Run an Individual GROMACS Command
GROMACS command line utilities can be called directly using our docker container in the following way:

```zsh
docker run -v $INPUT_GMX_CPU:/input -v $OUTPUT_GMX_CPU:/output -it grms_dock gmx <command> 
```
**EXAMPLE**
```zsh
docker run -v $INPUT_GMX_CPU:/input -v $OUTPUT_GMX_CPU:/output -it grms_dock gmx pdb2gmx 
```
For detailed command usage, visit (https://manual.gromacs.org/documentation/current/onlinehelp/).

## 4. Access Simulation Results                            
Once the simulation completes, access the results:
```zsh
ls $OUTPUT_GMX_CPU
```

This directory contains all output data for the **md01** stage of the complete workflow simulation. For individual command runs, everything will be stored in the `grms_input` directory.

---

**The Original README for Gromacs starts here:**


 Welcome to the official version of GROMACS!

If you are familiar with Unix, it should be fairly trivial to compile and
install GROMACS. GROMACS uses only the CMake build system, and our
installation guide can be found at
http://manual.gromacs.org/documentation/current/install-guide/index.html

Visit http://forums.gromacs.org/ for discussions and advice.
Report bugs at https://gitlab.com/gromacs/gromacs/-/issues
Of course we will do our utmost to help you with any problems, but PLEASE
READ THE INSTALLATION INSTRUCTIONS BEFORE CONTACTING US!

There are also several other online resources available from the homepage,
and special information for developers.

If you use GROMACS for research (or other purposes) make sure to reference
the correct version of the code. The CITATION.cff file next to this README
contains all metadata useful for citation.

If you are a developer, or change the source for any other reason, check
out http://www.gromacs.org/development.


GROMACS is free software, distributed under the GNU Lesser General
Public License, version 2.1 However, scientific software is a little
special compared to most other programs. Both you, we, and all other
GROMACS users depend on the quality of the code, and when we find bugs
(every piece of software has them) it is crucial that we can correct
it and say that it was fixed in version X of the file or package
release. For the same reason, it is important that you can reproduce
other people's result from a certain GROMACS version.

The easiest way to avoid this kind of problems is to get your modifications
included in the main distribution. We'll be happy to consider any decent
code. If it's a separate program it can probably be included in the contrib
directory straight away (not supported by us), but for major changes in the
main code we appreciate if you first test that it works with (and without)
MPI, threads, double precision, etc.

If you still want to distribute a modified version or use part of GROMACS
in your own program, remember that the entire project must be licensed
according to the requirements of the LGPL v2.1 license under which you
received this copy of GROMACS. We request that it must clearly be labeled as
derived work. It should not use the name "official GROMACS", and make
sure support questions are directed to you instead of the GROMACS developers.
Sorry for the hard wording, but it is meant to protect YOUR research results!



The development of GROMACS is mainly funded by academic research grants.
To help us fund development, we humbly ask that you cite the GROMACS papers:

* GROMACS: A message-passing parallel molecular dynamics implementation
  H.J.C. Berendsen, D. van der Spoel and R. van Drunen
  Comp. Phys. Comm. 91, 43-56 (1995)
  DOI: https://doi.org/10.1016/0010-4655(95)00042-E

* GROMACS 4: Algorithms for highly efficient, load-balanced, and scalable
  molecular simulation
  B. Hess and C. Kutzner and D. van der Spoel and E. Lindahl
  J. Chem. Theory Comput. 4 (2008) pp. 435-447
  DOI: https://doi.org/10.1021/ct700301q

* GROMACS 4.5: a high-throughput and highly parallel open source
  molecular simulation toolkit
  Sander Pronk, Szilárd Páll, Roland Schulz, Per Larsson, Pär Bjelkmar,
  Rossen Apostolov, Michael R. Shirts, Jeremy C. Smith, Peter M. Kasson,
  David van der Spoel, Berk Hess, Erik Lindahl.
  Bioinformatics 29 (2013) pp. 845-54
  DOI: https://doi.org/10.1093/bioinformatics/btt055

* Tackling Exascale Software Challenges in Molecular Dynamics Simulations
  with GROMACS
  Szilárd Páll, Mark J. Abraham, Carsten Kutzner, Berk Hess, Erik Lindahl
  In S. Markidis & E. Laure (Eds.), Solving Software Challenges for Exascale,
  Lecture Notes for Computer Science, 8759 (2015) pp. 3–27
  DOI: https://doi.org/10.1007/978-3-319-15976-8_1

* GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers
  M. J. Abraham, T. Murtola, R. Schulz, S. Páll, J. C. Smith, B. Hess, E. Lindahl,
  SoftwareX, 1, (2015), 19-25
  DOI: https://doi.org/10.1016/j.softx.2015.06.001

There are a lot of cool features we'd like to include in future versions,
but our resources are limited. All kinds of donations are welcome, both in
form of code, hardware and funding! Industrial users who choose to pay
for a license pro bono (it is still LGPL and can be redistributed freely) or
contribute in other ways are listed as GROMACS supporters on our webpages.
Don't hesitate to contact us if you are interested.
