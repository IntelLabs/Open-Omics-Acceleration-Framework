#!/bin/bash
set -e

#SCRIPT_PATH="${BASH_SOURCE:-$0}"
#ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
##echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
#ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
##echo "Value of ABS_DIRECTORY: ${ABS_DIRECTORY}"

if [ "$#" -ne "1" ]
then
    echo "pls. provide args: cloud/onprem"
fi

if [ "$1" == "cloud" ]
then
echo "Installing pre-requisite tools.."
bash basic_setup_ubuntu.sh
echo "Done"
fi

echo "Downloading and setting up miniconda..."
[[ ! -e "Miniforge3-24.3.0-0-Linux-x86_64.sh" ]] && wget https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Miniforge3-24.3.0-0-Linux-x86_64.sh
bash ./Miniforge3-24.3.0-0-Linux-x86_64.sh -u -b -p ./miniforge3

echo "Setting up conda env named with given argument"
./miniforge3/bin/conda env create --name fq2bam -f environment.yml
echo "Setting up conda env named new_env...DONE"

echo "Activating conda env..."
source ./miniforge3/bin/activate fq2bam

#echo "Downloading and setting up miniconda..."
#if [ ! -e "Miniconda3-py39_23.3.1-0-Linux-x86_64.sh" ]
#then
#    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.3.1-0-Linux-x86_64.sh
#fi
#
#bash ./Miniconda3-py39_23.3.1-0-Linux-x86_64.sh -b -p ./miniconda3
#echo "Downloading and setting up miniconda...DONE"
#
#echo "Seeting up conda env named with given argument"
#miniconda3/bin/conda env create --name distbwa -f environment.yml
#echo "Seeting up conda env named new_env...DONE"
#
#echo "Activating conda env..."
#source miniconda3/bin/activate distbwa
#echo "localhost" > hostfile

## build tools
WDIR=../../
EXEDIR=`pwd`

# compile bwa-mem2
echo "Build bwa-mem2"
cd ${WDIR}/applications/bwa-mem2
make clean
make -j multi
bwainstall="SUCESS"
if [ -e "${WDIR}/applications/bwa-mem2/bwa-mem2" ]; then
    echo "bwa-mem2 build successful"
else
    bwainstall="FAILED"
    echo "Error!! bwa-mem2 build failed"
fi

cd ${WDIR}/applications/mm2-fast
make clean
make -j
mm2install="SUCESS"
if [ -e "${WDIR}/applications/mm2-fast/minimap2" ]; then
    echo "mm2-fast build successful"
else
    mm2install="FAILED"
    echo "Error!! mm2-fast build failed"
    exit
fi


#make install   #uncomment this for installation

# compile htslib
cd ${WDIR}/applications/htslib
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
#make install   #uncomment this for installation

# compile bcftools
## cd ${WDIR}/applications/bcftools
## # The following is optional:
## #   autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters
## make
## #make install   #uncomment this for installation

# compile samtools
cd ${WDIR}/applications/samtools
autoheader
autoconf -Wno-syntax
chmod 775 configure
./configure           # Needed for choosing optional functionality
make
saminstall="SUCESS"
if [ -e "${WDIR}/applications/samtools/samtools" ]; then
    echo "SAMTools build successful"
else
    saminstall="FAILED"
    echo "Error!! SAMTools build failed"
fi

if [ "$?" == "0" ]
then
    echo "Samtools installed successfully"
else
    echo "Samtools installation failed"
fi
#make install         #uncomment this for installation

cd $EXEDIR

#[[ ! -d warp-tools ]] && git clone --recursive https://github.com/broadinstitute/warp-tools.git -b develop
#cd warp-tools/tools/fastqpreprocessing/
#./fetch_and_make_dep_libs.sh && make
## make -j

#if [ "$?" == "0" ]
#then
#    echo "fqprocess installed successfully"
#else
#    echo "fqprocess installation failed"
#fi

echo "bwa compilation is "$bwainstall
echo "mm2-fast compilation is "$mm2install
echo "samtools compilation is "$saminstall

echo "Compelete installation done."
