#!/bin/bash

SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
#echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
#echo "Value of ABS_DIRECTORY: ${ABS_DIRECTORY}"


echo "Downloading and setting up miniconda..."
if [ ! -e "Miniconda3-py39_23.3.1-0-Linux-x86_64.sh" ]
then
    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.3.1-0-Linux-x86_64.sh
fi

bash ./Miniconda3-py39_23.3.1-0-Linux-x86_64.sh -b -p ./miniconda3
echo "Downloading and setting up miniconda...DONE"

echo "Seeting up conda env named with given argument"
miniconda3/bin/conda env create --name bwa -f environment.yml
echo "Seeting up conda env named new_env...DONE"

echo "Activating conda env..."
source /data/nfs_home/mvasimud/container/miniconda3/bin/activate bwa
echo "localhost" > hostfile



## build tools
WDIR=../../

# compile bwa-mem2
echo "Build bwa-mem2"
cd ${WDIR}/applications/bwa-mem2
make -j multi
bwainstall="SUCESS"
if [ -e "${WDIR}/applications/bwa-mem2/bwa-mem2" ]; then
    echo "bwa-mem2 build successful"
else
    bwainstall="FAILED"
    echo "Error!! bwa-mem2 build failed"
fi

#make install   #uncomment this for installation

# compile htslib
cd ${WDIR}/applications/htslib
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
#make install   #uncomment this for installation

# compile bcftools
cd ${WDIR}/applications/bcftools
# The following is optional:
#   autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters
make
#make install   #uncomment this for installation

# compile samtools
cd ${WDIR}/applications/samtools
autoheader
autoconf -Wno-syntax
chmod 775 configure
./configure           # Needed for choosing optional functionality
make

if [ "$?" == "0" ]
then
    echo "Samtools installed successfully"
else
    echo "Samtools installation failed"
fi
echo "bwa compilation is "$bwainstall
#make install         #uncomment this for installation
