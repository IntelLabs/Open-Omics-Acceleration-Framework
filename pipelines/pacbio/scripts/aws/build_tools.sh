# Prerequisites: conda env activated


# Clone the repo: https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git

# git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git

cd ../../../../../Open-Omics-Acceleration-Framework
WDIR=`pwd`

cd ${WDIR}/pipelines/pacbio/

# Pre-req: conda env
source setup_env.sh  dv_env


# compile mm2-fast
echo "Build mm2-fast"
cd ${WDIR}/applications/mm2-fast
make
if [ -e "${WDIR}/applications/mm2-fast/minimap2" ]; then
    echo "mm2-fast build successful"
else
    echo "Error!! mm2-fast build failed"
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
#make install         #uncomment this for installation


