# Prerequisites: conda env activated


# Clone the repo: https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git

# git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git

cd ../../../../../Open-Omics-Acceleration-Framework
WDIR=`pwd`

cd ${WDIR}/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/

ls


# Pre-req: conda env
source setup_env.sh  dv_env


# compile bwa-mem2
echo "Build bwa-mem2"
cd ${WDIR}/applications/bwa-mem2
make multi
if [ -e "${WDIR}/applications/bwa-mem2/bwa-mem2" ]; then
    echo "bwa-mem2 build successful"
else
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
#make install         #uncomment this for installation


