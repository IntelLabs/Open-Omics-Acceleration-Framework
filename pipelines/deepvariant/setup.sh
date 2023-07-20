#ENV=/nfs/work04/ashish/envs/new_env/bin/activate
#source $ENV

SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
#echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
#echo "Value of ABS_DIRECTORY: ${ABS_DIRECTORY}"

export LD_PRELOAD=$LD_PRELOAD:"${ABS_DIRECTORY}/libmimalloc.so.2.0"
#echo $LD_PRELOAD

Container=docker

if [ $# -gt 0 ]
then
    Container="$1"
fi

# This will save deepvariant images
cd ${ABS_DIRECTORY}/../../applications/deepvariant
$Container build -t deepvariant .
# docker build --build-arg http_proxy="http://proxy-us.abc.com:123" --build-arg https_proxy="http://proxy-us.abc.com:123" --build-arg no_proxy="127.0.0.1,localhost"  -t deepvariant .


#save image(~7 GB) to tar file if you are using multiple nodes.
cd ${ABS_DIRECTORY}
$Container save -o deepvariant.tar deepvariant:latest


cd ${ABS_DIRECTORY}/../../applications/bwa-mem2
#make CXX=icpc multi
make
#make install   #uncomment this for installation

cd ${ABS_DIRECTORY}/../../applications/htslib
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
#make install   #uncomment this for installation

cd ${ABS_DIRECTORY}/../../applications/bcftools
# The following is optional:
#   autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters
make
#make install   #uncomment this for installation


cd ${ABS_DIRECTORY}/../../applications/samtools
autoheader
autoconf -Wno-syntax
chmod 775 configure
./configure           # Needed for choosing optional functionality
make
#make install         #uncomment this for installation
cd ${ABS_DIRECTORY}

bash load_deepvariant.sh $Container
