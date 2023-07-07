#ENV=/nfs/work04/ashish/envs/new_env/bin/activate
#source $ENV

SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
#echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
#echo "Value of ABS_DIRECTORY: ${ABS_DIRECTORY}"

export LD_PRELOAD=$LD_PRELOAD:"${ABS_DIRECTORY}/libmimalloc.so.2.0"
#echo $LD_PRELOAD


##############Podman settings ##################
#use can skip this if you don't need this. This need 
#mkdir /tmp/${USER}
#chmod 777 -R ~/.local/
#mkdir -p ~/.local/share
#rm -rf  ~/.local/share/containers
#ln -s /tmp/${USER} ~/.local/share/containers
################################################

cd ${ABS_DIRECTORY}/../../../applications/bwa-mem2
make CXX=icpc multi
#make
#make install   #uncomment this for installation

cd ${ABS_DIRECTORY}/../../../applications/htslib
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
#make install   #uncomment this for installation

cd ${ABS_DIRECTORY}/../../../applications/samtools
autoheader
autoconf -Wno-syntax
chmod 775 configure
./configure           # Needed for choosing optional functionality
make
#make install         #uncomment this for installation
cd ${ABS_DIRECTORY}

bash load_deepvariant.sh
