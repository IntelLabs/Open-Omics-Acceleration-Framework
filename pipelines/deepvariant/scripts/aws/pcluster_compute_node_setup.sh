WDIR=`pwd`
num_compute_nodes=$1
# Allocate compute nodes
echo "Allocating the compute nodes.."
salloc --nodes=${num_compute_nodes} --ntasks-per-node=1 --time=01:00:00
srun hostname > ../../hostfile

# Assumption --  hostfile is populated
# exiting from here. Need to wait

echo "Cluster alloccation done!!"
cat ../../hostfile

for i in `cat ../../hostfile`
do
  echo $i
  ssh $i "bash ${WDIR}/basic_setup.sh & sudo docker load -i ${WDIR}/deepvariant.tar & sudo docker images" &
done

