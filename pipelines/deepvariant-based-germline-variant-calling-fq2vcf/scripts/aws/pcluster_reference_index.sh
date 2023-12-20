WDIR=`pwd`
num_compute_nodes=1
allocation_time="02:00:00"

if [ -z $2 ]
then 
	echo "Allocating compute nodes by default for 2 hours"
else
	echo "Allocating compute nodes for $2 hours"
	allocation_time=$2
fi

# Allocate compute nodes
salloc --nodes=${num_compute_nodes} --ntasks-per-node=1 --wait-all-nodes=1 --time=${allocation_time} --no-shell &> tmp_salloc && grep "Granted job allocation" tmp_salloc | cut -d" " -f5 &> tmp_jobid

jid=`cat tmp_jobid | head -n 1`

rm tmp_salloc tmp_jobid

srun --jobid=$jid hostname > ../../hostfile

echo "Cluster alloccation done!!"
cat ../../hostfile

for i in `cat ../../hostfile`
do
  echo $i
  ssh $i "bash ${WDIR}/create_reference_index.sh"
done
scancel $jid

