WDIR=`pwd`
num_compute_nodes=$1
# Allocate compute nodes
echo "Allocating the compute nodes.."
salloc --nodes=${num_compute_nodes} --ntasks-per-node=1 --wait-all-nodes=1 --time=02:00:00 --no-shell &> tmp_salloc && cat tmp_salloc | cut -d" " -f5 &> tmp_jobid

jid=`cat tmp_jobid | head -n 1`

rm tmp_salloc tmp_jobid

srun --jobid=$jid hostname > ../../hostfile

echo "Cluster alloccation done!!"
cat ../../hostfile

for i in `cat ../../hostfile`
do
  echo $i
  ssh $i "bash ${WDIR}/basic_setup.sh && sudo docker load -i ${WDIR}/deepvariant.tar && sudo docker images && echo \"setup done for $i. Press enter to continue..\" " &
done


