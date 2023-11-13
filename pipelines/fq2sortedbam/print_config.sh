num_nodes=1
lscpu > compute_config


num_cpus_per_node=$(cat compute_config | grep -E '^CPU\(s\)' | awk  '{print $2}')
num_socket=$(cat compute_config | grep -E '^Socket'| awk  '{print $2}')
num_numa=$(cat compute_config | grep '^NUMA node(s)' | awk '{print $3}')
num_cpus_all_node=`expr ${num_cpus_per_node} \* ${num_nodes}`
threads_per_core=$(cat compute_config | grep -E '^Thread' | awk  '{print $4}')
echo "#############################################"
echo "Number of sockets: "$num_socket
echo "Number of NUMA domains: "$num_numa
echo "Number of CPUs: $num_cpus_all_node"

num_physical_cores_all_nodes=`expr ${num_cpus_all_node} / ${threads_per_core}`
num_physical_cores_per_nodes=`expr ${num_cpus_per_node} / ${threads_per_core}`
num_physical_cores_per_socket=`expr ${num_physical_cores_all_nodes} / ${num_socket}`
num_physical_cores_per_numa=`expr ${num_physical_cores_all_nodes} / ${num_numa}`
echo "Number physical cores: "$num_physical_cores_per_nodes
echo "Number physical cores per socket: "$num_physical_cores_per_socket
echo "Number physical cores per numa: "$num_physical_cores_per_numa

th=`expr ${num_physical_cores_per_numa} / 2`  #${num_physical_cores_per_numa}  ##20
if [ $th -le 10 ]
then
    th=${num_physical_cores_per_numa}
fi

while [ $num_physical_cores_per_nodes -gt $th ]
do
    num_physical_cores_per_nodes=`expr $num_physical_cores_per_nodes / 2`
done

num_physical_cores_per_rank=$num_physical_cores_per_nodes
total_num_ranks=`expr ${num_physical_cores_all_nodes} / ${num_physical_cores_per_rank}`

ranks_per_node=`expr ${total_num_ranks} / ${num_nodes}`
echo "Number of MPI ranks: "${total_num_ranks}
echo "Number of cores per MPI rank: "$num_physical_cores_per_nodes
echo "#############################################"
echo "Note: Each MPI rank runs a bwa-mem2 process on its input fastq files produced by fqprocess. Please ensure that the number of files created due to bam_size parameter to fqprocess (in config file) creates number of fastq files equal to ${total_num_ranks}"
echo "#############################################"
