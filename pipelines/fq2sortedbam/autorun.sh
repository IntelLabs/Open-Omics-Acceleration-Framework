source config

mode="pragzip"
if [ "$#" == "1" ]
    then
    mode=$1
else
    echo "<exec> <pragzip/flatmode/fqprocess>"
    exit
fi
echo "mode: "$mode
source miniconda3/bin/activate distbwa

echo "localhost" > hostfile
num_nodes=1
#num_nodes=`cat hostfile | wc -l`
#first_ip=`head -n 1 hostfile`
#echo $first_ip
#ssh ${first_ip} 'lscpu' > compute_config
lscpu > compute_config


num_cpus_per_node=$(cat compute_config | grep -E '^CPU\(s\)' | awk  '{print $2}')
num_socket=$(cat compute_config | grep -E '^Socket'| awk  '{print $2}')
num_numa=$(cat compute_config | grep '^NUMA node(s)' | awk '{print $3}')
num_cpus_all_node=`expr ${num_cpus_per_node} \* ${num_nodes}`
threads_per_core=$(cat compute_config | grep -E '^Thread' | awk  '{print $4}')
echo "Total number of CPUs: $num_cpus_all_node"
echo "Number of sockets: "$num_socket
echo "Number of NUMA domains: "$num_numa

num_physical_cores_all_nodes=`expr ${num_cpus_all_node} / ${threads_per_core}`
num_physical_cores_per_nodes=`expr ${num_cpus_per_node} / ${threads_per_core}`
num_physical_cores_per_socket=`expr ${num_physical_cores_all_nodes} / ${num_socket}`
num_physical_cores_per_numa=`expr ${num_physical_cores_all_nodes} / ${num_numa}`
echo "#############################################"
echo "Num physical cores per nodes: "$num_physical_cores_per_nodes
echo "Num physical cores per socket: "$num_physical_cores_per_socket
echo "Num physical cores per numa: "$num_physical_cores_per_numa

th=`expr ${num_physical_cores_per_numa} / 2`  #${num_physical_cores_per_numa}  ##20
if [ $th -gt ${num_physical_cores_per_numa} ]
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
echo "Cores per node: "$num_physical_cores_per_nodes
echo "Total number of ranks: "${total_num_ranks}
echo "Ranks per node: "${ranks_per_node}
echo "#############################################"
if [ "$R3" == "" ]
then
    R3="NONE"
fi
sh run_bwa.sh  ${total_num_ranks} ${ranks_per_node} ${REF} ${R1} ${R2} ${R3} ${mode}
#echo "Pipeline finished. Output vcf can be found at: $OUTPUT_DIR/output.vcf.gz"
