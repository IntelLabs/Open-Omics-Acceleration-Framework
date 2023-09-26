source config

cd ../..

source miniconda3/bin/activate dv_env

num_nodes=`cat hostfile | wc -l`

first_ip=`head -n 1 hostfile`

#ssh ${first_ip} lscpu > compute_config
lscpu > compute_config


num_cpus_per_node=$(cat compute_config | grep -E '^CPU\(s\)' | awk  '{print $2}')
num_cpus_all_node=`expr ${num_cpus_per_node} \* ${num_nodes}`
threads_per_core=$(cat compute_config | grep -E '^Thread' | awk  '{print $4}')
echo "Total number of CPUs across all nodes: $num_cpus_all_node"


num_physical_cores_all_nodes=`expr ${num_cpus_all_node} / ${threads_per_core}`

num_physical_cores_per_nodes=`expr ${num_cpus_per_node} / ${threads_per_core}`


while [ $num_physical_cores_per_nodes -ge 20 ]
do
   num_physical_cores_per_nodes=`expr $num_physical_cores_per_nodes / 2`
done

num_physical_cores_per_rank=$num_physical_cores_per_nodes

total_num_ranks=`expr ${num_physical_cores_all_nodes} / ${num_physical_cores_per_rank}`

ranks_per_node=`expr ${total_num_ranks} / ${num_nodes}`

sh run_pipeline.sh  ${total_num_ranks} ${ranks_per_node} ${REF} ${R1} ${R2} "sudo docker"

echo "Pipeline finished. Output vcf can be found at: $OUTPUT_DIR/output.vcf.gz"

