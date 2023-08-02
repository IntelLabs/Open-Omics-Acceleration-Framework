source config

cd ../..

source miniconda3/bin/activate dv_env




num_nodes=`cat hostfile | wc -l`
num_cpus=`nproc`  
num_cpus=`expr ${num_cpus} \* ${num_nodes}`


echo $num_cpus


n=`expr ${num_cpus} / 2`

num_ranks=`expr $n / 14`



echo "$num_ranks $n $num_ranks"

echo "${num_ranks} `expr ${num_ranks} / ${num_nodes}`"


sh run_pipeline.sh  ${num_ranks} `expr $num_ranks} / ${num_nodes}` ${REF} ${R1} ${R2} "sudo docker"
