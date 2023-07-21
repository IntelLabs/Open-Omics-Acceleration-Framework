

source config

cd ../..

source miniconda3/bin/activate dv_env

hostname > hostfile


num_cpus=`nproc`

n=`expr ${num_cpus} / 2`

num_ranks=`expr $n / 14`

echo "$num_ranks $n $num_ranks"


sh run_pipeline.sh  ${num_ranks} ${num_ranks} ${REF} ${R1} ${R2} "sudo docker"
