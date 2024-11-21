#!/bin/sh

if [ "$1" = "multimer" ]; then
  echo "Running command for multimer"
  mkdir weights && mkdir weights/extracted && python extract_params.py --input /data/params/params_model_1_multimer_v3.npz --output_dir ./weights/extracted/model_1_multimer_v3 \
    && python extract_params.py --input /data/params/params_model_2_multimer_v3.npz --output_dir ./weights/extracted/model_2_multimer_v3 \
    && python extract_params.py --input /data/params/params_model_3_multimer_v3.npz --output_dir ./weights/extracted/model_3_multimer_v3 \
    && python extract_params.py --input /data/params/params_model_4_multimer_v3.npz --output_dir ./weights/extracted/model_4_multimer_v3 \
    && python extract_params.py --input /data/params/params_model_5_multimer_v3.npz --output_dir ./weights/extracted/model_5_multimer_v3 \
    && python run_multiprocess_pre_multimer.py --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --data_dir=/data --input_dir=/samples --output_dir=/output
else
  echo "Running command for monomer"
  mkdir weights && mkdir weights/extracted && python extract_params.py --input /data/params/params_model_1.npz --output_dir ./weights/extracted/model_1 \
  && python run_multiprocess_pre.py --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --data_dir=/data --input_dir=/samples --output_dir=/output --model_name=model_1
fi