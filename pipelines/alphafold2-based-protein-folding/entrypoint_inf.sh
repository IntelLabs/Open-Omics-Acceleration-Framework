#!/bin/sh

if [ "$1" = "multimer" ]; then
  echo "Running command for multimer"
  mkdir weights && mkdir weights/extracted && python extract_params.py --input /data/params/params_model_1_multimer_v3.npz --output_dir ./weights/extracted/model_1_multimer_v3 \
  && python extract_params.py --input /data/params/params_model_2_multimer_v3.npz --output_dir ./weights/extracted/model_2_multimer_v3 \
  && python extract_params.py --input /data/params/params_model_3_multimer_v3.npz --output_dir ./weights/extracted/model_3_multimer_v3 \
  && python extract_params.py --input /data/params/params_model_4_multimer_v3.npz --output_dir ./weights/extracted/model_4_multimer_v3 \
  && python extract_params.py --input /data/params/params_model_5_multimer_v3.npz --output_dir ./weights/extracted/model_5_multimer_v3 \
  && LD_PRELOAD=/opt/conda/lib/libjemalloc.so:$LD_PRELOAD \
  MALLOC_CONF="oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" \
  python run_multiprocess_infer_multimer.py --root_condaenv=/opt/conda --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --input_dir=/samples --output_dir=/output --model_names="model_1_multimer_v3,model_2_multimer_v3,model_3_multimer_v3,model_4_multimer_v3,model_5_multimer_v3" --num_multimer_predictions_per_model=5
  if [ "$2" = "relax" ]; then
    echo "Running command for relaxation"
    python run_multiprocess_relax.py --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --input_dir=/samples --output_dir=/output --model_names="model_1_multimer_v3,model_2_multimer_v3,model_3_multimer_v3,model_4_multimer_v3,model_5_multimer_v3" --model_preset=multimer --num_multimer_predictions_per_model=5
  fi

elif [ "$1" = "monomer" ]; then
  echo "Running command for monomer"
  mkdir weights && mkdir weights/extracted && python extract_params.py --input /data/params/params_model_1.npz --output_dir ./weights/extracted/model_1 \
  && python extract_params.py --input /data/params/params_model_2.npz --output_dir ./weights/extracted/model_2 \
  && python extract_params.py --input /data/params/params_model_3.npz --output_dir ./weights/extracted/model_3 \
  && python extract_params.py --input /data/params/params_model_4.npz --output_dir ./weights/extracted/model_4 \
  && python extract_params.py --input /data/params/params_model_5.npz --output_dir ./weights/extracted/model_5 \
  && LD_PRELOAD=/opt/conda/lib/libjemalloc.so:$LD_PRELOAD \
  MALLOC_CONF="oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" \
  python run_multiprocess_infer.py --root_condaenv=/opt/conda --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --input_dir=/samples --output_dir=/output --model_names="model_1,model_2,model_3,model_4,model_5"

  if [ "$2" = "relax" ]; then
    echo "Running command for relaxation"
    python run_multiprocess_relax.py --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --input_dir=/samples --output_dir=/output --model_names="model_1,model_2,model_3,model_4,model_5" --model_preset=monomer
  fi

else
  echo "Running command for monomer"
  mkdir weights && mkdir weights/extracted && python extract_params.py --input /data/params/params_model_1.npz --output_dir ./weights/extracted/model_1 \
  && python extract_params.py --input /data/params/params_model_2.npz --output_dir ./weights/extracted/model_2 \
  && python extract_params.py --input /data/params/params_model_3.npz --output_dir ./weights/extracted/model_3 \
  && python extract_params.py --input /data/params/params_model_4.npz --output_dir ./weights/extracted/model_4 \
  && python extract_params.py --input /data/params/params_model_5.npz --output_dir ./weights/extracted/model_5 \
  && LD_PRELOAD=/opt/conda/lib/libjemalloc.so:$LD_PRELOAD \
  MALLOC_CONF="oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" \
  python run_multiprocess_infer.py --root_condaenv=/opt/conda --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --input_dir=/samples --output_dir=/output --model_names="model_1,model_2,model_3,model_4,model_5"
fi
