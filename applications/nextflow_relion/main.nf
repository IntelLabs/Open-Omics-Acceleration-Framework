nextflow.enable.dsl = 2

/*
 * Params
 */
params.data_root   = params.data_root    ?: 'relion_data'
params.dataset_glob= params.dataset_glob ?: 'data_*'
params.out_dir     = params.out_dir      ?: 's3://srilekha-new-s3/relion-outputs'

params.star_file   = params.star_file    ?: 'particles.star'
params.iter        = params.iter         ?: 25
params.pool        = params.pool         ?: 12
params.pad         = params.pad          ?: 2
params.tau2_fudge  = params.tau2_fudge   ?: 2
params.particle_diameter = params.particle_diameter ?: 200
params.K           = params.K            ?: 50
params.oversampling= params.oversampling ?: 1
params.psi_step    = params.psi_step     ?: 12
params.offset_range= params.offset_range ?: 5
params.offset_step = params.offset_step  ?: 2

params.mpi_ranks   = params.mpi_ranks    ?: 4
params.j_threads   = params.j_threads    ?: 6

params.run_stamp   = params.run_stamp    ?: new Date().format("yyyyMMdd_HHmmss")
params.run_id      = params.run_id       ?: "${params.run_stamp}_${workflow.runName}"
params.publish_dir = params.publish_dir  ?: "${params.out_dir}/${params.run_id}"

/*
 * Input channel: one task per dataset dir (data_1..data_8)
 */
relion_jobs = Channel
  .fromPath("${params.data_root}/${params.dataset_glob}", type: 'dir', checkIfExists: true)
  .map { d -> tuple(d.baseName, d) }   // (dataset_id, dataset_dir)


process RELION_2D_CLASSIFICATION {

  tag { dataset_id }
  label 'relion2d'

  publishDir params.publish_dir, mode: 'copy', overwrite: false

  input:
    tuple val(dataset_id), path(dataset_dir)

  output:
    path "${dataset_id}/Class2D", emit: class2d_out
    path "${dataset_id}/relion_stdout.log", optional: true, emit: relion_log

  script:
  """
  set -euo pipefail

  echo "=========== RELION 2D CLASSIFICATION =========="
  echo "Dataset    : ${dataset_id}"
  echo "Task PWD   : \$(pwd)"
  echo "Staged dir : ${dataset_dir}"
  echo "STAR file  : ${params.star_file}"
  echo "MPI ranks  : ${params.mpi_ranks}"
  echo "Threads -j : ${params.j_threads}"
  echo "RELION     : \$(which relion_refine_mpi || true)"
  echo "MPI        : \$(which mpirun || true)"
  echo "==============================================="

  mkdir -p workdir
  cp -R "${dataset_dir}/." workdir/
  cd workdir

  echo "workdir contents:"
  ls -la

  if [[ ! -f "${params.star_file}" ]]; then
    echo "ERROR: Missing STAR file '${params.star_file}' in dataset ${dataset_id}"
    exit 2
  fi

  mkdir -p Class2D
  OUTDIR="Class2D/run_${dataset_id}"

  set +e
  mpirun -n ${params.mpi_ranks} \$(which relion_refine_mpi) \\
    --o "\${OUTDIR}" \\
    --iter ${params.iter} \\
    --i ${params.star_file} \\
    --dont_combine_weights_via_disc \\
    --preread_images \\
    --pool ${params.pool} \\
    --pad ${params.pad} \\
    --ctf \\
    --tau2_fudge ${params.tau2_fudge} \\
    --particle_diameter ${params.particle_diameter} \\
    --K ${params.K} \\
    --flatten_solvent \\
    --zero_mask \\
    --center_classes \\
    --oversampling ${params.oversampling} \\
    --psi_step ${params.psi_step} \\
    --offset_range ${params.offset_range} \\
    --offset_step ${params.offset_step} \\
    --norm \\
    --scale \\
    --j ${params.j_threads} \\
    --cpu 2>&1 | tee relion_stdout.log
  rc=\${PIPESTATUS[0]:-0}
  set -e

  echo "RELION exit code: \$rc"
  if [[ \$rc -ne 0 ]]; then exit \$rc; fi

  cd ..
  mkdir -p "${dataset_id}/Class2D"
  cp -R workdir/Class2D/* "${dataset_id}/Class2D/" || true
  cp workdir/relion_stdout.log "${dataset_id}/relion_stdout.log" || true
  """
}


workflow {
  RELION_2D_CLASSIFICATION(relion_jobs)
}

