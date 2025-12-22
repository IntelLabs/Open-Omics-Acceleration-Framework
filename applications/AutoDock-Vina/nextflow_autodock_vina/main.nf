nextflow.enable.dsl = 2

params.list     = params.list     ?: 'complexes_16.txt'
params.data_dir = params.data_dir ?: 'data'
params.out_s3   = params.out_s3   ?: 's3://srilekha-new-s3/outputs-vina-test'

Channel
  .fromPath(params.list)
  .ifEmpty { error "List not found: ${params.list}" }
  .splitText()
  .map { it.trim() }
  .filter { it && !it.startsWith('#') }
  .map { complex ->
    def dir = file("${params.data_dir}/${complex}")
    tuple(complex, dir)
  }
  .set { jobs }


process runVinaRand0 {

  tag { complex }

  publishDir params.out_s3, mode: 'copy', overwrite: true

  input:
  tuple val(complex), path(complex_dir)

  output:
  path "${complex}", emit: vina_results

  script:
  """
  set -euo pipefail

  echo "================ runVinaRand0 ================"
  echo "Complex     : ${complex}"
  echo "Complex dir : ${complex_dir}"
  echo "CPUS        : ${task.cpus}"
  echo "=============================================="

  ROOTDIR=\$(pwd)

  # Make a local working copy (robust on AWS Batch staging)
  mkdir -p workdir
  cp -R "${complex_dir}/." workdir/
  cd workdir

  echo "Listing workdir contents:"
  ls -lah

  # Sanity checks
  if [[ ! -f protein.pdbqt ]]; then
    echo "ERROR: protein.pdbqt missing"
    exit 1
  fi
  if [[ ! -f rand-0.pdbqt ]]; then
    echo "ERROR: rand-0.pdbqt missing"
    exit 1
  fi
  if [[ ! -f vina_rand-0_config.txt ]]; then
    echo "ERROR: vina_rand-0_config.txt missing"
    exit 1
  fi

  echo "Running Vina..."
  vina --config vina_rand-0_config.txt --cpu ${task.cpus}

  cd "\$ROOTDIR"

  echo "Collecting outputs..."
  mkdir -p "${complex}"
  cp -R workdir/* "${complex}/" || true

  echo "Done for ${complex}"
  """
}

workflow {
  runVinaRand0(jobs)
}

