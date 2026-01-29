nextflow.enable.dsl = 2

// ----------------------
// Parameters
// ----------------------

params.tsv       = params.tsv       ?:'pmpnn_corrected.tsv'
params.input_dir = params.input_dir ?:'pdbs'
params.out_s3    = params.out_s3    ?:'s3://rahamath-new-s3/outputs-proteinmpnn'
params.num_seq   = params.num_seq   ?:8
params.temp      = params.temp      ?:0.1

// ----------------------
// Create jobs channel
// ----------------------
//println "DEBUG: params.input_dir = ${params.input_dir}"
//Channel
//  .fromPath("${params.input_dir}/**/*.pdb")
// .ifEmpty { error "No PDBs found in ${params.input_dir}" }
//  .map { pdb ->
//      def base = pdb.baseName
//      tuple(pdb, base)
//  }
//  .set { mpnn_jobs }

Channel
  .fromPath(params.tsv)
  .ifEmpty { error "TSV not found: ${params.tsv}" }
  .splitText()
  .map { it.trim() }
  .filter { it && !it.startsWith('#') }
  .map { line ->

    def pdb = line

    def base = pdb.replaceAll(/\.pdb$/, '')
    def prefix = base   // or "mpnn_${base}"

    tuple(
      file("${params.input_dir}/${pdb}"),
      prefix
    )
  }
  .set { mpnn_jobs }

// ----------------------
// ProteinMPNN process
// ----------------------

process runProteinMPNN {

  tag { prefix }

  publishDir params.out_s3, mode: 'copy', overwrite: false

  input:
  tuple path(pdb_file), val(prefix)

  output:
  path "out/${prefix}",        emit: mpnn_out
  path "logs/${prefix}.log",   emit: joblog

  shell:
  '''
  set -euo pipefail

  mkdir -p logs
  LOGFILE="logs/!{prefix}.log"

  {
    echo "=== ProteinMPNN run ==="
    echo "PDB: !{pdb_file}"
    echo "PREFIX: !{prefix}"
    echo

    mkdir -p run
    cp -f "!{pdb_file}" run/
    pdb_local="run/$(basename "!{pdb_file}")"

    python /ProteinMPNN/protein_mpnn_run.py \
        --pdb_path "$pdb_local" \
        --out_folder run/mpnn_out \
        --num_seq_per_target !{params.num_seq} \
        --sampling_temp !{params.temp} \
        --seed 37 \
        --batch_size 1

    mkdir -p "out/!{prefix}"
    cp -r run/mpnn_out/* "out/!{prefix}/"

  } &> "$LOGFILE"
  '''
}

// ----------------------
// Workflow entrypoint
// ----------------------

workflow {
  runProteinMPNN(mpnn_jobs)
}

