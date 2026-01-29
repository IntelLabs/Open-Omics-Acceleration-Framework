nextflow.enable.dsl = 2

// ----------------------
// Parameters
// ----------------------
params.tsv       = params.tsv       ?: 'contigs_8.tsv'
params.input_dir = params.input_dir ?: 'pdbs'
params.out_s3    = params.out_s3    ?: 's3://srilekha-new-s3/outputs-motif-demo'

// ----------------------
// Create jobs channel from TSV
// ----------------------
Channel
  .fromPath(params.tsv)
  .ifEmpty { error "TSV not found: ${params.tsv}" }
  .splitText()
  .map { it.trim() }
  .filter { it && !it.startsWith('#') }     // skip empty & comment lines
  .map { line ->
    def parts = line.split(/\s+/)
    if (parts.size() < 2)
      throw new IllegalArgumentException("Bad TSV line: '${line}' (need: PDB CONTIG)")
    def pdb    = parts[0]
    def contig = parts[1]
    def base   = pdb.replaceAll(/\.pdb$/, '')
    def prefix = "design_${base}"
    tuple(pdb, contig, prefix)
  }
  .map { pdbName, contig, prefix ->
    // Build full path to PDB
    def full = params.input_dir.endsWith('/') ?
               "${params.input_dir}${pdbName}" :
               "${params.input_dir}/${pdbName}"
    tuple(file(full), contig, prefix)
  }
  .set { jobs }

// ----------------------
// RFdiffusion process
// ----------------------
process runMotifScaffolding {

  tag { prefix }

  publishDir params.out_s3, mode: 'copy', overwrite: false

  input:
  tuple path(pdb_file), val(contig), val(prefix)

  output:
  path "out/${prefix}",        emit: designed_dir
  path "logs/${prefix}.log",   emit: joblog

  shell:
  '''
  set -euo pipefail

  mkdir -p logs
  LOGFILE="logs/!{prefix}.log"

  {
    echo "=== RFdiffusion motif run ==="
    echo "PDB:    !{pdb_file}"
    echo "CONTIG: !{contig}"
    echo "PREFIX: !{prefix}"
    echo

    # Working directory for RFdiffusion
    mkdir -p run
    cp -f "!{pdb_file}" "run/"
    pdb_local="run/$(basename "!{pdb_file}")"

    # Run RFdiffusion
    python /app/RFdiffusion/scripts/run_inference.py \
        inference.output_prefix=./run/test/!{prefix} \
        inference.input_pdb="$pdb_local" \
        "contigmap.contigs=[\"!{contig}\"]" \
        inference.num_designs=1 \
        inference.precision=bfloat16 \
        inference.deterministic=True

    # Collect outputs in a clean folder
    mkdir -p "out/!{prefix}"
    if [ -d run/test ]; then
        cp -r run/test/* "out/!{prefix}/" || true
    fi

  } &> "$LOGFILE"
  '''
}

// ----------------------
// Workflow entrypoint
// ----------------------
workflow {
  runMotifScaffolding(jobs)
}

