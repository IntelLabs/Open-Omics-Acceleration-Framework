nextflow.enable.dsl = 2

/*
 * PARAMETERS (define FIRST)
 */
params.protein_list = params.protein_list ?: 'protein_list.txt'
params.input_dir    = params.input_dir    ?: 'inputs'
params.out_dir      = params.out_dir      ?: 's3://srilekha-new-s3/autodock-outputs'

params.ligand_file  = params.ligand_file  ?: 'rand-0.pdbqt'
params.nrun         = params.nrun         ?: 100
params.lsmet        = params.lsmet        ?: 'sw'
params.seed         = params.seed         ?: '11,23'
params.nev          = params.nev          ?: 2048000

// Run-level unique folder
params.run_stamp    = params.run_stamp    ?: new Date().format("yyyyMMdd_HHmmss")
params.run_id       = params.run_id       ?: "${params.run_stamp}_${workflow.runName}"
params.publish_dir  = params.publish_dir  ?: "${params.out_dir}/${params.run_id}"


/*
 * Build jobs with a stable index (job_idx)
 * We do: list lines -> indexed -> convert to tuples
 */
Channel
    .fromPath(params.protein_list)
    .ifEmpty { error "Protein list file not found: ${params.protein_list}" }
    .splitText()
    .map { it.trim() }
    .filter { it }
    .toList()
    .map { lines ->
        // lines is a List<String>
        lines.indexed().collect { idx, proteinId ->
            def dirPath = params.input_dir.endsWith('/') ?
                "${params.input_dir}${proteinId}" :
                "${params.input_dir}/${proteinId}"

            tuple(proteinId, idx as int, file(dirPath))
        }
    }
    .flatMap { it }     // emit each tuple one by one
    .set { dock_jobs }


process RUN_AUTODOCK {

    label 'autodock'
    tag { "${protein_id}_job${job_idx}" }

    publishDir params.publish_dir, mode: 'copy', overwrite: true

    input:
    tuple val(protein_id), val(job_idx), path(protein_dir)

    output:
    path "${protein_id}_job${job_idx}", emit: docking_results

    script:
    """
    set -euo pipefail

    echo "================ RUN_AUTODOCK ================"
    echo "Protein ID   : ${protein_id}"
    echo "Job index    : ${job_idx}"
    echo "Protein dir  : ${protein_dir}"
    echo "Ligand       : ${params.ligand_file}"
    echo "nrun         : ${params.nrun}"
    echo "lsmet        : ${params.lsmet}"
    echo "seed         : ${params.seed}"
    echo "nev          : ${params.nev}"
    echo "Publish dir  : ${params.publish_dir}"
    echo "Task PWD     : \$(pwd)"
    echo "Sandbox ls  :"
    ls -la
    echo "Will create workdir at: \$(pwd)/workdir"
    echo "============================================="

    ROOTDIR=\$(pwd)

    mkdir -p workdir
    cp -R "${protein_dir}/." workdir/
    cd workdir

    if [[ ! -f protein.maps.fld ]]; then
        echo "ERROR: protein.maps.fld missing for ${protein_id} (job ${job_idx})"
        ls -R .
        exit 1
    fi

    if [[ ! -f "${params.ligand_file}" ]]; then
        echo "ERROR: Ligand file '${params.ligand_file}' not found for ${protein_id} (job ${job_idx})"
        ls -R .
        exit 1
    fi

    set +e
    autodock_cpu_64wi \\
      --ffile protein.maps.fld \\
      --lfile ${params.ligand_file} \\
      --nrun ${params.nrun} \\
      --lsmet ${params.lsmet} \\
      --seed ${params.seed} \\
      --nev ${params.nev} \\
      --resnam rand-0
    ad_exit=\$?
    set -e

    echo "AutoDock exit code: \$ad_exit"

    if [[ \$ad_exit -ne 0 ]]; then
        if [[ -f rand-0.dlg || -f rand-0.pdbqt ]]; then
            echo "WARNING: autodock exited with code \$ad_exit but output files exist; treating as success."
        else
            echo "ERROR: autodock failed with exit code \$ad_exit and no output files found."
            exit \$ad_exit
        fi
    fi

    cd "\$ROOTDIR"

    OUTDIR="${protein_id}_job${job_idx}"
    mkdir -p "\$OUTDIR"
    cp -R workdir/* "\$OUTDIR/" || true
    """
}

workflow {
    RUN_AUTODOCK(dock_jobs)
}

