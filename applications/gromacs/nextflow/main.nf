nextflow.enable.dsl = 2

params.pdb_list      = params.pdb_list      ?: 'pdb_list.txt'
params.pdb_dir       = params.pdb_dir       ?: 'pdbs'

params.custom_script = params.custom_script ?: 'run_commands.sh'

params.mdp_minim     = params.mdp_minim     ?: 'mdtut_minim.mdp'
params.mdp_nvt       = params.mdp_nvt       ?: 'mdtut_nvt.mdp'
params.mdp_npt       = params.mdp_npt       ?: 'mdtut_npt.mdp'
params.mdp_ions      = params.mdp_ions      ?: 'mdtut_ions.mdp'
params.mdp_md        = params.mdp_md        ?: 'mdtut_md.mdp'

params.out_dir       = params.out_dir       ?: 's3://<bucket-name>/gromacs-outputs'


Channel
    .fromPath(params.pdb_list)
    .ifEmpty { error "PDB list file not found: ${params.pdb_list}" }
    .splitText()
    .map { it.trim() }
    .filter { line -> line }     
    .map { pdbName ->
        def fullPath = params.pdb_dir.endsWith('/') ?
                "${params.pdb_dir}${pdbName}" :
                "${params.pdb_dir}/${pdbName}"

        def prefix = pdbName.replaceFirst(/\.pdb$/, '')
        tuple(pdbName, prefix, file(fullPath))
    }
    .set { pdb_jobs }

workflow {

    def custom_script_ch = Channel.value(file(params.custom_script))

    def mdp_minim_ch = Channel.value(file(params.mdp_minim))
    def mdp_nvt_ch   = Channel.value(file(params.mdp_nvt))
    def mdp_npt_ch   = Channel.value(file(params.mdp_npt))
    def mdp_ions_ch  = Channel.value(file(params.mdp_ions))
    def mdp_md_ch    = Channel.value(file(params.mdp_md))

    RUN_GROMACS(
        pdb_jobs,
        custom_script_ch,
        mdp_minim_ch,
        mdp_nvt_ch,
        mdp_npt_ch,
        mdp_ions_ch,
        mdp_md_ch
    )
}

process RUN_GROMACS {

    label 'gromacs'
    tag { prefix }

    publishDir params.out_dir, mode: 'copy', overwrite: true

    input:
    tuple val(pdb_name), val(prefix), path(pdb_file)

    path custom_script
    path mdp_minim
    path mdp_nvt
    path mdp_npt
    path mdp_ions
    path mdp_md

    output:
    path "${prefix}", emit: gromacs_results

    script:
    """
    set -euo pipefail

    echo "PDB name   : ${pdb_name}"
    echo "PDB file   : ${pdb_file}"
    echo "Prefix     : ${prefix}"
    echo "WorkDir    : \$PWD"
    echo "Script     : ${custom_script}"
    echo "MDP minim  : ${mdp_minim}"
    echo "MDP nvt    : ${mdp_nvt}"
    echo "MDP npt    : ${mdp_npt}"
    echo "MDP ions   : ${mdp_ions}"
    echo "MDP md     : ${mdp_md}"

    mkdir -p /input /output

    echo "Staging PDB into /input"
    cp "${pdb_file}" "/input/${pdb_name}"

    echo "Staging MDP files into /input"
    cp "${mdp_minim}" /input/mdtut_minim.mdp
    cp "${mdp_nvt}"   /input/mdtut_nvt.mdp
    cp "${mdp_npt}"   /input/mdtut_npt.mdp
    cp "${mdp_ions}"  /input/mdtut_ions.mdp
    cp "${mdp_md}"    /input/mdtut_md.mdp

    echo "Staging custom script into /input/custom_script.sh"
    cp "${custom_script}" /input/custom_script.sh
    chmod +x /input/custom_script.sh
    export CUSTOM_SCRIPT="custom_script.sh"

    echo "Calling container entrypoint for PDB: ${pdb_name}"
    PWD_BEFORE="\$PWD"
    cd /input

    /entrypoint.sh "${pdb_name}"

    cd "\$PWD_BEFORE"

    echo "Collecting results from /output"
    mkdir -p "${prefix}"

    latest_dir=\$(ls -dt /output/output_* 2>/dev/null | head -n 1 || true)
    if [[ -n "\$latest_dir" ]]; then
        echo "[INFO] Latest output directory: \$latest_dir"
        cp -r "\$latest_dir"/* "${prefix}/" || true
    else
        echo "Warning: No /output/output_* directories found"
    fi

    echo "[INFO] Done for ${prefix}"
    """
}

