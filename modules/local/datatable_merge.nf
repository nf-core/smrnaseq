// Import generic module functions
include { saveFiles; getProcessName } from './functions'

params.options = [:]

process TABLE_MERGE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"mirtop", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'conda-base::r-data.table=1.12.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-data.table:1.12.2"
    } else {
        container "quay.io/biocontainers/r-data.table:1.12.2"
    }

    input:
    path mirtop

    output:
    path "mirna.tsv"   , emit: mirna_tsv
    path "versions.yml", emit: versions

    script:
    """
    collapse_mirtop.r ${mirtop}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
