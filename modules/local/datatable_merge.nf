process TABLE_MERGE {
    label 'process_medium'

    conda (params.enable_conda ? 'conda-base::r-data.table=1.12.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-data.table:1.12.2' :
        'quay.io/biocontainers/r-data.table:1.12.2' }"

    input:
    path mirtop

    output:
    path "mirna.tsv"   , emit: mirna_tsv
    path "versions.yml", emit: versions

    script:
    """
    collapse_mirtop.r ${mirtop}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
