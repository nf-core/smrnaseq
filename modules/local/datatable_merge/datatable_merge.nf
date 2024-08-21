process TABLE_MERGE {
    label 'process_medium'

    conda 'conda-forge::r-data.table=1.12.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-data.table:1.12.2' :
        'biocontainers/r-data.table:1.12.2' }"

    input:
    path mirtop

    output:
    path "mirna.tsv"   , emit: mirna_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    collapse_mirtop.r ${mirtop}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
