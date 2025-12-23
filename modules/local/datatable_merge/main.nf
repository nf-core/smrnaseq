process DATATABLE_MERGE {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bcfd41c6d43be719236b5701a481d72b88fb562731c4afb3b56f43bc40d8ef83/data' :
        'community.wave.seqera.io/library/r-base_r-data.table:409a65bae991099c' }"

    input:
    tuple val(meta), path(mirtop)

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

    stub:
    """
    touch "mirna.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
