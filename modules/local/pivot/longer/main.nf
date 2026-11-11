process PIVOT_LONGER {
    tag"$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b2b51579f52b51605ebceff6ca4b0edbda2e244e3ee7d826ab85e315f608dd82/data' :
        'community.wave.seqera.io/library/r-base_r-optparse_r-tidyr_r-vroom:ae58a487c93865f0' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*_long.csv") , emit: csv
    path "versions.yml"                 , emit: versions

    script:
    """
    pivot_longer.R \\
        --input ${tsv} \\
        --output ${meta.id}_long.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        optparse: \$(Rscript -e "library(optparse); cat(as.character(packageVersion('optparse')))")
    END_VERSIONS
    """

    stub:
    """
     touch "${meta.id}_long.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        optparse: \$(Rscript -e "library(optparse); cat(as.character(packageVersion('optparse')))")
    END_VERSIONS
    """

}
