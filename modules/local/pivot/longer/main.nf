process PIVOT_LONGER {
    tag"$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/r-optparse_r-tidyverse_r-vroom:3cbb224fea84a0e1"

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
