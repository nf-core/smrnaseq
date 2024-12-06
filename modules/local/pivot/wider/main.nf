process PIVOT_WIDER {
    tag"$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/r-optparse_r-tidyverse_r-vroom:3cbb224fea84a0e1"

    input:
    tuple val(meta), path(csvs)

    output:
    tuple val(meta), path("*joined_samples_mirtop.csv") , emit: csv
    path "versions.yml"                         , emit: versions

    script:
    """
    awk 'NR == 1 || FNR > 1' ${csvs.join(' ')} > final_long_results_temp.csv

    pivot_wider.R \\
        --input final_long_results_temp.csv \\
        --output ${meta.id}_concatenated_temp.csv

    sort -t\$'\t' -k1,1 ${meta.id}_concatenated_temp.csv > joined_samples_mirtop.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyr: \$(Rscript -e "library(limma); cat(as.character(packageVersion('tidyr')))")
        dplyr: \$(Rscript -e "library(limma); cat(as.character(packageVersion('dplyr')))")
        optparse: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('optparse')))")
        vroom: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('vroom')))")
    END_VERSIONS
    """

}
