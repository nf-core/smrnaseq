process PIVOT_WIDER {
    tag"$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f0/f096fe9683b943ce1869bf984e9d240a364bc73ae5647abeae4ff0c3b26ef011/data' :
        'community.wave.seqera.io/library/r-base_r-dplyr_r-optparse_r-tidyr_r-vroom:c983cf5da1eb6e00' }"

    input:
    tuple val(meta), path(csvs)

    output:
    tuple val(meta), path("*joined_samples_mirtop.tsv") , emit: csv
    path "versions.yml"                                 , emit: versions

    script:
    """
    awk 'NR == 1 || FNR > 1' ${csvs.join(' ')} > final_long_results_temp.csv

    pivot_wider.R \\
        --input final_long_results_temp.csv \\
        --output ${meta.id}_concatenated_temp.csv

    sort -t\$'\t' -k1,1 ${meta.id}_concatenated_temp.csv > joined_samples_mirtop.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        optparse: \$(Rscript -e "library(optparse); cat(as.character(packageVersion('optparse')))")
        vroom: \$(Rscript -e "library(vroom); cat(as.character(packageVersion('vroom')))")
    END_VERSIONS
    """

    stub:
    """
    touch "joined_samples_mirtop.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        optparse: \$(Rscript -e "library(optparse); cat(as.character(packageVersion('optparse')))")
        vroom: \$(Rscript -e "library(vroom); cat(as.character(packageVersion('vroom')))")
    END_VERSIONS
    """
}
