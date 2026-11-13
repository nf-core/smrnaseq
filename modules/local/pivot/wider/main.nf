nextflow.enable.types = true

process PIVOT_WIDER {
    tag"$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f0/f096fe9683b943ce1869bf984e9d240a364bc73ae5647abeae4ff0c3b26ef011/data' :
        'community.wave.seqera.io/library/r-base_r-dplyr_r-optparse_r-tidyr_r-vroom:c983cf5da1eb6e00' }"

    input:
    tuple(meta: Map, csvs: Set<Path>)

    output:
    csv = tuple(meta, file('*joined_samples_mirtop.tsv'))

    topic:
    tuple(task.process, 'r-base', eval("Rscript -e 'cat(as.character(getRversion()))'")) >> 'versions'
    tuple(task.process, 'tidyr', eval("Rscript -e \"cat(as.character(packageVersion('tidyr')))\"")) >> 'versions'
    tuple(task.process, 'dplyr', eval("Rscript -e \"cat(as.character(packageVersion('dplyr')))\"")) >> 'versions'
    tuple(task.process, 'optparse', eval("Rscript -e \"cat(as.character(packageVersion('optparse')))\"")) >> 'versions'
    tuple(task.process, 'vroom', eval("Rscript -e \"cat(as.character(packageVersion('vroom')))\"")) >> 'versions'

    script:
    """
    awk 'NR == 1 || FNR > 1' ${csvs.join(' ')} > final_long_results_temp.csv

    pivot_wider.R \\
        --input final_long_results_temp.csv \\
        --output ${meta.id}_concatenated_temp.csv

    head -n 1 ${meta.id}_concatenated_temp.csv > joined_samples_mirtop.tsv
    tail -n +2 ${meta.id}_concatenated_temp.csv | sort -t\$'\t' -k1,1 >> joined_samples_mirtop.tsv
    """

    stub:
    """
    touch "joined_samples_mirtop.tsv"
    """
}
