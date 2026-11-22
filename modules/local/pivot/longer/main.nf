nextflow.enable.types = true

process PIVOT_LONGER {
    tag"$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b2b51579f52b51605ebceff6ca4b0edbda2e244e3ee7d826ab85e315f608dd82/data' :
        'community.wave.seqera.io/library/r-base_r-optparse_r-tidyr_r-vroom:ae58a487c93865f0' }"

    input:
    tuple(meta: Map, tsv: Path)

    output:
    csv = tuple(meta, file('*_long.csv'))

    topic:
    tuple(task.process, 'r-base', eval("Rscript -e 'cat(as.character(getRversion()))'")) >> 'versions'
    tuple(task.process, 'tidyr', eval("Rscript -e \"cat(as.character(packageVersion('tidyr')))\"")) >> 'versions'
    tuple(task.process, 'optparse', eval("Rscript -e \"cat(as.character(packageVersion('optparse')))\"")) >> 'versions'

    script:
    """
    pivot_longer.R \\
        --input ${tsv} \\
        --output ${meta.id}_long.csv
    """

    stub:
    """
     touch "${meta.id}_long.csv"
    """

}
