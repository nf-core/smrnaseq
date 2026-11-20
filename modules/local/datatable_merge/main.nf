nextflow.enable.types = true

process DATATABLE_MERGE {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bcfd41c6d43be719236b5701a481d72b88fb562731c4afb3b56f43bc40d8ef83/data' :
        'community.wave.seqera.io/library/r-base_r-data.table:409a65bae991099c' }"

    input:
    tuple(meta: Map, mirtop: Path)

    output:
    mirna_tsv = file('mirna.tsv')

    topic:
    tuple(task.process, 'r-base', eval("Rscript -e 'cat(as.character(getRversion()))'")) >> 'versions'

    script:
    """
    collapse_mirtop.r ${mirtop}
    """

    stub:
    """
    touch "mirna.tsv"
    """
}
