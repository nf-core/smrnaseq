nextflow.enable.types = true

process MIRNA_IDS {
    tag "$gff"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bcfd41c6d43be719236b5701a481d72b88fb562731c4afb3b56f43bc40d8ef83/data' :
        'community.wave.seqera.io/library/r-base_r-data.table:409a65bae991099c' }"

    input:
    tuple(meta: Map, gff: Path)
    tuple(meta2: Map, mature: Path)
    tuple(meta3: Map, hairpin: Path)
    species: String
    database: String

    output:
    mirna_ids = file('mirna_ids.tsv')

    topic:
    tuple(task.process, 'r-base', eval("Rscript -e 'cat(as.character(getRversion()))'")) >> 'versions'

    script:
    """
    mirna_ids.r \\
        --gff ${gff} \\
        --mature ${mature} \\
        --hairpin ${hairpin} \\
        --species '${species}' \\
        --database '${database}' \\
        --output mirna_ids.tsv
    """

    stub:
    """
    touch mirna_ids.tsv
    """
}
