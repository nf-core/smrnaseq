nextflow.enable.types = true

process MIRNA_IDS {
    tag "${gff ? gff[0].name : 'fasta'}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bcfd41c6d43be719236b5701a481d72b88fb562731c4afb3b56f43bc40d8ef83/data' :
        'community.wave.seqera.io/library/r-base_r-data.table:409a65bae991099c' }"

    input:
    // GFF used for the run: a list with one file, or empty to build the table from the FASTA files only
    gff: List<Path>
    // Formatted FASTA files used for quantification (their names are the ones reported by edgeR)
    tuple(meta2: Map, mature: Path)
    tuple(meta3: Map, hairpin: Path)
    // Raw reference FASTA files given to the pipeline (mature plus optional MirGeneDB star FASTA,
    // and hairpin). Their miRBase headers carry the accessions used when there is no GFF.
    mature_reference: List<Path>
    hairpin_reference: Path
    species: String
    database: String

    output:
    mirna_ids = file('mirna_ids.tsv')

    topic:
    tuple(task.process, 'r-base', eval("Rscript -e 'cat(as.character(getRversion()))'")) >> 'versions'

    script:
    def gff_arg = gff ? "--gff ${gff[0]}" : ''
    """
    mirna_ids.r \\
        ${gff_arg} \\
        --mature ${mature} \\
        --hairpin ${hairpin} \\
        --mature_reference ${mature_reference.join(',')} \\
        --hairpin_reference ${hairpin_reference} \\
        --species '${species}' \\
        --database '${database}' \\
        --output mirna_ids.tsv
    """

    stub:
    """
    touch mirna_ids.tsv
    """
}
