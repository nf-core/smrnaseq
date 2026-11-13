nextflow.enable.types = true

process FORMAT_FASTA_MIRNA {
    tag "$fasta"
    label 'process_medium'

    conda 'bioconda::fastx_toolkit=0.0.14'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e4dce7124230c2e3aa2782710246b68b7e5606a1fdafd29fe9d4aaffa2190a9/data' :
        'community.wave.seqera.io/library/fastx_toolkit:0.0.14--2d5a3f28610ed585' }"

    input:
    tuple(meta2: Map, fasta: Path)

    output:
    formatted_fasta = tuple(meta2, file('*_idx.fa'))
    // Not emitted via the versions topic: Nextflow 26.04 never closes a topic fed by a
    // typed process that is included under several aliases (FORMAT_MATURE/FORMAT_HAIRPIN).
    versions        = tuple(task.process, 'fastx_toolkit', '0.0.14')

    script:
    """
    fasta_formatter -w 0 -i $fasta -o ${fasta}_idx.fa
    """

    stub:
    """
    touch ${fasta}_idx.fa
    """
}
