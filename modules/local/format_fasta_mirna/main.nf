def VERSION = '0.0.14'

process FORMAT_FASTA_MIRNA {
    tag "$fasta"
    label 'process_medium'

    conda 'bioconda::fastx_toolkit=0.0.14'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e4dce7124230c2e3aa2782710246b68b7e5606a1fdafd29fe9d4aaffa2190a9/data' :
        'community.wave.seqera.io/library/fastx_toolkit:0.0.14--2d5a3f28610ed585' }"

    input:
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta2), path('*_idx.fa') , emit: formatted_fasta
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    fasta_formatter -w 0 -i $fasta -o ${fasta}_idx.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastx_toolkit:  \$(echo "$VERSION")
    END_VERSIONS
    """
}
