def VERSION = '0.0.14'

process FORMAT_FASTA_MIRNA {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::fastx_toolkit=0.0.14-9' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastx_toolkit:0.0.14--he1b5a44_8' :
        'quay.io/biocontainers/fastx_toolkit:0.0.14--he1b5a44_8' }"

    input:
    path fasta

    output:
    path '*_idx.fa'    , emit: formatted_fasta
    path "versions.yml", emit: versions

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
