process INDEX_MIRNA {
    label 'process_medium'

    conda 'bioconda::bowtie=1.3.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie:1.3.1--py310h7b97f60_6' :
        'biocontainers/bowtie:1.3.1--py310h7b97f60_6' }"

    input:
    tuple val(meta2), path(fasta)

    output:
    path 'fasta_bidx*' , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bowtie-build ${fasta} fasta_bidx --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

}
