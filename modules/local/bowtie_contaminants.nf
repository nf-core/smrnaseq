process INDEX_CONTAMINANTS {
    label 'process_medium'

    conda (params.enable_conda ? 'bowtie2=2.4.5' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.5--py39hd2f7db1_2' :
        'quay.io/biocontainers/bowtie2:2.4.5--py36hfca12d5_2' }"

    input:
    path fasta

    output:
    path 'fasta_bidx*'  , emit: bt_indices
    path "versions.yml" , emit: versions

    script:
    """
    bowtie2-build ${fasta} fasta_bidx --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

}