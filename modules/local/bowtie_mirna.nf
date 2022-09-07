process INDEX_MIRNA {
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::bowtie=1.3.0-2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie:1.3.0--py38hcf49a77_2' :
        'quay.io/biocontainers/bowtie:1.3.0--py38hcf49a77_2' }"

    input:
    path fasta

    output:
    path 'fasta_bidx*' , emit: bowtie_indices
    path "versions.yml", emit: versions

    script:
    """
    bowtie-build ${fasta} fasta_bidx --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

}
