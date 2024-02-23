process SEQCLUSTER_SEQUENCES {
    label 'process_medium'
    tag "$meta.id"

    conda 'bioconda::seqcluster=1.2.9'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqcluster:1.2.9--pyh5e36f6f_0' :
        'biocontainers/seqcluster:1.2.9--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("final/*.fastq.gz"), emit: collapsed
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    seqcluster collapse -f $reads -m 1 --min_size 15 -o collapsed
    gzip collapsed/*_trimmed.fastq
    mkdir final
    mv collapsed/*.fastq.gz final/.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqcluster: \$(echo \$(seqcluster --version 2>&1) | sed 's/^.*seqcluster //')
    END_VERSIONS
    """

}
