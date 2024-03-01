process MIRDEEP2_PIGZ {
    label 'process_low'
    tag "$meta.id"

    // TODO maybe create a mulled container and uncompress within mirdeep2_mapper?
    conda 'bioconda::bioconvert=1.1.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconvert:1.1.1--pyhdfd78af_0' :
        'biocontainers/bioconvert:1.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.{fastq,fq}"), emit: reads
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pigz -f -d -p $task.cpus $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

}
