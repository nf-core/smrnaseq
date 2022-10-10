process MIRDEEP2_PIGZ {
    label 'process_low'
    tag "$meta.id"

    // TODO maybe create a mulled container and uncompress within mirdeep2_mapper?
    conda (params.enable_conda ? 'bioconda::bioconvert=0.4.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconvert:0.4.3--py_0' :
        'quay.io/biocontainers/bioconvert:0.4.3--py_0' }"

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
    ${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

}
