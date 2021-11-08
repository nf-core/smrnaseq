// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]

process MIRDEEP2_PIGZ {
    label 'process_low'
    tag "$meta.id"

    // TODO maybe create a mulled container and uncompress within mirdeep2_mapper?
    conda (params.enable_conda ? 'bioconda::bioconvert=0.4.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconvert:0.4.3--py_0"
    } else {
        container "quay.io/biocontainers/bioconvert:0.4.3--py_0"
    }
    when:
    !params.skip_mirdeep  // TODO ? I think it would be better to have this logic outside the module

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq"), emit: reads
    path "versions.yml"          , emit: versions

    script:
    def unzip = reads.toString() - '.gz'
    """
    pigz -f -d -p $task.cpus $reads

    ${getProcessName(task.process)}:
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

}
