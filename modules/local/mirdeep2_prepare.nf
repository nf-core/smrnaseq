// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MIRDEEP2_PIGZ {
    label 'process_low'
    tag "$meta.id"

    conda (params.enable_conda ? 'bioconda::bioconvert=0.4.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconvert:0.4.3--py_0"
    } else {
        container "quay.io/biocontainers/bioconvert:0.4.3--py_0"
    }
    when:
    !params.skip_mirdeep

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq"), emit: reads
    path "*.version.txt" , emit: versions

    script:
    def unzip = reads.toString() - '.gz'
    def software = getSoftwareName(task.process)
    """
    pigz -f -d -p $task.cpus $reads
    echo "2.0.1" > ${software}.version.txt
    """

}
