// Import generic module functions
include { saveFiles; initOptions; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQCLUSTER_SEQUENCES {
    label 'process_medium'
    tag "$meta.id"

    conda (params.enable_conda ? 'bioconda::seqcluster=1.2.8-0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqcluster:1.2.8--pyh5e36f6f_0"
    } else {
        container "quay.io/biocontainers/seqcluster:1.2.8--pyh5e36f6f_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("final/*.fastq.gz"), emit: collapsed
    path "versions.yml"                      , emit: versions

    script:
    """
    seqcluster collapse -f $reads -m 1 --min_size 15 -o collapsed
    gzip collapsed/*_trimmed.fastq
    mkdir final
    mv collapsed/*.fastq.gz final/.

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(seqcluster --version 2>&1) | sed 's/^.*seqcluster //')
    END_VERSIONS
    """

}
