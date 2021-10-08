// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQCLUSTER_SEQUENCES {
    label 'process_medium'
    tag "$meta.id"

    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::seqcluster=1.2.8-0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqcluster:1.2.8--pyh5e36f6f_0"
    } else {
        container "quay.io/biocontainers/seqcluster:1.2.8--pyh5e36f6f_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("final/*.fastq.gz")    , emit: collapsed
    path "*.version.txt" , emit: versions

    script:
    def software = getSoftwareName(task.process)
    """
    seqcluster collapse -f $reads -m 1 --min_size 15 -o collapsed
    gzip collapsed/*_trimmed.fastq
    mkdir final
    mv collapsed/*.fastq.gz final/.
    echo \$(seqcluster --version 2>&1) | sed 's/^.*seqcluster //' > ${software}.version.txt
    """

}
