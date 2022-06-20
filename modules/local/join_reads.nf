process JOIN_FASTQS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::samtools=1.13' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:40128b496751b037e2bd85f6789e83d4ff8a4837-0' :
        'quay.io/biocontainers/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:40128b496751b037e2bd85f6789e83d4ff8a4837-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(unmapped_meta), path(unmapped_reads)

    output:
    tuple val(meta), path('*_merged.fq.gz'), emit: merged
    script:
    """
    cat ${reads} ${unmapped_reads} > ${meta.id}_merged.fq.gz
    """

}
