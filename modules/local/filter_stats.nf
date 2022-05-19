process FILTER_STATS {
//    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bowtie2=2.4.5' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.5--py39hd2f7db1_2' :
        'quay.io/biocontainers/bowtie2:2.4.5--py36hfca12d5_2' }"

    input:
    path reads
    file stats

    output:
    path "*_mqc.yaml"                                         , emit: stats

    script:
    """
    echo filtered*
    touch test_mqc.yaml
    """

}