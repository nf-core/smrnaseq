process BOWTIE_MAP_CONTAMINANTS {
//    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bowtie2=2.4.5' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.5--py39hd2f7db1_2' :
        'quay.io/biocontainers/bowtie2:2.4.5--py36hfca12d5_2' }"

    input:
    tuple val(meta), path(reads) 
    path index
    val contaminant_type

    output:
    tuple val(meta), path("*sam")                               , emit: bam
    tuple val(meta), path('*.filter.unmapped.contaminant.fastq'), emit: unmapped
    path "versions.yml"                                         , emit: versions
    path "filtered.*.stats"                                   , emit: stats

    script:
    def index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    """
    bowtie2 \\
        --threads ${task.cpus} \\
        --very-sensitive-local \\
        -k 1 \\
        -x $index_base \\
        --un ${meta.id}.filter.unmapped.contaminant.fastq \\
        ${reads} \\
        -S ${meta.id}.filter.contaminant.sam > ${meta.id}.contaminant_bowtie.log 2>&1

    # extracting number of reads from bowtie logs
    awk -v type=${contaminant_type} 'BEGIN{tot=0} {if(NR==4 || NR == 5){tot += \$1}} END {print type": "tot }' ${meta.id}.contaminant_bowtie.log | tr -d , > filtered.${meta.id}_${contaminant_type}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' | tr -d '\0')
    END_VERSIONS
    """

}