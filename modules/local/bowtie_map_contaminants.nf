process BOWTIE_MAP_CONTAMINANTS {
    label 'process_medium'

    conda 'bowtie2=2.4.5'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.5--py39hd2f7db1_2' :
        'biocontainers/bowtie2:2.4.5--py39hd2f7db1_2' }"

    input:
    tuple val(meta), path(reads)
    path index
    val contaminant_type

    output:
    tuple val(meta), path("*sam")                               , emit: bam
    tuple val(meta), path('*.filter.unmapped.contaminant.fastq'), emit: unmapped
    path "versions.yml"                                         , emit: versions
    path "filtered.*.stats"                                     , emit: stats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""

    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    bowtie2 \\
        -x \$INDEX \\
        -U ${reads} \\
        --threads ${task.cpus} \\
        --un ${meta.id}.${contaminant_type}.filter.unmapped.contaminant.fastq \\
        --very-sensitive-local \\
        -k 1 \\
        -S ${meta.id}.filter.contaminant.sam \\
        ${args} \\
        > ${meta.id}.contaminant_bowtie.log 2>&1

    # extracting number of reads from bowtie logs
    awk -v type=${contaminant_type} 'BEGIN{tot=0} {if(NR==4 || NR == 5){tot += \$1}} END {print "\\""type"\\": "tot }' ${meta.id}.contaminant_bowtie.log | tr -d , > filtered.${meta.id}_${contaminant_type}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' | tr -d '\0')
    END_VERSIONS
    """

}
