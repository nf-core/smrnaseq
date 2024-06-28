process BOWTIE_MAP_SEQ {
    tag "$meta.id"
    label 'process_medium'

    conda 'bowtie=1.3.0 bioconda::samtools=1.13'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:40128b496751b037e2bd85f6789e83d4ff8a4837-0' :
        'biocontainers/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:40128b496751b037e2bd85f6789e83d4ff8a4837-0' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*bam")           , emit: bam
    tuple val(meta), path('unmapped/*fq.gz'), emit: unmapped
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    INDEX=`find -L ./ -name "*.3.ebwt" | sed 's/.3.ebwt//'`
    bowtie \\
        -x \$INDEX \\
        -q <(zcat $reads) \\
        -p ${task.cpus} \\
        -t \\
        -k 50 \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        --un ${meta.id}_unmapped.fq -S > ${meta.id}.sam

    samtools view -bS ${meta.id}.sam > ${meta.id}.bam

    if [ ! -f  "${meta.id}_unmapped.fq" ]
    then
        touch ${meta.id}_unmapped.fq
    fi
    gzip ${meta.id}_unmapped.fq
    mkdir unmapped
    mv  ${meta.id}_unmapped.fq.gz  unmapped/.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}
