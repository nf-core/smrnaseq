process FILTER_STATS {
    label 'process_medium'
    tag "$meta.id"

    conda 'bowtie2=2.4.5'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.5--py39hd2f7db1_2' :
        'biocontainers/bowtie2:2.4.5--py39hd2f7db1_2' }"

    input:
    tuple val(meta), path(reads), path (stats_files)

    output:
    path "*_mqc.yaml"                           , emit: stats
    tuple val(meta), path('*.filtered.fastq.gz'), emit: reads, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    if [[ ${reads} == *.gz ]]; then
        readnumber=\$(zcat ${reads} | wc -l | awk '{ print \$1/4 }')
    else
        readnumber=\$(wc -l ${reads} | awk '{ print \$1/4 }')
    fi

    cat ./*${meta.id}*.stats | \\
    tr '\\n' ', ' | \\
    awk -v sample=${meta.id} -v readnumber=\$readnumber '{ print "id: \\"my_pca_section\\"\\nsection_name: \\"Contamination Filtering\\"\\ndescription: \\"This plot shows the amount of reads filtered by contaminant type.\\"\\nplot_type: \\"bargraph\\"\\npconfig:\\n  id: \\"contamination_filter_plot\\"\\n  title: \\"Contamination Plot\\"\\n  ylab: \\"Number of reads\\"\\ndata:\\n    "sample": {"\$0"\\"remaining reads\\": "readnumber"}" }' > ${meta.id}.contamination_mqc.yaml

    if [[ ${reads} == *.gz ]]; then
        cp ${reads} ${meta.id}.filtered.fastq.gz
    else
        gzip -c ${reads} > ${meta.id}.filtered.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BusyBox: \$(busybox | sed -n -E 's/.*v([[:digit:].]+)\\s\\(.*/\\1/p')
    END_VERSIONS
    """
}
