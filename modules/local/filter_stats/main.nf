process FILTER_STATS {
    label 'process_medium'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b244720eef0bd28a41d4f26e33d3800e75d9fc87f080e81d42d5d676b4960dc/data' :
        'community.wave.seqera.io/library/gawk:5.3.0--180f75ae8b0ce739' }"

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
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        gzip: \$(gzip --version | head -1 | cut -d ' ' -f 2)
        GNU coreutils: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
