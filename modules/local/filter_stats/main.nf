nextflow.enable.types = true

process FILTER_STATS {
    label 'process_medium'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b244720eef0bd28a41d4f26e33d3800e75d9fc87f080e81d42d5d676b4960dc/data' :
        'community.wave.seqera.io/library/gawk:5.3.0--180f75ae8b0ce739' }"

    input:
    tuple(meta: Map, reads: Path, stats_files: Set<Path>)

    output:
    stats = file('*_mqc.yaml')
    reads = tuple(meta, file('*.filtered.fastq.gz'))

    topic:
    tuple(task.process, 'gawk', eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'")) >> 'versions'
    tuple(task.process, 'gzip', eval("gzip --version | head -1 | cut -d ' ' -f 2")) >> 'versions'
    tuple(task.process, 'GNU coreutils', eval("cat --version | head -1 | sed 's/^.*coreutils) //; s/ .*\$//'")) >> 'versions'

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
    """

    stub:
    """
    cat <<-END_YAML > ${meta.id}.contamination_mqc.yaml
    id: "my_pca_section"
    section_name: "Contamination Filtering"
    description: "This plot shows the amount of reads filtered by contaminant type."
    plot_type: "bargraph"
    pconfig:
      id: "contamination_filter_plot"
      title: "Contamination Plot"
      ylab: "Number of reads"
    data:
        "${meta.id}": {"remaining reads": 0}
    END_YAML

    echo -n | gzip -c > ${meta.id}.filtered.fastq.gz
    """
}
