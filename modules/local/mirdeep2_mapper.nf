def VERSION = '2.0.1'

process MIRDEEP2_MAPPER {
    label 'process_medium'
    tag "$meta.id"

    conda 'bioconda::mirdeep2=2.0.1.3'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mirdeep2:2.0.1.3--hdfd78af_1' :
        'biocontainers/mirdeep2:2.0.1.3--hdfd78af_1' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple path('*_collapsed.fa'), path('*reads_vs_refdb.arf'), emit: mirdeep2_inputs
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    """
    mapper.pl \\
    $reads \\
        -e \\
        -h \\
        -i \\
        -j \\
        -m \\
        -p $index_base \\
        -s ${meta.id}_collapsed.fa \\
        -t ${meta.id}_reads_vs_refdb.arf \\
        -o 4

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapper: \$(echo "$VERSION")
    END_VERSIONS
    """
}
