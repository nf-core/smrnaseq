process MIRTRACE_RUN {
    label 'process_medium'

    conda 'bioconda::mirtrace=1.0.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mirtrace:1.0.1--hdfd78af_1' :
        'biocontainers/mirtrace:1.0.1--hdfd78af_1' }"

    input:
    tuple val(adapter), val(ids), val(reads)

    output:
    path "mirtrace/*"  , emit: mirtrace
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // mirtrace protocol defaults to 'params.protocol' if not set
    def primer = adapter ? "--adapter ${adapter}" : ""
    def protocol = params.protocol == 'custom' ? '' : "--protocol $params.protocol"
    def java_mem = ''
    if(task.memory){
        tmem = task.memory.toBytes()
        java_mem = "-Xms${tmem} -Xmx${tmem}"
    }
    def config_lines = [ids,reads]
    .transpose()
    .collect({ id, path -> "echo '${path},${id}' >> mirtrace_config" })
    """
    export mirtracejar=\$(dirname \$(which mirtrace))

    ${config_lines.join("\n    ")}

    java $java_mem -jar \$mirtracejar/mirtrace.jar --mirtrace-wrapper-name mirtrace qc  \\
        --species $params.mirtrace_species \\
        $primer \\
        $protocol \\
        --config mirtrace_config \\
        --write-fasta \\
        --output-dir mirtrace \\
        --force

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtrace: \$(echo \$(mirtrace -v))
    END_VERSIONS
    """

}
