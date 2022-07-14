process MIRTRACE_RUN {
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::mirtrace=1.0.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mirtrace:1.0.1--hdfd78af_1' :
        'quay.io/biocontainers/mirtrace:1.0.1--hdfd78af_1' }"

    input:
    path reads

    output:
    path "mirtrace/*"  , emit: mirtrace
    path "versions.yml", emit: versions

    script:
    def three_prime_adapter = params.three_prime_adapter
    // Presets
    if (params.protocol == "illumina"){
        three_prime_adapter = "TGGAATTCTCGGGTGCCAAGG"
    } else if (params.protocol == "nextflex"){
        three_prime_adapter = "TGGAATTCTCGGGTGCCAAGG"
    } else if (params.protocol == "qiaseq"){
        three_prime_adapter = "AACTGTAGGCACCATCAAT"
    } else if (params.protocol == "cats"){
        three_prime_adapter = "AAAAAAAA"
    } else if (params.protocol == "custom"){
        three_prime_adapter = params.three_prime_adapter
    }

    // mirtrace protocol defaults to 'params.protocol' if not set
    def primer = (params.protocol=="cats") ? " " : " --adapter $three_prime_adapter "
    def protocol = (params.protocol=="custom") ? " " : "--protocol $params.protocol"
    def java_mem = ''
    if(task.memory){
        tmem = task.memory.toBytes()
        java_mem = "-Xms${tmem} -Xmx${tmem}"
    }
    """
    export mirtracejar=\$(dirname \$(which mirtrace))
    for i in $reads
    do
        path=\$(realpath \${i})
        prefix=\$(echo \${i} | sed -e 's/.gz//' -e 's/.fastq//' -e 's/.fq//' -e 's/_val_1//' -e 's/_trimmed//' -e 's/_R1//' -e 's/.R1//')
        echo \$path","\$prefix
    done > mirtrace_config

    java $java_mem -jar \$mirtracejar/mirtrace.jar --mirtrace-wrapper-name mirtrace qc  \\
        --species $params.mirtrace_species \\
        $primer \\
        $protocol \\
        --config mirtrace_config \\
        --write-fasta \\
        --output-dir mirtrace \\
        --force

    cat <<-END_VERSIONS > versions.yml
    ${task.process}":
        mirtrace: \$(echo \$(mirtrace -v 2>&1))
    END_VERSIONS
    """

}
