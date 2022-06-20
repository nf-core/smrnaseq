def VERSION = '2.0.1'

process MIRDEEP2_RUN {
    label 'process_medium'
    errorStrategy 'ignore'

    conda (params.enable_conda ? 'bioconda::mirdeep2:2.0.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mirdeep2:2.0.1.3--hdfd78af_1' :
        'quay.io/biocontainers/mirdeep2:2.0.1.3--hdfd78af_1' }"

    input:
    path fasta
    tuple path(reads), path(arf)
    path hairpin
    path mature

    output:
    path 'result*.{bed,csv,html}', emit: result
    path "versions.yml"          , emit: versions

    script:
    """
    miRDeep2.pl \\
    $reads \\
    $fasta \\
    $arf \\
    $mature \\
    none \\
    $hairpin \\
    -d \\
    -z _${reads.simpleName}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}":
        mirdeep2: \$(echo "$VERSION")
    END_VERSIONS
    """
}
