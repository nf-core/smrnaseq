// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MIRDEEP2_RUN {
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::mirdeep2:2.0.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mirdeep2:2.0.1.3--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/mirdeep2:2.0.1.3--hdfd78af_1"
    }

    when:
    !params.skip_mirdeep

    input:
    path fasta
    tuple path(reads), path(arf)
    path hairpin
    path mature

    output:
    path 'result*.{bed,csv,html}'
    path "*.version.txt" , emit: versions

    script:
    def software = getSoftwareName(task.process)
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
    echo "2.0.1" > ${software}.version.txt
    """
}

