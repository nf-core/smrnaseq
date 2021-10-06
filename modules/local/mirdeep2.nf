// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MIRDEEP2_PIGZ {
    label 'process_low'
    tag "$meta.id"

    conda (params.enable_conda ? 'bioconda::bioconvert:0.4.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconvert:0.4.3--py_0"
    } else {
        // container "nfcore/smrnaseq:dev"
        container "quay.io/biocontainers/bioconvert:0.4.3--py_0"
    }
    when:
    !params.skip_mirdeep

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq"), emit: reads

    script:
    def unzip = reads.toString() - '.gz'
    """
    pigz -f -d -p $task.cpus $reads
    """

}

process MIRDEEP2_MAPPER {
    label 'process_medium'
    tag "$meta.id"
    //publishDir "${params.outdir}/mirdeep/${meta.id}"

    conda (params.enable_conda ? 'bioconda::mirdeep2:2.0.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mirdeep2:2.0.1.3--hdfd78af_1"
    } else {
        // container "nfcore/smrnaseq:dev"
        container "quay.io/biocontainers/mirdeep2:2.0.1.3--hdfd78af_1"
    }

    when:
    !params.skip_mirdeep

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple path('*_collapsed.fa'), path('*reads_vs_refdb.arf'), emit: mirdeep2_inputs
    path "*.version.txt" , emit: version

    script:
    def index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    def software = getSoftwareName(task.process)
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
    echo "2.0.1" > ${software}.version.txt
    """
}

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
        // container "nfcore/smrnaseq:dev"
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
    """
}

