// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MIRDEEP2_MAPPER {
    label 'process_medium'
    tag "$meta.id"
    //publishDir "${params.outdir}/mirdeep/${meta.id}"

    conda (params.enable_conda ? 'bioconda::mirdeep2=2.0.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mirdeep2:2.0.1.3--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/mirdeep2:2.0.1.3--hdfd78af_1"
    }

    when:
    !params.skip_mirdeep

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple path('*_collapsed.fa'), path('*reads_vs_refdb.arf'), emit: mirdeep2_inputs
    path "*.version.txt" , emit: versions

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
