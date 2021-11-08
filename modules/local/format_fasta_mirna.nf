// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FORMAT_FASTA_MIRNA {
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::fastx_toolkit=0.0.14-9' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fastx_toolkit:0.0.14--he1b5a44_8"
    } else {
        container "quay.io/biocontainers/fastx_toolkit:0.0.14--he1b5a44_8"
    }

    input:
    path fasta

    output:
    path '*_idx.fa' , emit: formatted_fasta
    //path "*.version.txt" , emit: versions

    script:
    def software = getSoftwareName(task.process)
    """
    fasta_formatter -w 0 -i $fasta -o ${fasta}_idx.fa

    #echo \$(fasta_formatter --version 2>&1) | sed 's/^.*version //; s/Last.*\$//' > ${software}.version.txt
    """

}
