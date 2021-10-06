// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process INDEX_MIRNA {
    label 'process_medium'

    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::bowtie=1.3.0-2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bowtie:1.3.0--py38hcf49a77_2"
    } else {
        container "quay.io/biocontainers/bowtie:1.3.0--py38hcf49a77_2"
    }

    input:
    path fasta

    output:
    path 'fasta_bidx*' , emit: bt_indeces
    path "*.version.txt" , emit: version

    script:
    def software = getSoftwareName(task.process)
    """

    bowtie-build ${fasta} fasta_bidx --threads ${task.cpus}

    bowtie --version 2>&1 | head -1 | sed 's/^.*version //' > ${software}.version.txt
    """

}