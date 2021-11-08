// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process INDEX_GENOME {
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::bowtie=1.3.0-2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bowtie:1.3.0--py38hcf49a77_2"
    } else {
        container "quay.io/biocontainers/bowtie:1.3.0--py38hcf49a77_2"
    }

    input:
    path fasta

    output:
    path 'genome*ebwt' , emit: bt_indices
    path 'genome.edited.fa' , emit: fasta

    script:
    def software = getSoftwareName(task.process)

    """

    # Remove any special base characters from reference genome FASTA file
    sed '/^[^>]/s/[^ATGCatgc]/N/g' $fasta > genome.edited.fa
    sed -i 's/ .*//' genome.edited.fa

    # Build bowtie index
    bowtie-build genome.edited.fa genome --threads ${task.cpus}

    """

}
