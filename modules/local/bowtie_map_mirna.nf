// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BOWTIE_MAP_SEQ {
    label 'process_medium'
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"getSoftwareName(task.process)/${suffix}", meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::bowtie=1.3.0-2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bowtie:1.3.0--py38hcf49a77_2"
    } else {
        container "quay.io/biocontainers/bowtie:1.3.0--py38hcf49a77_2"
    }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*sam"), emit: sam
    tuple val(meta), path('unmapped/*fq.gz') , emit: unmapped
    path "*.version.txt" , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def process_name = task.process.tokenize(':')[-1]
    def index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    """

    bowtie \\
        -x $index_base \\
        -q <(zcat $reads) \\
        -p ${task.cpus} \\
        -t \\
        -k 50 \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        --un ${meta.id}_unmapped.fq -S > ${meta.id}.sam

    if [ ! -f  "${meta.id}_unmapped.fq" ]
    then
        touch ${meta.id}_unmapped.fq
    fi
    gzip ${meta.id}_unmapped.fq
    mkdir unmapped
    mv  ${meta.id}_unmapped.fq.gz  unmapped/.
    bowtie --version 2>&1 | head -1 | sed 's/^.*version //' > ${software}.version.txt

    """

}
