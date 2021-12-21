process TRIMGALORE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::trim-galore=0.6.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.7--hdfd78af_0' :
        'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz")    , emit: reads
    tuple val(meta), path("*report.txt"), emit: log
    path "versions.yml"                 , emit: versions

    tuple val(meta), path("*.html"), emit: html optional true
    tuple val(meta), path("*.zip") , emit: zip optional true

    script:
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }
    def prefix   = "${meta.id}"
    // Define regular variables so that they can be overwritten
    def clip_r1 = params.clip_r1
    def three_prime_clip_r1 = params.three_prime_clip_r1
    def three_prime_adapter = params.three_prime_adapter
    def protocol = params.protocol
    // Presets
    if (params.protocol == "illumina"){
        clip_r1 = 0
        three_prime_clip_r1 = 0
        three_prime_adapter = "TGGAATTCTCGGGTGCCAAGG"
    } else if (params.protocol == "nextflex"){
        clip_r1 = 4
        three_prime_clip_r1 = 4
        three_prime_adapter = "TGGAATTCTCGGGTGCCAAGG"
    } else if (params.protocol == "qiaseq"){
        clip_r1 = 0
        three_prime_clip_r1 = 0
        three_prime_adapter = "AACTGTAGGCACCATCAAT"
    } else if (params.protocol == "cats"){
        clip_r1 = 3
        three_prime_clip_r1 = 0
        // three_prime_adapter = "GATCGGAAGAGCACACGTCTG"
        three_prime_adapter = "AAAAAAAA"
    } else {
        //custom protocol
        clip_r1 = params.clip_r1
        three_prime_clip_r1 = params.three_prime_clip_r1
        three_prime_adapter = params.three_prime_adapter
        protocol = params.protocol
    }
    def tg_length = "--length ${params.min_length}"
    def c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    def tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    """
    [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
    trim_galore --cores $cores --adapter ${three_prime_adapter} $tg_length $c_r1 $tpc_r1 --max_length ${params.trim_galore_max_length} --gzip ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    ${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
