
process MIRTOP_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/28/28ece5ab35c2432bf6f360682f58d4245aec76a0cbab3879478f44d248df0205/data' :
        'community.wave.seqera.io/library/pybedtools_pysam_samtools_biopython_pruned:8e04862200c8021c'}"

    input:
    tuple val(meta), path(mirtop_gff)

    output:
    tuple val(meta), path("stats/*.txt")        , emit: txt
    tuple val(meta), path("stats/*_stats.log")  , emit: log
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mirtop \\
        stats \\
        $args \\
        --out stats \\
        $mirtop_gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir stats
    touch stats/${prefix}.txt
    touch stats/${prefix}_stats.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """
}
