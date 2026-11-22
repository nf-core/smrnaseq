nextflow.enable.types = true

// PROTOTYPE (#535): there is no nf-core/modules fgumi/copyumi module yet.
// Copies the UMI that `umi_tools extract` appended to the read name (`<name>_<UMI>`)
// into the RX tag, which `fgumi sort`/`fgumi dedup` require.
process FGUMI_COPYUMI {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8b/8b1b1f653fc1fef5eb7d31f72d7cc60e29ff0eec5552802523ce97b1b5d0f920/data' :
        'community.wave.seqera.io/library/fgumi:0.9.0--1f0a2f2a31a251a3' }"

    input:
    tuple(meta: Map, bam: Path)

    output:
    bam     = tuple(meta, file("${task.ext.prefix ?: "${meta.id}.rx"}.bam"))
    metrics = tuple(meta, file('*.copy_umi.tsv'))

    topic:
    tuple(task.process, 'fgumi', eval('fgumi --version | sed "s/^fgumi //"')) >> 'versions'

    script:
    def args   = task.ext.args ?: '--field-delimiter _'
    def prefix = task.ext.prefix ?: "${meta.id}.rx"
    """
    fgumi copy-umi \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --metrics ${prefix}.copy_umi.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.rx"
    """
    touch ${prefix}.bam
    printf 'total_records\\trx_written\\trx_overwritten\\tnames_trimmed\\n0\\t0\\t0\\t0\\n' > ${prefix}.copy_umi.tsv
    """
}
