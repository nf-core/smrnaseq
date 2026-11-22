nextflow.enable.types = true

//
// Move the UMI that umi_tools extract wrote into the read name (`@<name>_<UMI> ...`) to the start of
// the read sequence, so that UMICollapse in FastQ mode deduplicates on UMI + insert.
// The UMI bases get the quality 'I'. UMI_FROM_SEQUENCE reverses this.
//
process UMI_TO_SEQUENCE {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b244720eef0bd28a41d4f26e33d3800e75d9fc87f080e81d42d5d676b4960dc/data' :
        'community.wave.seqera.io/library/gawk:5.3.0--180f75ae8b0ce739' }"

    input:
    tuple(meta: Map, reads: Path)

    output:
    reads = tuple(meta, file('*.umi_in_sequence.fastq.gz'))

    topic:
    tuple(task.process, 'gawk', eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cd ${reads} \\
        | awk '
            NR % 4 == 1 { n = split(\$1, name, "_"); umi = name[n]; print; next }
            NR % 4 == 2 { print umi \$0; next }
            NR % 4 == 3 { print; next }
            { qual = umi; gsub(/./, "I", qual); print qual \$0 }
        ' \\
        | gzip -c > ${prefix}.umi_in_sequence.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -n | gzip -c > ${prefix}.umi_in_sequence.fastq.gz
    """
}
