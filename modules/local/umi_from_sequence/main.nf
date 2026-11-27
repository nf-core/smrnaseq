nextflow.enable.types = true

//
// Remove the UMI that UMI_TO_SEQUENCE put at the start of each read sequence (and quality string)
// again, using the UMI in the read name (`@<name>_<UMI> ...`) for its length.
//
// Also writes UMI saturation statistics for the deduplicated reads: the number of reads (molecules)
// per insert sequence, and how much of the UMI space (4^UMI length) the most frequent sequence uses.
//
process UMI_FROM_SEQUENCE {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b244720eef0bd28a41d4f26e33d3800e75d9fc87f080e81d42d5d676b4960dc/data' :
        'community.wave.seqera.io/library/gawk:5.3.0--180f75ae8b0ce739' }"

    input:
    tuple(meta: Map, reads: Path)

    output:
    reads      = tuple(meta, file('*.umi_dedup.fastq.gz'))
    saturation = tuple(meta, file('*.umi_saturation.tsv'))

    topic:
    tuple(task.process, 'gawk', eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cd ${reads} \\
        | awk -v sample='${meta.id}' -v stats='${prefix}.umi_saturation.tsv' '
            NR % 4 == 1 { n = split(\$1, name, "_"); len = length(name[n]); if (len > umi_len) umi_len = len; print; next }
            NR % 4 == 2 { seq = substr(\$0, len + 1); count[seq]++; reads++; print seq; next }
            NR % 4 == 3 { print; next }
            { print substr(\$0, len + 1) }
            END {
                top = ""; top_count = 0; distinct = 0
                for (s in count) { distinct++; if (count[s] > top_count) { top = s; top_count = count[s] } }
                umi_space = 4 ^ umi_len
                fraction = umi_space > 0 ? top_count / umi_space : 0
                print "sample\\tumi_length\\treads\\tdistinct_sequences\\ttop_sequence\\ttop_sequence_reads\\tumi_space\\ttop_sequence_umi_fraction" > stats
                printf "%s\\t%d\\t%d\\t%d\\t%s\\t%d\\t%.0f\\t%.6f\\n", sample, umi_len, reads, distinct, top, top_count, umi_space, fraction > stats
            }
        ' \\
        | gzip -c > ${prefix}.umi_dedup.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -n | gzip -c > ${prefix}.umi_dedup.fastq.gz
    printf "sample\\tumi_length\\treads\\tdistinct_sequences\\ttop_sequence\\ttop_sequence_reads\\tumi_space\\ttop_sequence_umi_fraction\\n" > ${prefix}.umi_saturation.tsv
    """
}
