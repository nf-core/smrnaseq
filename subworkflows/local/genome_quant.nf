//
// Quantify mirna with bowtie and mirtop
//

include { BAM_SORT_STATS_SAMTOOLS           } from '../nf-core/bam_sort_stats_samtools'
include { BOWTIE_ALIGN as BOWTIE_MAP_GENOME } from '../../modules/nf-core/bowtie/align/main'

workflow GENOME_QUANT {
    take:
    ch_bowtie_index // channel: [ val(meta), path(directory_index) ]
    ch_fasta        // channel: [ val(meta), path(fasta) ]
    ch_reads        // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    BOWTIE_MAP_GENOME ( ch_reads, ch_bowtie_index, true )
    ch_versions = ch_versions.mix(BOWTIE_MAP_GENOME.out.versions)

    BAM_SORT_STATS_SAMTOOLS ( BOWTIE_MAP_GENOME.out.bam,  ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    versions = ch_versions                       // channel: [ versions.yml ]
}
