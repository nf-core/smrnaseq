//
// Quantify mirna with bowtie and mirtop
//

include { BAM_SORT_STATS_SAMTOOLS           } from '../nf-core/bam_sort_stats_samtools'
include { BOWTIE_ALIGN as BOWTIE_MAP_GENOME } from '../../modules/nf-core/bowtie/align/main'

workflow GENOME_QUANT {
    take:
    ch_bowtie_index // channel: [genome.1.ebwt, genome.2.ebwt, genome.3.ebwt, genome.4.ebwt, genome.rev.1.ebwt, genome.rev.2.ebwt]
    ch_fasta        // channel: [ val(meta), path(fasta) ]
    ch_reads        // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    BOWTIE_MAP_GENOME ( ch_reads, ch_bowtie_index, true )
    ch_versions = ch_versions.mix(BOWTIE_MAP_GENOME.out.versions)

    BAM_SORT_STATS_SAMTOOLS ( BOWTIE_MAP_GENOME.out.bam,  ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    fasta    = ch_fasta //TODO: This fasta is the same one that was used as input, ask the original developer, if they meant to have something else here
    index    = ch_bowtie_index //TODO: Same here, are we outputting the right files? We can remove these channels if we are.
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats

    versions = ch_versions
}
