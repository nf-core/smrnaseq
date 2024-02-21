//
// Quantify mirna with bowtie and mirtop
//

include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools'
include { BOWTIE_MAP_SEQ as BOWTIE_MAP_GENOME } from '../../modules/local/bowtie_map_mirna'

workflow GENOME_QUANT {
    take:
    bowtie_index
    fasta_formatted // fasta as generated by bowtie index step
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    BOWTIE_MAP_GENOME ( reads, bowtie_index.collect() )
    ch_versions = ch_versions.mix(BOWTIE_MAP_GENOME.out.versions)

    ch_fasta_formatted_for_sort = fasta_formatted .map { file -> tuple(file.baseName, file) }
    BAM_SORT_STATS_SAMTOOLS ( BOWTIE_MAP_GENOME.out.bam,  ch_fasta_formatted_for_sort )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    fasta    = fasta_formatted
    index    = bowtie_index
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats

    versions = ch_versions
}
