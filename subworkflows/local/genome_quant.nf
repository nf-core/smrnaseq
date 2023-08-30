//
// Quantify mirna with bowtie and mirtop
//

include { INDEX_GENOME      } from '../../modules/local/bowtie_genome'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools'
include { BOWTIE_MAP_SEQ as BOWTIE_MAP_GENOME } from '../../modules/local/bowtie_map_mirna'

workflow GENOME_QUANT {
    take:
    fasta
    index
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    if (!index){
        INDEX_GENOME ( fasta )
        bowtie_index    = INDEX_GENOME.out.index
        fasta_formatted = INDEX_GENOME.out.fasta
        ch_versions     = ch_versions.mix(INDEX_GENOME.out.versions)
    } else {
        bowtie_index    = Channel.fromPath("${index}**ebwt", checkIfExists: true).ifEmpty { exit 1, "Bowtie1 index directory not found: ${index}" }
        fasta_formatted = fasta
    }

    if (bowtie_index){
        BOWTIE_MAP_GENOME ( reads, bowtie_index.collect() )
        ch_versions = ch_versions.mix(BOWTIE_MAP_GENOME.out.versions)

        BAM_SORT_STATS_SAMTOOLS ( BOWTIE_MAP_GENOME.out.bam, Channel.empty()  )
        ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)
    }

    emit:
    fasta    = fasta_formatted
    index    = bowtie_index
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats

    versions = ch_versions
}
