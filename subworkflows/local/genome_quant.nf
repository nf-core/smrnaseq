//
// Quantify mirna with bowtie and mirtop
//

// params.samtools_options       = [:]
// params.map_options            = [:]
// params.samtools_sort_options  = [:]
// params.samtools_index_options = [:]
// params.samtools_stats_options = [:]

include { INDEX_GENOME      } from '../../modules/local/bowtie_genome'
include { BAM_SORT_SAMTOOLS } from '../nf-core/bam_sort_samtools'      //addParams( sort_options: params.samtools_sort_options, index_options: params.samtools_index_options, stats_options: params.samtools_stats_options )
include { BOWTIE_MAP_SEQ as BOWTIE_MAP_GENOME } from '../../modules/local/bowtie_map_mirna' //addParams( options: params.map_options)

workflow GENOME_QUANT {
    take:
    fasta
    bt_index
    reads      // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()

    if (!bt_index){
        INDEX_GENOME ( fasta )
        bt_indices      = INDEX_GENOME.out.bt_indices
        fasta_formatted = INDEX_GENOME.out.fasta
        ch_versions     = ch_versions.mix(INDEX_GENOME.out.versions)
    } else {
        bt_indices      = Channel.fromPath("${bt_index}**ebwt", checkIfExists: true).ifEmpty { exit 1, "Bowtie1 index directory not found: ${bt_index}" }
        fasta_formatted = fasta
    }

    if (bt_indices){
        BOWTIE_MAP_GENOME ( reads, bt_indices.collect() )
        ch_versions = ch_versions.mix(BOWTIE_MAP_GENOME.out.versions)

        BAM_SORT_SAMTOOLS ( BOWTIE_MAP_GENOME.out.bam, Channel.empty()  )
        ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)
    }

    emit:
    fasta    = fasta_formatted
    indices  = bt_indices
    stats    = BAM_SORT_SAMTOOLS.out.stats

    versions = ch_versions
}
