//
// Quantify mirna with bowtie and mirtop
//

params.samtools_options = [:]
params.map_options = [:]
params.samtools_sort_options  = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { INDEX_GENOME } from '../../modules/local/bowtie_genome'
include { MAP_MIRNA as MAP_GENOME } from '../../modules/local/bowtie_map_mirna' addParams( options: params.map_options)
include { SAMTOOLS_VIEW  as SAMTOOLS_VIEW_GENOME } from '../../modules/nf-core/modules/samtools/view/main' addParams( options: params.samtools_options )
include { BAM_SORT_SAMTOOLS as BAM_STATS_GENOME } from './bam_sort' addParams( sort_options: params.samtools_sort_options, index_options: params.samtools_index_options, stats_options: params.samtools_stats_options )

workflow GENOME_QUANT {
    take:
    fasta
    bt_index
    reads      // channel: [ val(meta), [ reads ] ]

    main:

    if (!bt_index){
        INDEX_GENOME( fasta )
        bt_indices = INDEX_GENOME.out.bt_indeces
        fasta_formatted = INDEX_GENOME.out.fasta
    } else {
        bt_indices = Channel.fromPath("${bt_index}**ebwt", checkIfExists: true).ifEmpty { exit 1, "Bowtie1 index directory not found: ${bt_index}" }
        fasta_formatted = fasta
    }
    // else {
    //   bt_indeces = Channel.empty()
    //}

    if (bt_indices){
        MAP_GENOME ( reads, bt_indices.collect() )
        SAMTOOLS_VIEW_GENOME ( MAP_GENOME.out.sam )
        BAM_STATS_GENOME ( SAMTOOLS_VIEW_GENOME.out.bam )

    }

    emit:
        fasta   = fasta_formatted
        indeces = bt_indices


}
