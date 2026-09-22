//
// Quantify mirna with bowtie and mirtop
//

include { BAM_SORT_STATS_SAMTOOLS           } from '../../nf-core/bam_sort_stats_samtools'
include { BOWTIE_ALIGN as BOWTIE_MAP_GENOME } from '../../../modules/nf-core/bowtie/align/main'

workflow GENOME_QUANT {
    take:
    ch_bowtie_index // channel: [ val(meta), path(directory_index) ]
    ch_fasta        // channel: [ val(meta), path(fasta) ]
    ch_reads        // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = channel.empty()
    ch_fasta_fai = ch_fasta.map { row ->
        row.size() == 2 ? [ row[0], row[1], [] ] : row
    }

    BOWTIE_MAP_GENOME ( ch_reads, ch_bowtie_index, true )

    BAM_SORT_STATS_SAMTOOLS ( BOWTIE_MAP_GENOME.out.bam,  ch_fasta_fai )

    emit:
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    versions = ch_versions                       // channel: [ versions.yml ]
}
