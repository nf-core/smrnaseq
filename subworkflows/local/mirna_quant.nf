//
// Quantify mirna with bowtie and mirtop
//

params.samtools_options = [:]
params.map_options = [:]
params.samtools_sort_options  = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { PARSE_FASTA_MIRNA  as PARSE_MATURE
          PARSE_FASTA_MIRNA  as PARSE_HAIRPIN        } from '../../modules/local/parse_fasta_mirna'

include { FORMAT_FASTA_MIRNA  as FORMAT_MATURE
          FORMAT_FASTA_MIRNA  as FORMAT_HAIRPIN        } from '../../modules/local/format_fasta_mirna'

include { INDEX_MIRNA  as INDEX_MATURE
          INDEX_MIRNA  as INDEX_HAIRPIN        } from '../../modules/local/bowtie_mirna'

include { MAP_MIRNA  as MAP_MATURE
          MAP_MIRNA  as MAP_HAIRPIN
          MAP_MIRNA  as MAP_SEQCLUSTER        } from '../../modules/local/bowtie_map_mirna' addParams(options: params.map_options)

include { SAMTOOLS_VIEW  as SAMTOOLS_VIEW_MATURE
          SAMTOOLS_VIEW  as SAMTOOLS_VIEW_HAIRPIN
          SAMTOOLS_VIEW  as SAMTOOLS_VIEW_SEQCLUSTER  } from '../../modules/nf-core/modules/samtools/view/main' addParams( options: params.samtools_options )

include { BAM_SORT_SAMTOOLS as BAM_STATS_MATURE
          BAM_SORT_SAMTOOLS as BAM_STATS_HAIRPIN } from './bam_sort' addParams( sort_options: params.samtools_sort_options, index_options: params.samtools_index_options, stats_options: params.samtools_stats_options )

include { SEQCLUSTER_SEQUENCES } from '../../modules/local/seqcluster_collapse.nf'
include { MIRTOP_QUANT } from '../../modules/local/mirtop_quant.nf'


workflow MIRNA_QUANT {
    take:
    mature     // channel: fasta file
    hairpin    // channel: fasta file
    gtf        // channle: GTF file
    reads      // channel: [ val(meta), [ reads ] ]

    main:
    PARSE_MATURE ( mature ).parsed_fasta.set { mirna_parsed }
    FORMAT_MATURE ( mirna_parsed )

    PARSE_HAIRPIN ( hairpin ).parsed_fasta.set { hairpin_parsed }
    FORMAT_HAIRPIN ( hairpin_parsed )

    INDEX_MATURE ( FORMAT_MATURE.out.formatted_fasta ).bt_indeces.set { mature_bowtie }
    INDEX_HAIRPIN ( FORMAT_HAIRPIN.out.formatted_fasta ).bt_indeces.set { hairpin_bowtie }

    reads
        .map { add_suffix(it, "mature") }
        .dump (tag:'msux')
        .set { reads_mirna }

    MAP_MATURE ( reads_mirna, mature_bowtie.collect() )
    SAMTOOLS_VIEW_MATURE ( MAP_MATURE.out.sam )

    MAP_MATURE.out.unmapped
        .map { add_suffix(it, "hairpin") }
        .dump (tag:'hsux')
        .set { reads_hairpin }

    MAP_HAIRPIN ( reads_hairpin, hairpin_bowtie.collect() )
    SAMTOOLS_VIEW_HAIRPIN ( MAP_HAIRPIN.out.sam )

    BAM_STATS_MATURE ( SAMTOOLS_VIEW_MATURE.out.bam )
    BAM_STATS_HAIRPIN ( SAMTOOLS_VIEW_HAIRPIN.out.bam )

    reads
        .map { add_suffix(it, "seqcluster") }
        .dump (tag:'ssux')
        .set { reads_seqcluster }

    SEQCLUSTER_SEQUENCES ( reads_seqcluster ).collapsed.set { reads_collapsed }
    MAP_SEQCLUSTER ( reads_collapsed, hairpin_bowtie.collect() )
    SAMTOOLS_VIEW_SEQCLUSTER ( MAP_SEQCLUSTER.out.sam )

    if (params.mirtrace_species){
        MIRTOP_QUANT ( SAMTOOLS_VIEW_SEQCLUSTER.out.bam.collect{it[1]}, FORMAT_HAIRPIN.out.formatted_fasta, gtf )
    }
    MAP_HAIRPIN.out.unmapped
        .map { add_suffix(it, "genome") }
        .dump (tag:'gsux')
        .set { reads_genome }

    emit:
    fasta_mature     = FORMAT_MATURE.out.formatted_fasta
    fasta_hairpin    = FORMAT_HAIRPIN.out.formatted_fasta
    unmapped         = reads_genome
    bowtie_version   = INDEX_MATURE.out.version
    samtools_version = BAM_STATS_MATURE.out.version
    seqcluster_version     = SEQCLUSTER_SEQUENCES.out.version
    mirtop_version   = MIRTOP_QUANT.out.version

}


def add_suffix(row, suffix) {
    def meta = [:]
    meta.id           = "${row[0].id}_${suffix}"
    def array = []
    array = [ meta, row[1] ]
    return array
}