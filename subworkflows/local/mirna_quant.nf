//
// Quantify mirna with bowtie and mirtop
//

include {   PARSE_FASTA_MIRNA  as PARSE_MATURE
            PARSE_FASTA_MIRNA  as PARSE_HAIRPIN      } from '../../modules/local/parse_fasta_mirna'

include {   FORMAT_FASTA_MIRNA  as FORMAT_MATURE
            FORMAT_FASTA_MIRNA  as FORMAT_HAIRPIN    } from '../../modules/local/format_fasta_mirna'

include {   INDEX_MIRNA  as INDEX_MATURE
            INDEX_MIRNA  as INDEX_HAIRPIN            } from '../../modules/local/bowtie_mirna'

include {   BOWTIE_MAP_SEQ  as BOWTIE_MAP_MATURE
            BOWTIE_MAP_SEQ  as BOWTIE_MAP_HAIRPIN
            BOWTIE_MAP_SEQ  as BOWTIE_MAP_SEQCLUSTER } from '../../modules/local/bowtie_map_mirna'

include {   BAM_SORT_STATS_SAMTOOLS as BAM_STATS_MATURE
            BAM_SORT_STATS_SAMTOOLS as BAM_STATS_HAIRPIN   } from '../nf-core/bam_sort_stats_samtools'

include { SEQCLUSTER_SEQUENCES } from '../../modules/local/seqcluster_collapse.nf'
include { MIRTOP_QUANT         } from '../../modules/local/mirtop_quant.nf'
include { TABLE_MERGE          } from '../../modules/local/datatable_merge.nf'
include { EDGER_QC             } from '../../modules/local/edger_qc.nf'

workflow MIRNA_QUANT {
    take:
    mature     // channel: [ val(meta), fasta file]
    hairpin    // channel: [ val(meta), fasta file]
    gtf        // channle: GTF file
    reads      // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    PARSE_MATURE ( mature ).parsed_fasta.set { mirna_parsed }
    ch_versions = ch_versions.mix(PARSE_MATURE.out.versions)

    FORMAT_MATURE ( mirna_parsed )
    ch_versions = ch_versions.mix(FORMAT_MATURE.out.versions)

    INDEX_MATURE ( FORMAT_MATURE.out.formatted_fasta ).index.set { mature_bowtie }
    ch_versions = ch_versions.mix(INDEX_MATURE.out.versions)

    reads
        .map { add_suffix(it, "mature") }
        .dump (tag:'msux')
        .set { reads_mirna }

    BOWTIE_MAP_MATURE ( reads_mirna, mature_bowtie.collect() )
    ch_versions = ch_versions.mix(BOWTIE_MAP_MATURE.out.versions)

    BOWTIE_MAP_MATURE.out.unmapped
        .map { add_suffix(it, "hairpin") }
        .dump (tag:'hsux')
        .set { reads_hairpin }

    BAM_STATS_MATURE ( BOWTIE_MAP_MATURE.out.bam, FORMAT_MATURE.out.formatted_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_MATURE.out.versions)

    PARSE_HAIRPIN ( hairpin ).parsed_fasta.set { hairpin_parsed }
    ch_versions = ch_versions.mix(PARSE_HAIRPIN.out.versions)

    FORMAT_HAIRPIN ( hairpin_parsed )
    ch_versions = ch_versions.mix(FORMAT_HAIRPIN.out.versions)

    INDEX_HAIRPIN ( FORMAT_HAIRPIN.out.formatted_fasta ).index.set { hairpin_bowtie }
    ch_versions = ch_versions.mix(INDEX_HAIRPIN.out.versions)

    BOWTIE_MAP_HAIRPIN ( reads_hairpin, hairpin_bowtie.collect() )
    ch_versions = ch_versions.mix(BOWTIE_MAP_HAIRPIN.out.versions)

    BAM_STATS_HAIRPIN ( BOWTIE_MAP_HAIRPIN.out.bam, FORMAT_HAIRPIN.out.formatted_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_HAIRPIN.out.versions)

    BAM_STATS_MATURE.out.idxstats.collect{it[1]}
        .mix(BAM_STATS_HAIRPIN.out.idxstats.collect{it[1]})
        .dump(tag:'edger')
        .flatten()
        .collect()
        .set { edger_input }

    EDGER_QC ( edger_input )
    ch_versions.mix(EDGER_QC.out.versions)

    reads
        .map { add_suffix(it, "seqcluster") }
        .dump (tag:'ssux')
        .set { reads_seqcluster }

    SEQCLUSTER_SEQUENCES ( reads_seqcluster ).collapsed.set { reads_collapsed }
    ch_versions = ch_versions.mix(SEQCLUSTER_SEQUENCES.out.versions)

    BOWTIE_MAP_SEQCLUSTER ( reads_collapsed, hairpin_bowtie.collect() )
    ch_versions = ch_versions.mix(BOWTIE_MAP_SEQCLUSTER.out.versions)

    ch_mirtop_logs = Channel.empty()
    if (params.mirtrace_species){
        MIRTOP_QUANT ( BOWTIE_MAP_SEQCLUSTER.out.bam.collect{it[1]}, FORMAT_HAIRPIN.out.formatted_fasta.collect{it[1]}, gtf )
        ch_mirtop_logs = MIRTOP_QUANT.out.logs
        ch_versions = ch_versions.mix(MIRTOP_QUANT.out.versions)

        TABLE_MERGE ( MIRTOP_QUANT.out.mirtop_table )
        ch_versions = ch_versions.mix(TABLE_MERGE.out.versions)
    }
    BOWTIE_MAP_HAIRPIN.out.unmapped
        .map { add_suffix(it, "genome") }
        .dump (tag:'gsux')
        .set { reads_genome }

    emit:
    fasta_mature        = FORMAT_MATURE.out.formatted_fasta
    fasta_hairpin       = FORMAT_HAIRPIN.out.formatted_fasta
    unmapped            = reads_genome
    mature_stats        = BAM_STATS_MATURE.out.stats
    hairpin_stats       = BAM_STATS_HAIRPIN.out.stats
    mirtop_logs         = ch_mirtop_logs

    versions            = ch_versions
}

def add_suffix(row, suffix) {
    def meta = [:]
    meta.id = "${row[0].id}_${suffix}"
    def array = []
    array = [ meta, row[1] ]
    return array
}
