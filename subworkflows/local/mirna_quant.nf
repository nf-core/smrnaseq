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
include { SEQCLUSTER_COLLAPSE  } from '../../modules/nf-core/seqcluster/collapse/main'
include { MIRTOP_QUANT         } from '../../modules/local/mirtop_quant.nf'
include { TABLE_MERGE          } from '../../modules/local/datatable_merge/datatable_merge.nf'
include { EDGER_QC             } from '../../modules/local/edger_qc/edger_qc.nf'

workflow MIRNA_QUANT {
    take:
    ch_reference_mature  // channel: [ val(meta), fasta file]
    ch_reference_hairpin // channel: [ val(meta), fasta file]
    ch_mirna_gtf         // channel: [ path(GTF) ]
    ch_reads_for_mirna   // channel: [ val(meta), [ reads ] ]
    ch_mirtrace_species  // channel: [ val(string) ]

    main:
    ch_versions = Channel.empty()
    ch_parse_species_input = params.mirgenedb ? Channel.value(params.mirgenedb_species) : ch_mirtrace_species

    PARSE_MATURE ( ch_reference_mature, ch_parse_species_input )
    ch_mirna_parsed = PARSE_MATURE.out.parsed_fasta
    ch_versions = ch_versions.mix(PARSE_MATURE.out.versions)

    FORMAT_MATURE ( ch_mirna_parsed )
    ch_versions = ch_versions.mix(FORMAT_MATURE.out.versions)

    INDEX_MATURE ( FORMAT_MATURE.out.formatted_fasta )
    ch_mature_bowtie = INDEX_MATURE.out.index
    ch_versions = ch_versions.mix(INDEX_MATURE.out.versions)

    ch_reads_mirna = ch_reads_for_mirna
        .map { add_suffix(it, "mature") }

    BOWTIE_MAP_MATURE ( ch_reads_mirna, ch_mature_bowtie.collect() )
    ch_versions = ch_versions.mix(BOWTIE_MAP_MATURE.out.versions)

    ch_reads_hairpin = BOWTIE_MAP_MATURE.out.unmapped
        .map { add_suffix(it, "hairpin") }

    BAM_STATS_MATURE ( BOWTIE_MAP_MATURE.out.bam, FORMAT_MATURE.out.formatted_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_MATURE.out.versions)

    PARSE_HAIRPIN ( ch_reference_hairpin, ch_parse_species_input )
    ch_hairpin_parsed = PARSE_HAIRPIN.out.parsed_fasta
    ch_versions = ch_versions.mix(PARSE_HAIRPIN.out.versions)

    FORMAT_HAIRPIN ( ch_hairpin_parsed )
    ch_versions = ch_versions.mix(FORMAT_HAIRPIN.out.versions)

    INDEX_HAIRPIN ( FORMAT_HAIRPIN.out.formatted_fasta )
    hairpin_bowtie = INDEX_HAIRPIN.out.index
    ch_versions = ch_versions.mix(INDEX_HAIRPIN.out.versions)

    BOWTIE_MAP_HAIRPIN ( ch_reads_hairpin, hairpin_bowtie.collect() )
    ch_versions = ch_versions.mix(BOWTIE_MAP_HAIRPIN.out.versions)

    BAM_STATS_HAIRPIN ( BOWTIE_MAP_HAIRPIN.out.bam, FORMAT_HAIRPIN.out.formatted_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_HAIRPIN.out.versions)

    ch_edger_input = BAM_STATS_MATURE.out.idxstats.collect{it[1]}
        .mix(BAM_STATS_HAIRPIN.out.idxstats.collect{it[1]})
        .flatten()
        .collect()

    EDGER_QC ( ch_edger_input )
    ch_versions.mix(EDGER_QC.out.versions)

    SEQCLUSTER_COLLAPSE ( ch_reads_for_mirna )
    ch_versions = ch_versions.mix(SEQCLUSTER_COLLAPSE.out.versions)

    ch_reads_collapsed = SEQCLUSTER_COLLAPSE.out.fastq

    BOWTIE_MAP_SEQCLUSTER ( ch_reads_collapsed, hairpin_bowtie.collect() )
    ch_versions = ch_versions.mix(BOWTIE_MAP_SEQCLUSTER.out.versions)

    ch_mirtop_logs = Channel.empty()
    //Block wont run if ch_mirtrace_species is emtpy, no need for conditional statement, TODO remove
    MIRTOP_QUANT ( BOWTIE_MAP_SEQCLUSTER.out.bam.collect{it[1]}, FORMAT_HAIRPIN.out.formatted_fasta.collect{it[1]}, ch_mirna_gtf, ch_mirtrace_species )
    ch_mirtop_logs = MIRTOP_QUANT.out.logs
    ch_versions = ch_versions.mix(MIRTOP_QUANT.out.versions)
    TABLE_MERGE ( MIRTOP_QUANT.out.mirtop_table )
    ch_versions = ch_versions.mix(TABLE_MERGE.out.versions)

    ch_reads_genome = BOWTIE_MAP_HAIRPIN.out.unmapped
        .map { add_suffix(it, "genome") }

    emit:
    fasta_mature        = FORMAT_MATURE.out.formatted_fasta // channel: [ val(meta), path(fasta) ]
    fasta_hairpin       = FORMAT_HAIRPIN.out.formatted_fasta // channel: [ val(meta), path(fasta) ]
    unmapped            = ch_reads_genome // channel: [ val(meta), path(bam) ]
    mature_stats        = BAM_STATS_MATURE.out.stats //TODO not used for antything, should we remove them?
    hairpin_stats       = BAM_STATS_HAIRPIN.out.stats //TODO not used for antything, should we remove them?
    mirtop_logs         = ch_mirtop_logs //TODO not used for antything, should we remove them?
    versions            = ch_versions
}

def add_suffix(row, suffix) {
    def meta = [:]
    meta.id = "${row[0].id}_${suffix}"
    def array = []
    array = [ meta, row[1] ]
    return array
}
