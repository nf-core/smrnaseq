//
// Quantify mirna with bowtie and mirtop
//

include {   PARSE_FASTA_MIRNA  as PARSE_MATURE
            PARSE_FASTA_MIRNA  as PARSE_HAIRPIN      } from '../../modules/local/parse_fasta_mirna'

include {   FORMAT_FASTA_MIRNA  as FORMAT_MATURE
            FORMAT_FASTA_MIRNA  as FORMAT_HAIRPIN    } from '../../modules/local/format_fasta_mirna'

include { BOWTIE_BUILD  as INDEX_MATURE     } from '../../modules/nf-core/bowtie/build/main'
include { BOWTIE_BUILD  as INDEX_HAIRPIN    } from '../../modules/nf-core/bowtie/build/main'

include {   BOWTIE_ALIGN  as BOWTIE_MAP_MATURE
            BOWTIE_ALIGN  as BOWTIE_MAP_HAIRPIN
            BOWTIE_ALIGN  as BOWTIE_MAP_SEQCLUSTER } from '../../modules/nf-core/bowtie/align/main'

include {   BAM_SORT_STATS_SAMTOOLS as BAM_STATS_MATURE
            BAM_SORT_STATS_SAMTOOLS as BAM_STATS_HAIRPIN   } from '../nf-core/bam_sort_stats_samtools'

include { SEQCLUSTER_COLLAPSE  } from '../../modules/nf-core/seqcluster/collapse/main'
include { TABLE_MERGE          } from '../../modules/local/datatable_merge/datatable_merge.nf'
include { EDGER_QC             } from '../../modules/local/edger_qc/edger_qc.nf'
include { BAM_STATS_MIRNA_MIRTOP } from '../../subworkflows/nf-core/bam_stats_mirna_mirtop/main'
include { CSVTK_JOIN             } from '../../modules/nf-core/csvtk/join/main'

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
    ch_mature_bowtie = INDEX_MATURE.out.index.first()
    ch_versions = ch_versions.mix(INDEX_MATURE.out.versions)

    ch_reads_mirna = ch_reads_for_mirna
        .map { add_suffix(it, "mature") }

    BOWTIE_MAP_MATURE ( ch_reads_mirna, ch_mature_bowtie, true )
    ch_versions = ch_versions.mix(BOWTIE_MAP_MATURE.out.versions)

    ch_reads_hairpin = BOWTIE_MAP_MATURE.out.fastq
        .map { add_suffix(it, "hairpin") }

    BAM_STATS_MATURE ( BOWTIE_MAP_MATURE.out.bam, FORMAT_MATURE.out.formatted_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_MATURE.out.versions)

    PARSE_HAIRPIN ( ch_reference_hairpin, ch_parse_species_input )
    ch_hairpin_parsed = PARSE_HAIRPIN.out.parsed_fasta
    ch_versions = ch_versions.mix(PARSE_HAIRPIN.out.versions)

    FORMAT_HAIRPIN ( ch_hairpin_parsed )
    ch_versions = ch_versions.mix(FORMAT_HAIRPIN.out.versions)

    INDEX_HAIRPIN ( FORMAT_HAIRPIN.out.formatted_fasta )
    hairpin_bowtie = INDEX_HAIRPIN.out.index.first()
    ch_versions = ch_versions.mix(INDEX_HAIRPIN.out.versions)

    BOWTIE_MAP_HAIRPIN ( ch_reads_hairpin, hairpin_bowtie, true )
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

    BOWTIE_MAP_SEQCLUSTER ( ch_reads_collapsed, hairpin_bowtie, true )
    ch_versions = ch_versions.mix(BOWTIE_MAP_SEQCLUSTER.out.versions)

    ch_mirtop_logs = Channel.empty()

    // nf-core/mirtop

    ch_mirna_gtf_species = ch_mirna_gtf
            .combine(ch_mirtrace_species)
            .map{ gtf, species -> [ [id:species.toString()], gtf, species ] }
            .collect()

    BAM_STATS_MIRNA_MIRTOP(
            BOWTIE_MAP_SEQCLUSTER.out.bam,
            FORMAT_HAIRPIN.out.formatted_fasta,
            ch_mirna_gtf_species )

    ch_mirtop_logs = BAM_STATS_MIRNA_MIRTOP.out.stats_log
    ch_versions = ch_versions.mix(BAM_STATS_MIRNA_MIRTOP.out.versions)

    ch_tsvs = BAM_STATS_MIRNA_MIRTOP.out.counts
        .collect{it[1]}
        .map{it -> return [[id:"TSVs"], it]}

    CSVTK_JOIN ( ch_tsvs )

    TABLE_MERGE ( CSVTK_JOIN.out.csv )
    ch_versions = ch_versions.mix(TABLE_MERGE.out.versions)

    ch_reads_genome = BOWTIE_MAP_HAIRPIN.out.fastq
        .map { add_suffix(it, "genome") }

    emit:
    fasta_mature        = FORMAT_MATURE.out.formatted_fasta // channel: [ val(meta), path(fasta) ]
    fasta_hairpin       = FORMAT_HAIRPIN.out.formatted_fasta // channel: [ val(meta), path(fasta) ]
    unmapped            = ch_reads_genome // channel: [ val(meta), path(bam) ]
    mature_stats        = BAM_STATS_MATURE.out.stats //TODO not used for antything, should we remove them?
    hairpin_stats       = BAM_STATS_HAIRPIN.out.stats //TODO not used for antything, should we remove them?
    mirtop_logs         = ch_mirtop_logs // channel: [ val(meta), path(log) ]
    versions            = ch_versions
}

def add_suffix(row, suffix) {
    def meta_clone = row[0].clone()
    meta_clone.id = "${row[0].id}_${suffix}"
    return [ meta_clone, row[1] ]
}
