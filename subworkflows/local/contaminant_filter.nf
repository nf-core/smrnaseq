//
// Filter contamination by rrna, trna, cdna, ncma, pirna
//

include { BLAT_MIRNA  as BLAT_CDNA
          BLAT_MIRNA  as BLAT_NCRNA
          BLAT_MIRNA  as BLAT_PIRNA } from '../../modules/local/blat_mirna'

include { INDEX_CONTAMINANTS as INDEX_RRNA
          INDEX_CONTAMINANTS as INDEX_TRNA
          INDEX_CONTAMINANTS as INDEX_CDNA
          INDEX_CONTAMINANTS as INDEX_NCRNA
          INDEX_CONTAMINANTS as INDEX_PIRNA } from '../../modules/local/bowtie_contaminants'

include { BOWTIE_MAP_CONTAMINANTS as MAP_RRNA
          BOWTIE_MAP_CONTAMINANTS as MAP_TRNA
          BOWTIE_MAP_CONTAMINANTS as MAP_CDNA
          BOWTIE_MAP_CONTAMINANTS as MAP_NCRNA
          BOWTIE_MAP_CONTAMINANTS as MAP_PIRNA } from '../../modules/local/bowtie_map_contaminants'

include { FILTER_STATS } from '../../modules/local/filter_stats'

workflow CONTAMINANT_FILTER {
    take:
    mirna
    rrna
    trna
    cdna
    ncrna
    pirna
    reads      // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()
    ch_filter_stats = Channel.empty()

    rrna_reads = reads

    reads
        .map { add_suffix(it, "rrna") }
        .dump (tag:'rrna')
        .set { rrna_reads }


    if (params.rrna) {
        // Index DB and filter $reads emit: $rrna_reads
        INDEX_RRNA ( rrna )
        ch_versions = ch_versions.mix(INDEX_RRNA.out.versions)
        MAP_RRNA ( rrna_reads, INDEX_RRNA.out.bt_indices, '"rRNA"' )
        ch_versions = ch_versions.mix(MAP_RRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_RRNA.out.stats.ifEmpty(null))
        MAP_RRNA.out.unmapped
            .map { add_suffix(it, "rrna") }
            .dump (tag:'rrna')
            .set { rrna_reads }
    }

    rrna_reads
        .map { add_suffix(it, "rrna") }
        .dump (tag:'rrna')
        .set { trna_reads }

    if (params.trna) {
        // Index DB and filter $rrna_reads emit: $trna_reads
        INDEX_TRNA ( trna )
        ch_versions = ch_versions.mix(INDEX_TRNA.out.versions)
        MAP_TRNA ( rrna_reads, INDEX_TRNA.out.bt_indices, '"tRNA"')
        ch_versions = ch_versions.mix(MAP_TRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_TRNA.out.stats.ifEmpty(null))
        MAP_TRNA.out.unmapped
            .map { add_suffix(it, "trna") }
            .dump (tag:'trna')
            .set { trna_reads }
    }

    trna_reads
        .map { add_suffix(it, "cdna") }
        .dump (tag:'cdna')
        .set { cdna_reads }


    if (params.cdna) {
        BLAT_CDNA ( 'cdna', mirna, cdna )
        ch_versions = ch_versions.mix(BLAT_CDNA.out.versions)
        INDEX_CDNA ( BLAT_CDNA.out.filtered_set )
        ch_versions = ch_versions.mix(INDEX_CDNA.out.versions)
        MAP_CDNA (  trna_reads, INDEX_CDNA.out.bt_indices, '"cDNA"' )
        ch_versions = ch_versions.mix(MAP_CDNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_CDNA.out.stats.ifEmpty(null))
        MAP_CDNA.out.unmapped
            .map { add_suffix(it, "cdna") }
            .dump (tag:'cdna')
            .set { cdna_reads }
    }

    cdna_reads
        .map { add_suffix(it, "ncrna") }
        .dump (tag:'ncrna')
        .set { ncrna_reads }

    if (params.ncrna) {
        BLAT_NCRNA ( 'ncrna', mirna, ncrna )
        ch_versions = ch_versions.mix(BLAT_NCRNA.out.versions)
        INDEX_NCRNA ( BLAT_NCRNA.out.filtered_set )
        ch_versions = ch_versions.mix(INDEX_NCRNA.out.versions)
        MAP_NCRNA ( cdna_reads, INDEX_NCRNA.out.bt_indices, '"ncRNA"' )
        ch_versions = ch_versions.mix(MAP_NCRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_NCRNA.out.stats.ifEmpty(null))
        MAP_NCRNA.out.unmapped
            .map { add_suffix(it, "ncrna") }
            .dump (tag:'ncrna')
            .set { ncrna_reads }
    }

    ncrna_reads
        .map { add_suffix(it, "pirna") }
        .dump (tag:'pirna')
        .set { pirna_reads }

    if (params.pirna) {
        BLAT_PIRNA ( 'other', mirna, pirna )
        ch_versions = ch_versions.mix(BLAT_PIRNA.out.versions)
        INDEX_PIRNA ( BLAT_PIRNA.out.filtered_set )
        ch_versions = ch_versions.mix(INDEX_PIRNA.out.versions)
        MAP_PIRNA (ncrna_reads, INDEX_PIRNA.out.bt_indices, '"piRNA"' )
        ch_versions = ch_versions.mix(MAP_PIRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_PIRNA.out.stats.ifEmpty(null))
        MAP_PIRNA.out.unmapped
            .map { add_suffix(it, "pirna") }
            .dump (tag:'pirna')
            .set { pirna_reads }
    }

    FILTER_STATS ( pirna_reads, ch_filter_stats.collect() )


    emit:
    filtered_reads = pirna_reads
    versions = ch_versions
    filter_stats = FILTER_STATS.out.stats
}

def add_suffix(row, suffix) {
    def meta = [:]
    meta.id           = "${row[0].id}_${suffix}"
    def array = []
    array = [ meta, row[1] ]
    return array
}