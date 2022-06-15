// 
// Deduplicate the UMI reads by mapping them to the complete genome.
//

include { INDEX_GENOME                        } from '../../modules/local/bowtie_genome'
include { BOWTIE_MAP_SEQ as BOWTIE_MAP_GENOME } from '../../modules/local/bowtie_map_mirna'
include { BAM_SORT_SAMTOOLS                   } from '../../subworkflows/nf-core/bam_sort_samtools'
include { UMITOOLS_DEDUP                      } from '../../modules/nf-core/modules/umitools/dedup/main'
include { SAMTOOLS_BAM2FQ                     } from '../../modules/nf-core/modules/samtools/bam2fq/main'

workflow DEDUPLICATE_UMIS {
    take:
    fasta
    bt_index
    reads      // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()
    ch_dedup_stats = Channel.empty()

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

        reads.view()
        
        BOWTIE_MAP_GENOME ( reads, bt_indices.collect() )
        ch_versions = ch_versions.mix(BOWTIE_MAP_GENOME.out.versions)

        BAM_SORT_SAMTOOLS ( BOWTIE_MAP_GENOME.out.bam, Channel.empty() )
        ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

        BAM_SORT_SAMTOOLS.out.bam.view()
        ch_umi_dedup = BAM_SORT_SAMTOOLS.out.bam.join(BAM_SORT_SAMTOOLS.out.bai)

        ch_umi_dedup.view()

        UMITOOLS_DEDUP ( ch_umi_dedup )
        ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)
        ch_dedup_stats = ch_dedup_stats.mix(UMITOOLS_DEDUP.out.tsv_edit_distance).join(UMITOOLS_DEDUP.out.tsv_per_umi).join(UMITOOLS_DEDUP.out.tsv_umi_per_position)

        SAMTOOLS_BAM2FQ ( UMITOOLS_DEDUP.out.bam, false )
        ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions)
    }

    emit:
    reads    = SAMTOOLS_BAM2FQ.out.reads // channel: [ val(meta), [ reads ] ]
    indices  = bt_indices
    stats    = ch_dedup_stats
    versions = ch_versions
}


def add_suffix(row, suffix) {
    def meta = [:]
    meta.id           = "${row[0].id}_${suffix}"
    def array = []
    array = [ meta, row[1] ]
    return array
}