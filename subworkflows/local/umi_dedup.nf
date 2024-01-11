//
// Deduplicate the UMI reads by mapping them to the complete genome.
//

include { INDEX_GENOME                        } from '../../modules/local/bowtie_genome'
include { BOWTIE_MAP_SEQ as UMI_MAP_GENOME    } from '../../modules/local/bowtie_map_mirna'
include { BAM_SORT_STATS_SAMTOOLS             } from '../../subworkflows/nf-core/bam_sort_stats_samtools'
include { UMITOOLS_DEDUP                      } from '../../modules/nf-core/modules/umitools/dedup/main'
include { SAMTOOLS_BAM2FQ                     } from '../../modules/nf-core/modules/samtools/bam2fq/main'
include { CAT_CAT                             } from '../../modules/nf-core/modules/cat/cat/main'

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
        bt_index      = INDEX_GENOME.out.bowtie_indices
        fasta_formatted = INDEX_GENOME.out.fasta
        ch_versions     = ch_versions.mix(INDEX_GENOME.out.versions)
    } else {
        bt_index     = Channel.fromPath("${bt_index}**ebwt", checkIfExists: true).ifEmpty { exit 1, "Bowtie1 index directory not found: ${bt_index}" }
        fasta_formatted = fasta
    }

    if (bt_index){

        UMI_MAP_GENOME ( reads, bt_index.collect() )
        ch_versions = ch_versions.mix(UMI_MAP_GENOME.out.versions)

        BAM_SORT_SAMTOOLS ( UMI_MAP_GENOME.out.bam, Channel.empty() )
        ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

        ch_umi_dedup = BAM_SORT_SAMTOOLS.out.bam.join(BAM_SORT_SAMTOOLS.out.bai)
        UMITOOLS_DEDUP ( ch_umi_dedup )
        ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)
        ch_dedup_stats = ch_dedup_stats.mix(UMITOOLS_DEDUP.out.tsv_edit_distance).join(UMITOOLS_DEDUP.out.tsv_per_umi).join(UMITOOLS_DEDUP.out.tsv_umi_per_position)

        SAMTOOLS_BAM2FQ ( UMITOOLS_DEDUP.out.bam, false )
        ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions)

        ch_dedup_reads = SAMTOOLS_BAM2FQ.out.reads

        if ( params.umi_merge_unmapped ) {

            SAMTOOLS_BAM2FQ.out.reads
                .join(UMI_MAP_GENOME.out.unmapped)
                .map { meta, file1, file2 -> [meta, [file1, file2]]}
                .set { ch_cat }

            CAT_CAT ( ch_cat )
            ch_dedup_reads = CAT_CAT.out.file_out
        }
    }

    emit:
    reads    = ch_dedup_reads
    indices  = bt_index
    stats    = ch_dedup_stats
    versions = ch_versions
}
