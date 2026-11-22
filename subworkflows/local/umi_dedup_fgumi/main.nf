//
// PROTOTYPE (#535): alignment-based UMI deduplication with fgumi.
//
// Input reads must already carry the UMI in the read name (`umi_tools extract`, `<name>_<UMI>`).
// Steps: bowtie1 (-k 1, one record per read, required by fgumi) -> fgumi copy-umi (name -> RX)
//        -> fgumi sort (template-coordinate) -> fgumi dedup (--remove-duplicates) -> fgumi fastq.
// The deduplicated FASTQ keeps the `_<UMI>` read-name suffix, so all downstream steps
// (mirtrace, MIRNA_QUANT, GENOME_QUANT, miRDeep2) can consume it unchanged.
//

include { BOWTIE_ALIGN as BOWTIE_MAP_UMI_DEDUP } from '../../../modules/nf-core/bowtie/align/main'
include { FGUMI_COPYUMI                        } from '../../../modules/local/fgumi_copyumi/main'
include { FGUMI_SORT                           } from '../../../modules/nf-core/fgumi/sort/main'
include { FGUMI_DEDUP                          } from '../../../modules/nf-core/fgumi/dedup/main'
include { FGUMI_FASTQ                          } from '../../../modules/nf-core/fgumi/fastq/main'

workflow UMI_DEDUP_FGUMI {
    take:
    ch_reads        // channel: [ val(meta), path(reads) ]  UMI in read name
    ch_bowtie_index // channel: [ val(meta2), path(index) ]  e.g. the formatted hairpin or mature index

    main:
    BOWTIE_MAP_UMI_DEDUP ( ch_reads, ch_bowtie_index.first(), false )

    FGUMI_COPYUMI ( BOWTIE_MAP_UMI_DEDUP.out.bam )

    FGUMI_SORT ( FGUMI_COPYUMI.out.bam )

    FGUMI_DEDUP ( FGUMI_SORT.out.bam )

    FGUMI_FASTQ ( FGUMI_DEDUP.out.bam )

    emit:
    reads     = FGUMI_FASTQ.out.reads        // channel: [ val(meta), path(fastq.gz) ]
    bam       = FGUMI_DEDUP.out.bam          // channel: [ val(meta), path(bam) ]
    metrics   = FGUMI_DEDUP.out.metrics      // channel: [ val(meta), path(metrics.txt) ]
    histogram = FGUMI_DEDUP.out.histogram    // channel: [ val(meta), path(family_size_histogram.txt) ]
    copy_umi  = FGUMI_COPYUMI.out.metrics    // channel: [ val(meta), path(copy_umi.tsv) ]
    align_log = BOWTIE_MAP_UMI_DEDUP.out.log // channel: [ val(meta), path(out) ]
}
