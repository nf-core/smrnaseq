//
// Alignment-based UMI deduplication with fgumi (`--umi_dedup_method fgumi`, experimental, #535).
//
// Input reads must already carry the UMI in the read name (`umi_tools extract`, `<name>_<UMI>`).
// Steps: Bowtie 1 (-k 1 --best: one record per read, required by fgumi) -> fgumi copy-umi (read name -> RX tag)
//        -> fgumi sort (template-coordinate) -> fgumi dedup (exact UMI matching, duplicates removed,
//        unaligned reads passed through) -> fgumi fastq.
// The deduplicated FASTQ keeps the `_<UMI>` read-name suffix, so all downstream steps
// (contamination filter, MIRNA_QUANT, GENOME_QUANT, miRDeep2) can use it like the UMICollapse output.
// The tool options are set in conf/modules.config.
//

include { BOWTIE_ALIGN as BOWTIE_MAP_UMI_DEDUP } from '../../../modules/nf-core/bowtie/align/main'
include { FGUMI_COPYUMI                        } from '../../../modules/nf-core/fgumi/copyumi/main'
include { FGUMI_SORT                           } from '../../../modules/nf-core/fgumi/sort/main'
include { FGUMI_DEDUP                          } from '../../../modules/nf-core/fgumi/dedup/main'
include { FGUMI_FASTQ                          } from '../../../modules/nf-core/fgumi/fastq/main'

workflow UMI_DEDUP_FGUMI {
    take:
    ch_reads        // channel: [ val(meta), path(reads) ]  UMI in read name
    ch_bowtie_index // channel: [ val(meta2), path(index) ]  index of the deduplication reference

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
    copy_umi  = FGUMI_COPYUMI.out.metrics    // channel: [ val(meta), path(metrics.txt) ]
    align_log = BOWTIE_MAP_UMI_DEDUP.out.log // channel: [ val(meta), path(out) ]
}
