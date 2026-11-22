//
// Build a cross-sample consensus BED of the miRNAs (known and novel) reported by miRDeep2
//

include { GAWK as GAWK_MIRDEEP2_SAMPLES     } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_SORT                     } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MERGE                    } from '../../../modules/nf-core/bedtools/merge/main'
include { GAWK as GAWK_MIRDEEP2_CONSENSUS   } from '../../../modules/nf-core/gawk/main'

workflow MIRDEEP2_CONSENSUS {
    take:
    ch_mirdeep2  // channel: [ val(meta), [ bed, csv, html ] ]
    min_samples  // integer: minimum number of samples supporting a locus

    main:
    // All per-sample miRDeep2 BED files in a single task. Samples where miRDeep2 failed
    // with exit status 255 (ignored) have no output and do not count.
    ch_beds = ch_mirdeep2
        .map { _meta, files -> [ files ].flatten().find { f -> f.name.endsWith('.bed') } }
        .collect(sort: true)
        .map { beds -> [ [ id: 'mirdeep2_consensus' ], beds, 'bed' ] }

    // Concatenate the BED files, prefix each miRDeep2 ID with its sample and add a sample column:
    // chrom, start, end, <sample>:<known|novel>:<id>, miRDeep2 score, strand, sample
    ch_samples_program = channel.value(
        'BEGIN { FS = OFS = "\t" } ' +
        'FNR == 1 { sample = FILENAME; sub(/^.*\\//, "", sample); sub(/^result_/, "", sample); sub(/\\.bed$/, "", sample) } ' +
        '/^(track|browser|#)/ || NF < 6 { next } ' +
        '{ print $1, $2, $3, sample ":" $4, $5, $6, sample }'
    ).collectFile(name: 'mirdeep2_samples.awk')

    GAWK_MIRDEEP2_SAMPLES ( ch_beds, ch_samples_program, false )

    BEDTOOLS_SORT ( GAWK_MIRDEEP2_SAMPLES.out.output, [], 'bed' )

    // Merge overlapping loci on the same strand (ext.args in conf/modules.config):
    // chrom, start, end, strand, IDs, max score, samples, number of samples
    BEDTOOLS_MERGE ( BEDTOOLS_SORT.out.sorted )

    // Keep loci found in at least min_samples samples and write BED6 plus support columns
    ch_consensus_program = channel.value(
        'BEGIN { FS = OFS = "\t"; print "#chrom", "start", "end", "name", "score", "strand", "n_samples", "samples", "mirdeep2_ids" } ' +
        '$8 >= ' + min_samples + ' { n++; type = ($5 ~ /:known:/) ? "known" : "novel"; ' +
        'print $1, $2, $3, type ":consensus_" n, $6, $4, $8, $7, $5 }'
    ).collectFile(name: 'mirdeep2_consensus.awk')

    GAWK_MIRDEEP2_CONSENSUS ( BEDTOOLS_MERGE.out.bed.map { meta, bed -> [ meta, bed, 'bed' ] }, ch_consensus_program, false )

    emit:
    bed = GAWK_MIRDEEP2_CONSENSUS.out.output // channel: [ val(meta), path(bed) ]
}
