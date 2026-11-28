//
// Build a cross-sample consensus (BED6 and TSV) of the miRNAs (known and novel) reported by miRDeep2
//

include { GAWK as GAWK_MIRDEEP2_SAMPLES     } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_SORT                     } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MERGE                    } from '../../../modules/nf-core/bedtools/merge/main'
include { GAWK as GAWK_MIRDEEP2_CONSENSUS_TSV } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_MIRDEEP2_CONSENSUS_BED } from '../../../modules/nf-core/gawk/main'

workflow MIRDEEP2_CONSENSUS {
    take:
    ch_mirdeep2      // channel: [ val(meta), [ bed, csv, html ] ]
    min_samples      // integer: minimum number of samples supporting a locus
    min_probability  // integer: minimum miRDeep2 true-positive probability in percent (0: no filter)

    main:
    // All per-sample miRDeep2 BED and CSV files in a single task, CSV files first so that
    // their probabilities are known when the BED files are read. Samples where miRDeep2
    // failed with exit status 255 (ignored) have no output and do not count.
    ch_files = ch_mirdeep2
        .flatMap { _meta, files -> [ files ].flatten().findAll { f -> f.name.endsWith('.bed') || f.name.endsWith('.csv') } }
        .collect(sort: true)
        .map { files ->
            def csvs = files.findAll { f -> f.name.endsWith('.csv') }
            def beds = files.findAll { f -> f.name.endsWith('.bed') }
            [ [ id: 'mirdeep2_consensus' ], csvs + beds, 'bed' ]
        }

    // Read the miRDeep2 true-positive probability of each miRNA from the CSV files: the
    // third column of the novel and the known ("mature miRBase miRNAs detected") sections,
    // e.g. "91 +/- 2%", or "-" when miRDeep2 gives no estimate. The BED name field is
    // <novel|known>:<first CSV column> (see mirdeep2bed.pl).
    // Then concatenate the BED files, drop entries below min_probability (with min_probability
    // above 0, entries without a probability are dropped too), prefix each miRDeep2 ID with
    // its sample and add columns:
    // chrom, start, end, <sample>:<known|novel>:<id>, miRDeep2 score, strand, sample, probability
    ch_samples_program = channel.value(
        'BEGIN { FS = OFS = "\t"; min_probability = ' + min_probability + ' } ' +
        'FNR == 1 { sample = FILENAME; sub(/^.*\\//, "", sample); sub(/^result_/, "", sample); ' +
        'is_csv = sub(/\\.csv$/, "", sample); sub(/\\.bed$/, "", sample); section = "" } ' +
        'is_csv { ' +
        'if ($0 ~ /novel miRNAs predicted by miRDeep2/) section = "novel"; ' +
        'else if ($0 ~ /mature miRBase miRNAs detected/) section = "known"; ' +
        'else if ($0 ~ /miRBase miRNAs not detected/) section = ""; ' +
        'else if (section != "" && NF >= 3 && $1 != "provisional id" && $1 != "tag id") ' +
        'prob[sample, section ":" $1] = match($3, /^[0-9]+(\\.[0-9]+)?/) ? substr($3, 1, RLENGTH) : "-"; ' +
        'next } ' +
        '/^(track|browser|#)/ || NF < 6 { next } ' +
        '{ key = sample SUBSEP $4; p = (key in prob) ? prob[key] : "-"; ' +
        'if (min_probability > 0 && (p == "-" || p + 0 < min_probability)) next; ' +
        'print $1, $2, $3, sample ":" $4, $5, $6, sample, p }'
    ).collectFile(name: 'mirdeep2_samples.awk')

    GAWK_MIRDEEP2_SAMPLES ( ch_files, ch_samples_program, false )

    BEDTOOLS_SORT ( GAWK_MIRDEEP2_SAMPLES.out.output, [], 'bed' )

    // Merge loci that overlap by at least 1 bp on the same strand (ext.args in conf/modules.config;
    // -d -1 keeps book-ended loci, where one ends at the start of the next, separate):
    // chrom, start, end, strand, IDs, max score, samples, number of samples, probabilities
    BEDTOOLS_MERGE ( BEDTOOLS_SORT.out.sorted )

    // Keep loci found in at least min_samples samples and write one row per locus with its
    // support. The probability is the highest one across the samples, "-" if none has one.
    ch_tsv_program = channel.value(
        'BEGIN { FS = OFS = "\t"; print "name", "chrom", "start", "end", "strand", "n_samples", "samples", "score", "probability", "mirdeep2_ids" } ' +
        '$8 >= ' + min_samples + ' { n++; type = ($5 ~ /:known:/) ? "known" : "novel"; ' +
        'best = "-"; k = split($9, ps, ","); ' +
        'for (i = 1; i <= k; i++) if (ps[i] != "-" && (best == "-" || ps[i] + 0 > best + 0)) best = ps[i]; ' +
        'print type ":consensus_" n, $1, $2, $3, $4, $8, $7, $6, best, $5 }'
    ).collectFile(name: 'mirdeep2_consensus_tsv.awk')

    GAWK_MIRDEEP2_CONSENSUS_TSV ( BEDTOOLS_MERGE.out.bed.map { meta, bed -> [ meta, bed, 'tsv' ] }, ch_tsv_program, false )

    // Standard BED6 of the same loci (no header, no extra columns), so that genome browsers
    // read it correctly. The BED score has to be an integer from 0 to 1000: it is the number
    // of supporting samples. The miRDeep2 score and all other support are in the TSV.
    ch_bed_program = channel.value(
        'BEGIN { FS = OFS = "\t" } ' +
        'FNR > 1 { print $2, $3, $4, $1, ($6 > 1000 ? 1000 : $6), $5 }'
    ).collectFile(name: 'mirdeep2_consensus_bed.awk')

    GAWK_MIRDEEP2_CONSENSUS_BED ( GAWK_MIRDEEP2_CONSENSUS_TSV.out.output.map { meta, tsv -> [ meta, tsv, 'bed' ] }, ch_bed_program, false )

    emit:
    bed = GAWK_MIRDEEP2_CONSENSUS_BED.out.output // channel: [ val(meta), path(bed) ]
    tsv = GAWK_MIRDEEP2_CONSENSUS_TSV.out.output // channel: [ val(meta), path(tsv) ]
}
