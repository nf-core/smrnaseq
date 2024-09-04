//
// Filter contamination by rrna, trna, cdna, ncma, pirna
//

include { BLAT as BLAT_CDNA                } from '../../../modules/nf-core/blat/main'
include { BLAT as BLAT_NCRNA               } from '../../../modules/nf-core/blat/main'
include { BLAT as BLAT_PIRNA               } from '../../../modules/nf-core/blat/main'
include { BLAT as BLAT_OTHER               } from '../../../modules/nf-core/blat/main'

include { GAWK as GAWK_CDNA                } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_NCRNA               } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_PIRNA               } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_OTHER               } from '../../../modules/nf-core/gawk/main'

include { SEQKIT_GREP as SEQKIT_GREP_CDNA  } from '../../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_NCRNA } from '../../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_PIRNA } from '../../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_OTHER } from '../../../modules/nf-core/seqkit/grep/main'

include { INDEX_CONTAMINANTS as INDEX_RRNA
        INDEX_CONTAMINANTS as INDEX_TRNA
        INDEX_CONTAMINANTS as INDEX_CDNA
        INDEX_CONTAMINANTS as INDEX_NCRNA
        INDEX_CONTAMINANTS as INDEX_PIRNA
        INDEX_CONTAMINANTS as INDEX_OTHER } from '../../../modules/local/bowtie_contaminants'

include { BOWTIE_MAP_CONTAMINANTS as MAP_RRNA
        BOWTIE_MAP_CONTAMINANTS as MAP_TRNA
        BOWTIE_MAP_CONTAMINANTS as MAP_CDNA
        BOWTIE_MAP_CONTAMINANTS as MAP_NCRNA
        BOWTIE_MAP_CONTAMINANTS as MAP_PIRNA
        BOWTIE_MAP_CONTAMINANTS as MAP_OTHER } from '../../../modules/local/bowtie_map_contaminants'

include { FILTER_STATS } from '../../../modules/local/filter_stats'

workflow CONTAMINANT_FILTER {
    take:
    ch_reference_hairpin   // channel: [ val(meta), path(fasta) ]
    ch_rrna                // channel: [ path(fasta) ]
    ch_trna                // channel: [ path(fasta) ]
    ch_cdna                // channel: [ val(meta), path(fasta) ]
    ch_ncrna               // channel: [ val(meta), path(fasta) ]
    ch_pirna               // channel: [ val(meta), path(fasta) ]
    ch_other_contamination // channel: [ val(meta), path(fasta) ]
    ch_reads_for_mirna     // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions     = Channel.empty()
    ch_filter_stats = Channel.empty()
    ch_mqc_results  = Channel.empty()

    ch_reads_for_mirna.set { rrna_reads }

    if (params.rrna) {
        // Index DB and filter $reads emit: $rrna_reads
        INDEX_RRNA ( ch_rrna )
        ch_versions = ch_versions.mix(INDEX_RRNA.out.versions)
        MAP_RRNA ( ch_reads_for_mirna, INDEX_RRNA.out.index.first(), Channel.value('rRNA') )
        ch_versions = ch_versions.mix(MAP_RRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_RRNA.out.stats.ifEmpty(null))
        MAP_RRNA.out.unmapped.set { rrna_reads }
    }

    rrna_reads.set { trna_reads }

    if (params.trna) {
        // Index DB and filter $rrna_reads emit: $trna_reads
        INDEX_TRNA ( ch_trna )
        ch_versions = ch_versions.mix(INDEX_TRNA.out.versions)
        MAP_TRNA ( rrna_reads, INDEX_TRNA.out.index.first(), Channel.value("tRNA") )
        ch_versions = ch_versions.mix(MAP_TRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_TRNA.out.stats.ifEmpty(null))
        MAP_TRNA.out.unmapped.set { trna_reads }
    }

    trna_reads.set { cdna_reads }

    // Define how to filter significant BLAT hits
    ch_program = Channel.value('BEGIN{FS="\t"}{if(\$11 < 1e-5) print \$2;}').collectFile(name:"program.txt")

    if (params.cdna) {
        // Search which hairpin miRNAs are present in the cDNA data
        BLAT_CDNA(ch_reference_hairpin, ch_cdna)
        ch_versions = ch_versions.mix(BLAT_CDNA.out.versions)

        // Extract the significant hits
        GAWK_CDNA(BLAT_CDNA.out.psl, ch_program)
        ch_versions = ch_versions.mix(GAWK_CDNA.out.versions)

        // Get only unique elements of the list
        ch_pattern = GAWK_CDNA.out.output
                .map { meta, file -> file.text.readLines() }
                .flatten()
                .unique()
                .collectFile(name: 'ch_hairpin_cDNA_unique.txt', newLine: true)

        // Remove the hairpin miRNAs from the cDNA data
        SEQKIT_GREP_CDNA(ch_cdna, ch_pattern)
        ch_versions = ch_versions.mix(SEQKIT_GREP_CDNA.out.versions)

        // Remove metamap to make it compatible with previous code
        ch_filtered_cdna = SEQKIT_GREP_CDNA.out.filter.map{meta, file -> [file]}

        // Previous original code:
        INDEX_CDNA ( ch_filtered_cdna )
        ch_versions = ch_versions.mix(INDEX_CDNA.out.versions)
        MAP_CDNA ( trna_reads, INDEX_CDNA.out.index.first(), Channel.value('cDNA'))
        ch_versions = ch_versions.mix(MAP_CDNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_CDNA.out.stats.ifEmpty(null))
        MAP_CDNA.out.unmapped.set { cdna_reads }
    }

    cdna_reads.set { ncrna_reads }

    if (params.ncrna) {
        // Search which hairpin miRNAs are present in the ncRNA data
        BLAT_NCRNA(ch_reference_hairpin, ch_ncrna)
        ch_versions = ch_versions.mix(BLAT_NCRNA.out.versions)

        // Extract the significant hits
        GAWK_NCRNA(BLAT_NCRNA.out.psl, ch_program)
        ch_versions = ch_versions.mix(GAWK_NCRNA.out.versions)

        // Get only unique elements of the list
        ch_pattern = GAWK_NCRNA.out.output
                .map { meta, file -> file.text.readLines() }
                .flatten()
                .unique()
                .collectFile(name: 'ch_hairpin_ncRNA_unique.txt', newLine: true)

        // Remove the hairpin miRNAs from the ncRNA data
        SEQKIT_GREP_NCRNA(ch_ncrna, ch_pattern)
        ch_versions = ch_versions.mix(SEQKIT_GREP_NCRNA.out.versions)

        // Remove metamap to make it compatible with previous code
        ch_filtered_ncrna = SEQKIT_GREP_NCRNA.out.filter.map{meta, file -> [file]}

        // Previous original code:
        INDEX_NCRNA ( ch_filtered_ncrna )
        ch_versions = ch_versions.mix(INDEX_NCRNA.out.versions)
        MAP_NCRNA ( cdna_reads, INDEX_NCRNA.out.index.first(), Channel.value('ncRNA') )
        ch_versions = ch_versions.mix(MAP_NCRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_NCRNA.out.stats.ifEmpty(null))
        MAP_NCRNA.out.unmapped.set { ncrna_reads }
    }

    ncrna_reads.set { pirna_reads }

    if (params.pirna) {
        // Search which hairpin miRNAs are present in the piRNA data
        BLAT_PIRNA(ch_reference_hairpin, ch_pirna)
        ch_versions = ch_versions.mix(BLAT_PIRNA.out.versions)

        // Extract the significant hits
        GAWK_PIRNA(BLAT_PIRNA.out.psl, ch_program)
        ch_versions = ch_versions.mix(GAWK_PIRNA.out.versions)

        // Get only unique elements of the list
        ch_pattern = GAWK_PIRNA.out.output
                .map { meta, file -> file.text.readLines() }
                .flatten()
                .unique()
                .collectFile(name: 'ch_hairpin_piRNA_unique.txt', newLine: true)

        // Remove the hairpin miRNAs from the piRNA data
        SEQKIT_GREP_PIRNA(ch_pirna, ch_pattern)
        ch_versions = ch_versions.mix(SEQKIT_GREP_PIRNA.out.versions)

        // Remove metamap to make it compatible with previous code
        ch_filtered_pirna = SEQKIT_GREP_PIRNA.out.filter.map{meta, file -> [file]}

        // Previous original code:
        INDEX_PIRNA ( ch_filtered_pirna )
        ch_versions = ch_versions.mix(INDEX_PIRNA.out.versions)
        MAP_PIRNA ( ncrna_reads, INDEX_PIRNA.out.index.first(), Channel.value('piRNA'))
        ch_versions = ch_versions.mix(MAP_PIRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_PIRNA.out.stats.ifEmpty(null))
        MAP_PIRNA.out.unmapped.set { pirna_reads }
    }

    pirna_reads.set { other_cont_reads }

    if (params.other_contamination) {
        // Search which hairpin miRNAs are present in the other data
        BLAT_OTHER(ch_reference_hairpin, ch_other_contamination)
        ch_versions = ch_versions.mix(BLAT_OTHER.out.versions)

        // Extract the significant hits
        GAWK_OTHER(BLAT_OTHER.out.psl, ch_program)
        ch_versions = ch_versions.mix(GAWK_OTHER.out.versions)

        // Get only unique elements of the list
        ch_pattern = GAWK_OTHER.out.output
                .map { meta, file -> file.text.readLines() }
                .flatten()
                .unique()
                .collectFile(name: 'ch_hairpin_other_unique.txt', newLine: true)

        // Remove the hairpin miRNAs from the other data
        SEQKIT_GREP_OTHER(ch_other_contamination, ch_pattern)
        ch_versions = ch_versions.mix(SEQKIT_GREP_OTHER.out.versions)

        // Remove metamap to make it compatible with previous code
        ch_filtered_other = SEQKIT_GREP_OTHER.out.filter.map{meta, file -> [file]}

        // Previous original code:
        INDEX_OTHER ( ch_filtered_other )
        ch_versions = ch_versions.mix(INDEX_OTHER.out.versions)
        MAP_OTHER ( ncrna_reads, INDEX_OTHER.out.index.first(), Channel.value('other'))
        ch_versions = ch_versions.mix(MAP_OTHER.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_OTHER.out.stats.ifEmpty(null))
        MAP_OTHER.out.unmapped.set { other_cont_reads }
    }

    FILTER_STATS ( other_cont_reads, ch_filter_stats.collect() )

    emit:
    filtered_reads = FILTER_STATS.out.reads
    versions = ch_versions.mix(FILTER_STATS.out.versions)
    filter_stats = FILTER_STATS.out.stats
}
