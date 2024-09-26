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

include { BOWTIE2_BUILD as INDEX_TRNA      } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as INDEX_CDNA      } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as INDEX_NCRNA     } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as INDEX_PIRNA     } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as INDEX_OTHER     } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as INDEX_RRNA      } from '../../../modules/nf-core/bowtie2/build/main'

include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_RRNA  } from '../../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_TRNA  } from '../../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_CDNA  } from '../../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_NCRNA } from '../../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_PIRNA } from '../../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_OTHER } from '../../../modules/nf-core/bowtie2/align/main'

include { GAWK as STATS_GAWK_RRNA                } from '../../../modules/nf-core/gawk/main'
include { GAWK as STATS_GAWK_TRNA                } from '../../../modules/nf-core/gawk/main'
include { GAWK as STATS_GAWK_CDNA                } from '../../../modules/nf-core/gawk/main'
include { GAWK as STATS_GAWK_NCRNA               } from '../../../modules/nf-core/gawk/main'
include { GAWK as STATS_GAWK_PIRNA               } from '../../../modules/nf-core/gawk/main'
include { GAWK as STATS_GAWK_OTHER               } from '../../../modules/nf-core/gawk/main'


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

        // Add meta.contaminant to input reads channel
        ch_reads_for_mirna = ch_reads_for_mirna.map{meta, fastq -> return [[id:meta.id, contaminant: "rRNA", single_end:meta.single_end], fastq]}

        // Map which reads are rRNAs
        BOWTIE2_ALIGN_RRNA(ch_reads_for_mirna, INDEX_RRNA.out.index.first(), [[],[]], true, false)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN_RRNA.out.versions)

        // Obtain how many hits were contaminants
        ch_bowtie = BOWTIE2_ALIGN_RRNA.out.log

        STATS_GAWK_RRNA(ch_bowtie, [])
        ch_versions = ch_versions.mix(STATS_GAWK_RRNA.out.versions)

        // Remove meta.contaminant and collect all contaminant stats in a single channel
        ch_filter_stats = ch_filter_stats
                .mix(STATS_GAWK_RRNA.out.output
                        .map{meta, stats -> return [[id:meta.id, single_end:meta.single_end], stats]}
                        .ifEmpty(null))

        // Assign clean reads to new channel
        rrna_reads = BOWTIE2_ALIGN_RRNA.out.fastq
    }

    rrna_reads.set { trna_reads }

    if (params.trna) {
        // Index DB and filter $rrna_reads emit: $trna_reads
        INDEX_TRNA ( ch_trna )
        ch_versions = ch_versions.mix(INDEX_TRNA.out.versions)

        // Add meta.contaminant to input reads channel
        rrna_reads = rrna_reads.map{meta, fastq -> return [[id:meta.id, contaminant: "tRNA", single_end:meta.single_end], fastq]}

        // Map which reads are tRNAs
        BOWTIE2_ALIGN_TRNA(rrna_reads, INDEX_TRNA.out.index.first(), [[],[]], true, false)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN_TRNA.out.versions)

        // Obtain how many hits were contaminants
        ch_bowtie = BOWTIE2_ALIGN_TRNA.out.log

        STATS_GAWK_TRNA(ch_bowtie, [])
        ch_versions = ch_versions.mix(STATS_GAWK_TRNA.out.versions)

        // Remove meta.contaminant and collect all contaminant stats in a single channel
        ch_filter_stats = ch_filter_stats
                .mix(STATS_GAWK_TRNA.out.output
                        .map{meta, stats -> return [[id:meta.id, single_end:meta.single_end], stats]}
                        .ifEmpty(null))

        // Assign clean reads to new channel
        trna_reads = BOWTIE2_ALIGN_TRNA.out.fastq
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

        // Previous original code:
        INDEX_CDNA ( SEQKIT_GREP_CDNA.out.filter )
        ch_versions = ch_versions.mix(INDEX_CDNA.out.versions)

        // Add meta.contaminant to input reads channel
        trna_reads = trna_reads.map{meta, fastq -> return [[id:meta.id, contaminant: "cDNA", single_end:meta.single_end], fastq]}

        // Map which reads are cDNA
        BOWTIE2_ALIGN_CDNA(trna_reads, INDEX_CDNA.out.index.first(), [[],[]], true, false)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN_CDNA.out.versions)

        // Obtain how many hits were contaminants
        STATS_GAWK_CDNA(BOWTIE2_ALIGN_CDNA.out.log, [])
        ch_versions = ch_versions.mix(STATS_GAWK_CDNA.out.versions)

        // Remove meta.contaminant and collect all contaminant stats in a single channel
        ch_filter_stats = ch_filter_stats
                .mix(STATS_GAWK_CDNA.out.output
                        .map{meta, stats -> return [[id:meta.id, single_end:meta.single_end], stats]}
                        .ifEmpty(null))

        // Assign clean reads to new channel
        cdna_reads = BOWTIE2_ALIGN_CDNA.out.fastq
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

        // Previous original code:
        INDEX_NCRNA ( SEQKIT_GREP_NCRNA.out.filter )
        ch_versions = ch_versions.mix(INDEX_NCRNA.out.versions)

        // Add meta.contaminant to input reads channel
        cdna_reads = cdna_reads.map{meta, fastq -> return [[id:meta.id, contaminant: "ncRNA", single_end:meta.single_end], fastq]}

        // Map which reads are ncRNA
        BOWTIE2_ALIGN_NCRNA(cdna_reads, INDEX_NCRNA.out.index.first(), [[],[]], true, false)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN_NCRNA.out.versions)

        // Obtain how many hits were contaminants
        STATS_GAWK_NCRNA(BOWTIE2_ALIGN_NCRNA.out.log, [])
        ch_versions = ch_versions.mix(STATS_GAWK_NCRNA.out.versions)

        // Remove meta.contaminant and collect all contaminant stats in a single channel
        ch_filter_stats = ch_filter_stats
                .mix(STATS_GAWK_NCRNA.out.output
                        .map{meta, stats -> return [[id:meta.id, single_end:meta.single_end], stats]}
                        .ifEmpty(null))

        // Assign clean reads to new channel
        ncrna_reads = BOWTIE2_ALIGN_NCRNA.out.fastq
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

        // Previous original code:
        INDEX_PIRNA ( SEQKIT_GREP_PIRNA.out.filter )
        ch_versions = ch_versions.mix(INDEX_PIRNA.out.versions)

        // Add meta.contaminant to input reads channel
        ncrna_reads = ncrna_reads.map{meta, fastq -> return [[id:meta.id, contaminant: "piRNA", single_end:meta.single_end], fastq]}

        // Map which reads are piRNA
        BOWTIE2_ALIGN_PIRNA(ncrna_reads, INDEX_PIRNA.out.index.first(), [[],[]], true, false)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN_PIRNA.out.versions)

        // Obtain how many hits were contaminants
        STATS_GAWK_PIRNA(BOWTIE2_ALIGN_PIRNA.out.log, [])
        ch_versions = ch_versions.mix(STATS_GAWK_PIRNA.out.versions)

        // Remove meta.contaminant and collect all contaminant stats in a single channel
        ch_filter_stats = ch_filter_stats
                .mix(STATS_GAWK_PIRNA.out.output
                        .map{meta, stats -> return [[id:meta.id, single_end:meta.single_end], stats]}
                        .ifEmpty(null))

        // Assign clean reads to new channel
        pirna_reads = BOWTIE2_ALIGN_PIRNA.out.fastq
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

        // Previous original code:
        INDEX_OTHER ( SEQKIT_GREP_OTHER.out.filter )
        ch_versions = ch_versions.mix(INDEX_OTHER.out.versions)

        // Map which reads are other
        BOWTIE2_ALIGN_OTHER(pirna_reads, INDEX_OTHER.out.index.first(), [[],[]], true, false)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN_OTHER.out.versions)

        // Obtain how many hits were contaminants
        STATS_GAWK_OTHER(BOWTIE2_ALIGN_OTHER.out.log, [])
        ch_versions = ch_versions.mix(STATS_GAWK_OTHER.out.versions)

        // Remove meta.contaminant and collect all contaminant stats in a single channel
        ch_filter_stats = ch_filter_stats
                .mix(STATS_GAWK_OTHER.out.output
                        .map{meta, stats -> return [[id:meta.id, single_end:meta.single_end], stats]}
                        .ifEmpty(null))

        // Assign clean reads to new channel
        other_cont_reads = BOWTIE2_ALIGN_OTHER.out.fastq
    }

    // Remove meta.contaminant from final set of reads
    other_cont_reads = other_cont_reads
                        .map{meta, reads -> return [[id:meta.id, single_end:meta.single_end], reads]}

    // Create channel with reads and contaminants
    ch_reads_contaminants = other_cont_reads.join(ch_filter_stats.groupTuple())

    // Filter all contaminant stats and create MultiQC file
    FILTER_STATS ( ch_reads_contaminants )
    FILTER_STATS.out.stats.dump(tag:"FILTER_STATS.out.stats")

    emit:
    filtered_reads  = other_cont_reads          // channel: [ val(meta), path(fastq) ]
    filter_stats    = FILTER_STATS.out.stats    // channel: [  path(stats) ]
    versions        = ch_versions.mix(FILTER_STATS.out.versions)
}
