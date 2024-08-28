//
// Filter contamination by rrna, trna, cdna, ncma, pirna
//

include { BLAT as BLAT_CDNA ; BLAT as BLAT_NCRNA; BLAT as BLAT_PIRNA ; BLAT as BLAT_OTHER} from '../../../modules/nf-core/blat/main'
include { GAWK as GAWK_CDNA ; GAWK as GAWK_NCRNA; GAWK as GAWK_PIRNA ; GAWK as GAWK_OTHER} from '../../../modules/nf-core/gawk/main'
include { SEQKIT_GREP as SEQKIT_GREP_CDNA ; SEQKIT_GREP as SEQKIT_GREP_NCRNA ; SEQKIT_GREP as SEQKIT_GREP_PIRNA ; SEQKIT_GREP as SEQKIT_GREP_OTHER} from '../../../modules/nf-core/seqkit/grep/main'

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
    mirna
    rrna
    trna
    cdna
    ncrna
    pirna
    other
    reads      // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()
    ch_filter_stats = Channel.empty()
    ch_mqc_results = Channel.empty()

    rrna_reads = reads

    reads.dump(tag:"reads")

    reads.set { rrna_reads }

    if (params.rrna) {
        // Index DB and filter $reads emit: $rrna_reads
        INDEX_RRNA ( rrna )
        ch_versions = ch_versions.mix(INDEX_RRNA.out.versions)
        MAP_RRNA ( reads, INDEX_RRNA.out.index, 'rRNA' )
        ch_versions = ch_versions.mix(MAP_RRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_RRNA.out.stats.ifEmpty(null))
        MAP_RRNA.out.unmapped.set { rrna_reads }
    }

    rrna_reads.set { trna_reads }

    if (params.trna) {
        // Index DB and filter $rrna_reads emit: $trna_reads
        INDEX_TRNA ( trna )
        ch_versions = ch_versions.mix(INDEX_TRNA.out.versions)
        MAP_TRNA ( rrna_reads, INDEX_TRNA.out.index, 'tRNA')
        ch_versions = ch_versions.mix(MAP_TRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_TRNA.out.stats.ifEmpty(null))
        MAP_TRNA.out.unmapped.set { trna_reads }
    }

    trna_reads.set { cdna_reads }


    if (params.cdna) {
        // Add metamap to input channels: cDNA and hairpin
        ch_cdna = Channel.value([[id:'cDNA'], cdna])
        ch_mirna = Channel.value([[id:'hairpin'], mirna])

        // Search which hairpin miRNAs are present in the cDNA data
        BLAT_CDNA(ch_mirna, ch_cdna)
        ch_versions = ch_versions.mix(BLAT_CDNA.out.versions)
        BLAT_CDNA.out.psl.dump(tag:"BLAT_CDNA")

        // Extract the significant hits
        ch_program = Channel.of('BEGIN{FS="\t"}{if(\$11 < 1e-5) print \$2;}').collectFile(name:"program.txt")
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
        // Index filtered cDNA
        INDEX_CDNA ( ch_filtered_cdna )
        ch_versions = ch_versions.mix(INDEX_CDNA.out.versions)

        // Map which input reads are not in the cDNA contaminants
        MAP_CDNA ( trna_reads, INDEX_CDNA.out.index, 'cDNA' )
        ch_versions = ch_versions.mix(MAP_CDNA.out.versions)

        // Extract number of reads aligning to contaminants
        ch_filter_stats = ch_filter_stats.mix(MAP_CDNA.out.stats.ifEmpty(null))

        // Create a channel with the set of input reads without cDNA
        MAP_CDNA.out.unmapped.set { cdna_reads }
    }

    cdna_reads.set { ncrna_reads }

    if (params.ncrna) {

        // Add metamap to input channels: ncRNA and hairpin
        ch_ncrna = Channel.value([[id:'ncRNA'], ncrna])
        ch_mirna = Channel.value([[id:'hairpin'], mirna])

        // Search which hairpin miRNAs are present in the ncRNA data
        BLAT_NCRNA(ch_mirna, ch_ncrna)
        BLAT_NCRNA.out.psl.dump(tag:"BLAT_NCRNA")

        // Extract the significant hits
        ch_program = Channel.of('BEGIN{FS="\t"}{if(\$11 < 1e-5) print \$2;}').collectFile(name:"program.txt")
        GAWK_NCRNA(BLAT_NCRNA.out.psl, ch_program)
        ch_versions = ch_versions.mix(GAWK_NCRNA.out.versions)
        GAWK_NCRNA.out.output.dump(tag:"GAWK_NCRNA")

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
        ch_filtered_ncrna.dump(tag:"ch_filtered_ncrna")

        // Previous original code:
        // Index filtered ncRNA
        INDEX_NCRNA ( ch_filtered_ncrna )
        ch_versions = ch_versions.mix(INDEX_NCRNA.out.versions)
        MAP_NCRNA ( cdna_reads, INDEX_NCRNA.out.index, 'ncRNA' )
        ch_versions = ch_versions.mix(MAP_NCRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_NCRNA.out.stats.ifEmpty(null))
        MAP_NCRNA.out.unmapped.set { ncrna_reads }
    }

    ncrna_reads.set { pirna_reads }

    if (params.pirna) {
        // Add metamap to input channels: piRNA and hairpin
        ch_pirna = Channel.value([[id:'piRNA'], pirna])
        ch_mirna = Channel.value([[id:'hairpin'], mirna])

        // Search which hairpin miRNAs are present in the piRNA data
        BLAT_PIRNA(ch_mirna, ch_pirna)
        BLAT_PIRNA.out.psl.dump(tag:"BLAT_PIRNA")

        // Extract the significant hits
        ch_program = Channel.of('BEGIN{FS="\t"}{if(\$11 < 1e-5) print \$2;}').collectFile(name:"program.txt")
        GAWK_PIRNA(BLAT_PIRNA.out.psl, ch_program)
        ch_versions = ch_versions.mix(GAWK_PIRNA.out.versions)
        GAWK_PIRNA.out.output.dump(tag:"GAWK_PIRNA")

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
        ch_filtered_pirna.dump(tag:"ch_filtered_pirna")

        // Previous original code:
        // Index filtered piRNA
        INDEX_PIRNA ( ch_filtered_pirna )
        ch_versions = ch_versions.mix(INDEX_PIRNA.out.versions)
        MAP_PIRNA ( ncrna_reads, INDEX_PIRNA.out.index, 'piRNA' )
        ch_versions = ch_versions.mix(MAP_PIRNA.out.versions)
        ch_filter_stats = ch_filter_stats.mix(MAP_PIRNA.out.stats.ifEmpty(null))
        MAP_PIRNA.out.unmapped.set { pirna_reads }
    }

    pirna_reads.set { other_cont_reads }

    if (other) {
        // Add metamap to input channels: other and hairpin
        ch_other = Channel.value([[id:'other'], other])
        ch_mirna = Channel.value([[id:'hairpin'], mirna])

        // Search which hairpin miRNAs are present in the other data
        BLAT_OTHER(ch_mirna, ch_other)
        BLAT_OTHER.out.psl.dump(tag:"BLAT_OTHER")

        // Extract the significant hits
        ch_program = Channel.of('BEGIN{FS="\t"}{if(\$11 < 1e-5) print \$2;}').collectFile(name:"program.txt")
        GAWK_OTHER(BLAT_OTHER.out.psl, ch_program)
        ch_versions = ch_versions.mix(GAWK_OTHER.out.versions)
        GAWK_OTHER.out.output.dump(tag:"GAWK_OTHER")

        // Get only unique elements of the list
        ch_pattern = GAWK_OTHER.out.output
                .map { meta, file -> file.text.readLines() }
                .flatten()
                .unique()
                .collectFile(name: 'ch_hairpin_other_unique.txt', newLine: true)

        // Remove the hairpin miRNAs from the other data
        SEQKIT_GREP_OTHER(ch_other, ch_pattern)
        ch_versions = ch_versions.mix(SEQKIT_GREP_OTHER.out.versions)

        // Remove metamap to make it compatible with previous code
        ch_filtered_other = SEQKIT_GREP_OTHER.out.filter.map{meta, file -> [file]}
        ch_filtered_other.dump(tag:"ch_filtered_other")

        // Previous original code:
        INDEX_OTHER ( ch_filtered_other )
        ch_versions = ch_versions.mix(INDEX_OTHER.out.versions)
        MAP_OTHER ( ncrna_reads, INDEX_OTHER.out.index, 'other' )
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
