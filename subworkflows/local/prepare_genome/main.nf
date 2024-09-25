//
// Uncompress and prepare reference genome files
//

// nf-core modules
include { UNTARFILES    as UNTAR_BOWTIE_INDEX } from '../../../modules/nf-core/untarfiles'
include { BOWTIE_BUILD  as INDEX_GENOME       } from '../../../modules/nf-core/bowtie/build/main'
include { BIOAWK        as CLEAN_FASTA        } from '../../../modules/nf-core/bioawk/main'


workflow PREPARE_GENOME {
    take:
    val_fasta                      // file: /path/to/genome.fasta
    val_bowtie_index               // file or directory: /path/to/bowtie/ or /path/to/bowtie.tar.gz
    val_mirtrace_species           // string: Species for miRTrace
    val_rrna                       // string: Path to the rRNA fasta file to be used as contamination database.
    val_trna                       // string: Path to the tRNA fasta file to be used as contamination database.
    val_cdna                       // string: Path to the cDNA fasta file to be used as contamination database.
    val_ncrna                      // string: Path to the ncRNA fasta file to be used as contamination database.
    val_pirna                      // string: Path to the piRNA fasta file to be used as contamination database.
    val_other_contamination        // string: Path to the additional fasta file to be used as contamination database.
    val_fastp_known_mirna_adapters // string: Path to Fasta with known miRNA adapter sequences for adapter trimming
    val_mirna_gtf                  // string: Path to GFF/GTF file with coordinates positions of precursor and miRNAs

    main:
    ch_versions = Channel.empty()

    // Parameter channel handling
    ch_fasta                  = val_fasta                     ? Channel.fromPath(val_fasta, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect()        : Channel.empty()
    ch_bowtie_index           = val_bowtie_index              ? Channel.fromPath(val_bowtie_index, checkIfExists: true).map{ it -> [ [], it ] }.collect() : Channel.empty()

    bool_mirtrace_species     = val_mirtrace_species          ? true : false
    bool_has_fasta            = val_fasta                     ? true : false

    ch_mirtrace_species       = val_mirtrace_species          ? Channel.value(val_mirtrace_species) : Channel.empty()
    mirna_gtf_from_species    = val_mirtrace_species          ? (val_mirtrace_species == 'hsa' ? "https://raw.githubusercontent.com/nf-core/test-datasets/smrnaseq/reference/hsa.gff3" : "https://mirbase.org/download/CURRENT/genomes/${val_mirtrace_species}.gff3") : false
    ch_mirna_gtf              = val_mirna_gtf                 ? Channel.fromPath(val_mirna_gtf, checkIfExists: true) : ( mirna_gtf_from_species ? Channel.fromPath(mirna_gtf_from_species, checkIfExists: true).collect() :  Channel.empty() )
    ch_mirna_adapters         = params.with_umi               ? [] : Channel.fromPath(val_fastp_known_mirna_adapters, checkIfExists: true).collect()

    ch_rrna                   = val_rrna                      ? Channel.fromPath(val_rrna).map{ it -> [ [id:'rRNA'], it ] }                : Channel.empty()
    ch_trna                   = val_trna                      ? Channel.fromPath(val_trna).map{ it -> [ [id:'tRNA'], it ] }.collect()                   : Channel.empty()
    ch_cdna                   = val_cdna                      ? Channel.fromPath(val_cdna).map{ it -> [ [id:'cDNA'], it ] }.collect()                : Channel.empty()
    ch_ncrna                  = val_ncrna                     ? Channel.fromPath(val_ncrna).map{ it -> [ [id:'ncRNA'], it ] }.collect()               : Channel.empty()
    ch_pirna                  = val_pirna                     ? Channel.fromPath(val_pirna).map{ it -> [ [id:'piRNA'], it ] }.collect()               : Channel.empty()
    ch_other_contamination    = val_other_contamination       ? Channel.fromPath(val_other_contamination).map{ it -> [ [id:'other'], it ] }.collect() : Channel.empty()

    // even if bowtie index is specified, there still needs to be a fasta.
    // without fasta, no genome analysis.
    if(val_fasta) {
        //Prepare bowtie index, unless specified
        //This needs to be done here as the index is used by GENOME_QUANT
        if(val_bowtie_index) {
            if (val_bowtie_index.endsWith(".tar.gz")) {
                UNTAR_BOWTIE_INDEX ( ch_bowtie_index )
                ch_bowtie_index = UNTAR_BOWTIE_INDEX.out.files
                ch_versions  = ch_versions.mix(UNTAR_BOWTIE_INDEX.out.versions)
            } else {
                ch_bowtie_index = Channel.fromPath("${val_bowtie_index}**ebwt", checkIfExists: true).map{it -> [ [id:it.baseName], it ] }.collect()
                    .ifEmpty{ error "Bowtie1 index directory not found: ${val_bowtie_index}" }
                    .filter { it != null }
            }
        } else {
            // Clean fasta (replace non-ATCGs with Ns, remove whitespaces from headers)
            CLEAN_FASTA ( ch_fasta )
            ch_versions      = ch_versions.mix(CLEAN_FASTA.out.versions)

            // Index FASTA with nf-core Bowtie1
            INDEX_GENOME ( CLEAN_FASTA.out.output )
            ch_versions      = ch_versions.mix(INDEX_GENOME.out.versions)

            // Set channels: clean fasta and its index
            ch_fasta         = CLEAN_FASTA.out.output
            ch_bowtie_index  = INDEX_GENOME.out.index.collect()
        }
    }

    //Config checks
    // Check optional parameters
    if (!params.mirgenedb && !val_mirtrace_species) {
            exit 1, "Reference species for miRTrace is not defined via the --mirtrace_species parameter."
        }

    // Genome options
    if (!params.mirgenedb) {
        ch_reference_mature  = params.mature  ? Channel.fromPath(params.mature, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect()  : { exit 1, "Mature miRNA fasta file not found: ${params.mature}" }
        ch_reference_hairpin = params.hairpin ? Channel.fromPath(params.hairpin, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() : { exit 1, "Hairpin miRNA fasta file not found: ${params.hairpin}" }
    } else {
        if (!params.mirgenedb_species) {
            exit 1, "MirGeneDB species not set, please specify via the --mirgenedb_species parameter"
        }
        ch_reference_mature  = params.mirgenedb_mature  ? Channel.fromPath(params.mirgenedb_mature, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect()  : { exit 1, "Mature miRNA fasta file not found via --mirgenedb_mature: ${params.mirgenedb_mature}" }
        ch_reference_hairpin = params.mirgenedb_hairpin ? Channel.fromPath(params.mirgenedb_hairpin, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() : { exit 1, "Hairpin miRNA fasta file not found via --mirgenedb_hairpin: ${params.mirgenedb_hairpin}" }
        ch_mirna_gtf         = params.mirgenedb_gff     ? Channel.fromPath(params.mirgenedb_gff, checkIfExists: true).collect()                                           : { exit 1, "MirGeneDB gff file not found via --mirgenedb_gff: ${params.mirgenedb_gff}"}
    }

    emit:
    fasta                 = ch_fasta               // channel: [ val(meta), path(fasta) ]
    has_fasta             = bool_has_fasta         // boolean
    bowtie_index          = ch_bowtie_index        // channel: [genome.1.ebwt, genome.2.ebwt, genome.3.ebwt, genome.4.ebwt, genome.rev.1.ebwt, genome.rev.2.ebwt]
    versions              = ch_versions            // channel: [ versions.yml ]
    mirtrace_species      = ch_mirtrace_species    // channel: [ val(string) ]
    has_mirtrace_species  = bool_mirtrace_species  // boolean
    reference_mature      = ch_reference_mature    // channel: [ val(meta), path(fasta) ]
    reference_hairpin     = ch_reference_hairpin   // channel: [ val(meta), path(fasta) ]
    mirna_gtf             = ch_mirna_gtf           // channel: [ path(GTF) ]
    rrna                  = ch_rrna                // channel: [ path(fasta) ]
    trna                  = ch_trna                // channel: [ path(fasta) ]
    cdna                  = ch_cdna                // channel: [ path(fasta) ]
    ncrna                 = ch_ncrna               // channel: [ path(fasta) ]
    pirna                 = ch_pirna               // channel: [ path(fasta) ]
    other_contamination   = ch_other_contamination // channel: [ path(fasta) ]
    mirna_adapters        = ch_mirna_adapters      // channel: [ path(fasta) ]
}
