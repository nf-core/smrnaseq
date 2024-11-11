//
// Uncompress and prepare reference genome files
//

// nf-core modules
include { UNTAR         as UNTAR_BOWTIE_INDEX } from '../../../modules/nf-core/untar'
include { BOWTIE_BUILD  as INDEX_GENOME       } from '../../../modules/nf-core/bowtie/build/main'
include { BIOAWK        as CLEAN_FASTA        } from '../../../modules/nf-core/bioawk/main'

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Extract prefix from bowtie index files
//
def extractFirstIndexPrefix(files_path) {
    def files = files_path.listFiles()
    if (files == null || files.length == 0) {
        throw new Exception("The provided bowtie_index path doesn't contain any files.")
    }
    def index_prefix = ''
    for (file_path in files) {
        def file_name = file_path.getName()
        if (file_name.endsWith(".1.ebwt") && !file_name.endsWith(".rev.1.ebwt")) {
            index_prefix = file_name.substring(0, file_name.lastIndexOf(".1.ebwt"))
            break
        }
    }
    if (index_prefix == '') {
        throw new Exception("Unable to extract the prefix from the Bowtie index files. No file with the '.1.ebwt' extension was found. Please ensure that the correct files are in the specified path.")
    }
    return index_prefix
}


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
    ch_fasta                  = val_fasta                     ? Channel.fromPath(val_fasta, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect()  : Channel.empty()
    ch_bowtie_index           = val_bowtie_index              ? Channel.fromPath(val_bowtie_index, checkIfExists: true).map{ it -> [ [id:""], it ] }.collect()    : Channel.empty()

    bool_mirtrace_species     = val_mirtrace_species          ? true : false
    bool_has_fasta            = val_fasta                     ? true : false

    ch_mirtrace_species       = val_mirtrace_species          ? Channel.value(val_mirtrace_species) : Channel.empty()
    mirna_gtf_from_species    = val_mirtrace_species          ? (val_mirtrace_species == 'hsa' ? "https://raw.githubusercontent.com/nf-core/test-datasets/smrnaseq/reference/hsa.gff3" : "https://mirbase.org/download/CURRENT/genomes/${val_mirtrace_species}.gff3") : false
    ch_mirna_gtf              = val_mirna_gtf                 ? Channel.fromPath(val_mirna_gtf, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() : ( mirna_gtf_from_species ? Channel.fromPath(mirna_gtf_from_species, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() :  Channel.empty() )

    // Add species of the gtf in the meta
    ch_mirna_gtf = ch_mirna_gtf
        .combine(ch_mirtrace_species.ifEmpty('unknown'))
        .map { meta, gtf, species ->
            def new_meta = meta.clone() + [species: species]
            [new_meta, gtf]
        }

    ch_mirna_adapters         = params.with_umi               ? [] : Channel.fromPath(val_fastp_known_mirna_adapters, checkIfExists: true).collect()

    ch_rrna                   = val_rrna                      ? Channel.fromPath(val_rrna, checkIfExists: true).map{ it -> [ [id:'rRNA'], it ] }.collect()                 : Channel.empty()
    ch_trna                   = val_trna                      ? Channel.fromPath(val_trna, checkIfExists: true).map{ it -> [ [id:'tRNA'], it ] }.collect()                 : Channel.empty()
    ch_cdna                   = val_cdna                      ? Channel.fromPath(val_cdna, checkIfExists: true).map{ it -> [ [id:'cDNA'], it ] }.collect()                 : Channel.empty()
    ch_ncrna                  = val_ncrna                     ? Channel.fromPath(val_ncrna, checkIfExists: true).map{ it -> [ [id:'ncRNA'], it ] }.collect()               : Channel.empty()
    ch_pirna                  = val_pirna                     ? Channel.fromPath(val_pirna, checkIfExists: true).map{ it -> [ [id:'piRNA'], it ] }.collect()               : Channel.empty()
    ch_other_contamination    = val_other_contamination       ? Channel.fromPath(val_other_contamination, checkIfExists: true).map{ it -> [ [id:'other'], it ] }.collect() : Channel.empty()

    // even if bowtie index is specified, there still needs to be a fasta.
    // without fasta, no genome analysis.
    if(val_fasta) {
        // Clean fasta (replace non-ATCGs with Ns, remove whitespaces from headers)
        // Note: CLEAN_FASTA runs even when a bowtie_index is provided, as cleaning doesn't affect it, making regeneration unnecessary.
        CLEAN_FASTA ( ch_fasta )
        ch_versions      = ch_versions.mix(CLEAN_FASTA.out.versions)
        ch_fasta         = CLEAN_FASTA.out.output

        //Prepare bowtie index, unless specified
        //This needs to be done here as the index is used by GENOME_QUANT
        if(val_bowtie_index) {
            if (val_bowtie_index.endsWith(".tar.gz")) {
                UNTAR_BOWTIE_INDEX ( ch_bowtie_index )
                ch_bowtie_index = UNTAR_BOWTIE_INDEX.out.untar
                    .map{ meta, index_dir ->
                        def index_prefix = extractFirstIndexPrefix(index_dir)
                        [[id:index_prefix], index_dir]
                    }
                ch_versions  = ch_versions.mix(UNTAR_BOWTIE_INDEX.out.versions)
            } else {
                ch_bowtie_index = Channel.fromPath(val_bowtie_index, checkIfExists: true)
                    .map{it ->
                        def index_prefix = extractFirstIndexPrefix(it)
                        [[id:index_prefix], it]
                    }
            }

        } else {

            // Index FASTA with nf-core Bowtie1
            INDEX_GENOME ( CLEAN_FASTA.out.output )
            ch_versions      = ch_versions.mix(INDEX_GENOME.out.versions)

            // Set channels: clean fasta and its index
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
        ch_reference_mature  = params.mature  ? Channel.fromPath(params.mature, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect()  : { exit 1, "Specify a Mature miRNA fasta file via '--mature'" }
        ch_reference_hairpin = params.hairpin ? Channel.fromPath(params.hairpin, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() : { exit 1, "Specify a Hairpin miRNA fasta file via '--hairpin'" }
    } else {
        if (!params.mirgenedb_species) {
            exit 1, "MirGeneDB species not set, please specify via the --mirgenedb_species parameter"
        }
        ch_reference_mature  = params.mirgenedb_mature  ? Channel.fromPath(params.mirgenedb_mature, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect()  : { exit 1, "Specify a mirgenedb Mature miRNA fasta file via '--mirgenedb_mature'" }
        ch_reference_hairpin = params.mirgenedb_hairpin ? Channel.fromPath(params.mirgenedb_hairpin, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() : { exit 1, "Specify a mirgenedb Hairpin miRNA fasta file via '--mirgenedb_hairpin'" }
        ch_mirna_gtf         = params.mirgenedb_gff     ? Channel.fromPath(params.mirgenedb_gff, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect()     : { exit 1, "Specify a MirGeneDB gff file via '--mirgenedb_gff'"}

        // Create a channel for mirgenedb_species
        ch_mirgenedb_species = Channel.value(params.mirgenedb_species)

        // Add species of the gtf
        // When mirgenedb workflow is not indicated, species defaults to val_mirtrace_species.
        // If mirgenedb workflow parameters are indicated, the params.mirgenedb_species is used instead.
        // If both mirgenedb workflow parameters and mirtrace_species (or mirna_gtf) are provided, params.mirgenedb_species is used as species value
        ch_mirna_gtf = ch_mirna_gtf
            .combine(ch_mirgenedb_species)
            .map { meta, gtf, species ->
                def new_meta = meta.clone() + [species: species]
                [new_meta, gtf]
            }
    }

    emit:
    fasta                 = ch_fasta               // channel: [ val(meta), path(fasta) ]
    has_fasta             = bool_has_fasta         // boolean
    bowtie_index          = ch_bowtie_index        // channel: [ val(meta), path(directory_index) ]
    versions              = ch_versions            // channel: [ versions.yml ]
    mirtrace_species      = ch_mirtrace_species    // channel: [ val(string) ]
    has_mirtrace_species  = bool_mirtrace_species  // boolean
    reference_mature      = ch_reference_mature    // channel: [ val(meta), path(fasta) ]
    reference_hairpin     = ch_reference_hairpin   // channel: [ val(meta), path(fasta) ]
    mirna_gtf             = ch_mirna_gtf           // channel: [ val(meta), path(gtf) ]
    rrna                  = ch_rrna                // channel: [ val(meta), path(fasta) ]
    trna                  = ch_trna                // channel: [ val(meta), path(fasta) ]
    cdna                  = ch_cdna                // channel: [ val(meta), path(fasta) ]
    ncrna                 = ch_ncrna               // channel: [ val(meta), path(fasta) ]
    pirna                 = ch_pirna               // channel: [ val(meta), path(fasta) ]
    other_contamination   = ch_other_contamination // channel: [ val(meta), path(fasta) ]
    mirna_adapters        = ch_mirna_adapters      // channel: [ val(meta), path(fasta) ]
}
