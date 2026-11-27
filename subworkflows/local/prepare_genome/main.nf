//
// Uncompress and prepare reference genome files
//

// nf-core modules
include { UNTAR         as UNTAR_BOWTIE_INDEX } from '../../../modules/nf-core/untar'
include { BOWTIE_BUILD  as INDEX_GENOME       } from '../../../modules/nf-core/bowtie/build/main'
include { BIOAWK        as CLEAN_FASTA        } from '../../../modules/nf-core/bioawk/main'

// local functions

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Extract prefix from bowtie index files
//
def extractFirstIndexPrefix(files_path) {
    def files = files_path.listDirectory().toList()
    if (!files) {
        throw new Exception("The provided bowtie_index path doesn't contain any files.")
    }
    def index_file = files.find { file_path ->
        def file_name = file_path.getName()
        file_name.endsWith(".1.ebwt") && !file_name.endsWith(".rev.1.ebwt")
    }
    if (index_file == null) {
        throw new Exception("Unable to extract the prefix from the Bowtie index files. No file with the '.1.ebwt' extension was found. Please ensure that the correct files are in the specified path.")
    }
    def index_file_name = index_file.getName()
    return index_file_name.substring(0, index_file_name.lastIndexOf(".1.ebwt"))
}

//
// Error message when the miRBase GFF3 derived from --mirtrace_species cannot be found
//
def mirbaseGffMissingMessage(species, url) {
    return "The miRBase GFF3 file for --mirtrace_species '${species}' could not be found at ${url}. " +
        "miRBase only provides GFF3 files for species with genome coordinates, so it is missing for many species (e.g. sheep, 'oar'), " +
        "or miRBase may be temporarily unavailable. The GFF3 file is only needed for the mirtop isomiR annotation. " +
        "Either provide a miRBase-style GFF3 file via '--mirna_gtf', or run with '--skip_mirtop' to skip mirtop " +
        "and keep the mature and hairpin miRNA quantification. See https://nf-co.re/smrnaseq/usage#species-without-a-mirbase-gff3-file"
}

//
// URL of the miRBase GFF3 file used when `--mirna_gtf` is not provided.
// Uses the raw `https://mirbase.org/download/<species>.gff3` endpoint: the
// `https://mirbase.org/download/CURRENT/genomes/<species>.gff3` page serves the file
// wrapped in HTML (`<p>##gff-version 3<br>...`), so mirtop cannot find the
// `# microRNAs: miRBase vXX` header line and fails with "Database not found" (#284).
//
def getMirbaseGffUrl(species) {
    if (!species) {
        return null
    }
    if (species == 'hsa') {
        return "https://raw.githubusercontent.com/nf-core/test-datasets/smrnaseq/reference/hsa.gff3"
    }
    return "https://mirbase.org/download/${species}.gff3"
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
    ch_versions = channel.empty()
    ch_bowtie_index_untar = channel.empty()
    ch_bowtie_index_built = channel.empty()

    // Parameter channel handling
    ch_fasta                  = val_fasta                     ? channel.fromPath(val_fasta, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect()  : channel.empty()
    ch_bowtie_index           = val_bowtie_index              ? channel.fromPath(val_bowtie_index, checkIfExists: true).map{ it -> [ [id:""], it ] }.collect()    : channel.empty()

    bool_mirtrace_species     = val_mirtrace_species          ? true : false
    bool_has_fasta            = val_fasta                     ? true : false

    ch_mirtrace_species       = val_mirtrace_species          ? channel.value(val_mirtrace_species) : channel.empty()
    // The miRNA GFF3 is only used by mirtop: it is not downloaded from miRBase with --skip_mirtop, --mirna_gtf or --mirgenedb
    mirna_gtf_from_species    = (!val_mirna_gtf && !params.skip_mirtop && !params.mirgenedb) ? getMirbaseGffUrl(val_mirtrace_species) : null
    if (mirna_gtf_from_species && !file(mirna_gtf_from_species).exists()) {
        error(mirbaseGffMissingMessage(val_mirtrace_species, mirna_gtf_from_species))
    }
    ch_mirna_gtf              = params.skip_mirtop            ? channel.empty() : ( val_mirna_gtf ? channel.fromPath(val_mirna_gtf, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() : ( mirna_gtf_from_species ? channel.fromPath(mirna_gtf_from_species, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() :  channel.empty() ) )

    // Add species of the gtf in the meta
    ch_mirna_gtf = ch_mirna_gtf
        .combine(ch_mirtrace_species.ifEmpty('unknown'))
        .map { meta, gtf, species ->
            def new_meta = meta.clone() + [species: species]
            [new_meta, gtf]
        }

    ch_mirna_adapters         = params.with_umi  ? [] : channel.fromPath(val_fastp_known_mirna_adapters, checkIfExists: true).collect()

    ch_rrna                   = val_rrna                      ? channel.fromPath(val_rrna, checkIfExists: true).map{ it -> [ [id:'rRNA'], it ] }.collect()                 : channel.empty()
    ch_trna                   = val_trna                      ? channel.fromPath(val_trna, checkIfExists: true).map{ it -> [ [id:'tRNA'], it ] }.collect()                 : channel.empty()
    ch_cdna                   = val_cdna                      ? channel.fromPath(val_cdna, checkIfExists: true).map{ it -> [ [id:'cDNA'], it ] }.collect()                 : channel.empty()
    ch_ncrna                  = val_ncrna                     ? channel.fromPath(val_ncrna, checkIfExists: true).map{ it -> [ [id:'ncRNA'], it ] }.collect()               : channel.empty()
    ch_pirna                  = val_pirna                     ? channel.fromPath(val_pirna, checkIfExists: true).map{ it -> [ [id:'piRNA'], it ] }.collect()               : channel.empty()
    ch_other_contamination    = val_other_contamination       ? channel.fromPath(val_other_contamination, checkIfExists: true).map{ it -> [ [id:'other'], it ] }.collect() : channel.empty()

    // even if bowtie index is specified, there still needs to be a fasta.
    // without fasta, no genome analysis.
    if(val_fasta) {
        // Clean fasta (replace non-ATCGs with Ns, remove whitespaces from headers)
        // Note: CLEAN_FASTA runs even when a bowtie_index is provided, as cleaning doesn't affect it, making regeneration unnecessary.
        CLEAN_FASTA ( ch_fasta, [], false, 'fa' )
        ch_fasta         = CLEAN_FASTA.out.output

        //Prepare bowtie index, unless specified
        //This needs to be done here as the index is used by GENOME_QUANT
        if(val_bowtie_index) {
            if (val_bowtie_index.toString().endsWith(".tar.gz")) {
                UNTAR_BOWTIE_INDEX ( ch_bowtie_index )
                ch_bowtie_index_untar = UNTAR_BOWTIE_INDEX.out.untar
                ch_bowtie_index = UNTAR_BOWTIE_INDEX.out.untar
                    .map{ _meta, index_dir ->
                        def index_prefix = extractFirstIndexPrefix(index_dir)
                        [[id:index_prefix], index_dir]
                    }
            } else {
                ch_bowtie_index = channel.fromPath(val_bowtie_index, checkIfExists: true)
                    .map{it ->
                        def index_prefix = extractFirstIndexPrefix(it)
                        [[id:index_prefix], it]
                    }
            }

        } else {

            // Index FASTA with nf-core Bowtie1
            INDEX_GENOME ( CLEAN_FASTA.out.output )

            // Set channels: clean fasta and its index
            ch_bowtie_index  = INDEX_GENOME.out.index.collect()
            ch_bowtie_index_built = INDEX_GENOME.out.index
        }
    }

    //Config checks
    // Check optional parameters
    if (!params.mirgenedb && !val_mirtrace_species && !(params.help || params.help_full)) {
            error("Reference species is not defined via the --mirtrace_species parameter. Set it to the miRBase species code (e.g. 'hsa', 'mmu', 'oar'), which is used to select the species' miRNAs from --mature and --hairpin, for miRTrace and for the miRBase GFF3 download.")
        }

    // Genome options
    if (!params.mirgenedb) {
        ch_reference_mature  = params.mature  ? channel.fromPath(params.mature, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect()  : error("Specify a Mature miRNA fasta file via '--mature'")
        ch_reference_hairpin = params.hairpin ? channel.fromPath(params.hairpin, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() : error("Specify a Hairpin miRNA fasta file via '--hairpin'")
    } else {
        if (!params.mirgenedb_species) {
            error("MirGeneDB species not set, please specify via the --mirgenedb_species parameter")
        }
        ch_reference_mature  = params.mirgenedb_mature  ? channel.fromPath(params.mirgenedb_mature, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect()  : error("Specify a mirgenedb Mature miRNA fasta file via '--mirgenedb_mature'")
        ch_reference_hairpin = params.mirgenedb_hairpin ? channel.fromPath(params.mirgenedb_hairpin, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() : error("Specify a mirgenedb Hairpin miRNA fasta file via '--mirgenedb_hairpin'")
        ch_mirna_gtf         = params.skip_mirtop       ? channel.empty() : ( params.mirgenedb_gff ? channel.fromPath(params.mirgenedb_gff, checkIfExists: true).map{ it -> [ [id:it.baseName], it ] }.collect() : error("Specify a MirGeneDB gff file via '--mirgenedb_gff', or skip mirtop with '--skip_mirtop'") )

        // Create a channel for mirgenedb_species
        ch_mirgenedb_species = channel.value(params.mirgenedb_species)

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
    bowtie_index_untar    = ch_bowtie_index_untar  // channel: [ val(meta), path(directory_index) ] (only when a .tar.gz index is provided)
    bowtie_index_built    = ch_bowtie_index_built  // channel: [ val(meta), path(directory_index) ] (only when the index is built by the pipeline)
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
