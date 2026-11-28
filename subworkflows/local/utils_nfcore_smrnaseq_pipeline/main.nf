//
// Subworkflow with functionality specific to the nf-core/smrnaseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsHelp           } from 'plugin/nf-schema'
include { paramsSummaryLog     } from 'plugin/nf-schema'
include { paramsSummaryMap     } from 'plugin/nf-schema'
include { samplesheetToList    } from 'plugin/nf-schema'
include { validateParameters   } from 'plugin/nf-schema'
include { checkCondaChannels   } from 'plugin/nf-core-utils'
include { checkConfigProvided  } from 'plugin/nf-core-utils'
include { checkProfileProvided } from 'plugin/nf-core-utils'
include { completionEmail      } from 'plugin/nf-core-utils'
include { completionSummary    } from 'plugin/nf-core-utils'
include { dumpParametersToJSON } from 'plugin/nf-core-utils'
include { getWorkflowVersion   } from 'plugin/nf-core-utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version                 // boolean: Display version and exit
    validate_params         // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs         // boolean: Do not use coloured log outputs
    nextflow_cli_args       //   array: List of positional nextflow CLI args
    outdir                  //  string: The output directory where the results will be saved
    input                   //  string: Path to input samplesheet
    val_three_prime_adapter //  string: Sequencing adapter sequence to use for trimming
    val_phred_offset        //  string: The PHRED quality offset to be used for any input fastq files
    help                    // boolean: Display help message and exit
    help_full               // boolean: Show the full help message
    show_hidden             // boolean: Show hidden parameters in the help message

    main:

    ch_versions            = channel.empty()
    ch_three_prime_adapter = channel.value(val_three_prime_adapter)
    ch_phred_offset        = channel.value(val_phred_offset)

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    if (version) {
        log.info("${workflow.manifest.name} ${getWorkflowVersion()}")
        System.exit(0)
    }
    if (outdir) {
        dumpParametersToJSON(outdir.toString(), params)
    }
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba'])) {
        checkCondaChannels()
    }

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/smrnaseq ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/smrnaseq/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    if (help || help_full) {
        log.info(
            paramsHelp(
                [beforeText: before_text, afterText: after_text, command: command, showHidden: show_hidden, fullHelp: help_full],
                (help instanceof String && help != 'true') ? help : ''
            )
        )
        System.exit(0)
    }

    log.info(before_text)
    log.info(paramsSummaryLog(workflow))
    log.info(after_text)
    if (validate_params) {
        validateParameters()
    }

    //
    // Check config provided to the pipeline
    //
    checkConfigProvided()
    checkProfileProvided(nextflow_cli_args, monochrome_logs)

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Validate the bowtie align arguments provided through params.bowtie_align_args
    //
    validateBowtieGenomeArgs()

    //
    // Warn if there is no genome FASTA, so genome quantification (and miRDeep2) are skipped
    //
    warnMissingGenomeFasta(params.bowtie_index, params.fasta, params.skip_mirdeep)

    //
    // Create channel from input file provided through params.input
    //

    ch_samplesheet = channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }

    emit:
    samplesheet         = ch_samplesheet         // channel: sample fastqs parsed from --input
    versions            = ch_versions            // channel: [ versions.yml ]
    three_prime_adapter = ch_three_prime_adapter // channel: [ val(string) ]
    phred_offset        = ch_phred_offset        // channel: [ val(string) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.isBound() ? multiqc_reports.getVal() : [],
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

//
// miRNA reference database version for the software versions YAML: read from the miRNA GFF3
// header, or, when no GFF3 is used (e.g. mirtop is skipped and the GFF channel is empty),
// derived from the database and the mature/hairpin FASTA files. Never fails the run for a
// missing GFF3.
//
workflow MIRNA_REFERENCE_VERSIONS {

    take:
    ch_mirna_gff  // channel: [ val(meta), path(gff) ], may be empty
    val_mirgenedb // boolean: MirGeneDB (true) or miRBase (false)
    val_mature    //  string: mature miRNA FASTA path or URL
    val_hairpin   //  string: hairpin miRNA FASTA path or URL

    main:
    def database = val_mirgenedb ? 'MirGeneDB' : 'miRBase'
    def ch_versions = ch_mirna_gff
        .map { _meta, gff -> mirnaReferenceVersions(gff) }
        .ifEmpty { mirnaReferenceFallbackVersions(database, val_mature, val_hairpin) }
        .flatMap { versions -> versions }

    emit:
    versions = ch_versions // channel: [ process, tool, version ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()

    if (!params.mirgenedb) {
        // Validate mature miRNA fasta file
        if (!params.mature) {
            error("Mature miRNA fasta file not found. Please specify using the `--mature` parameter.")
        }
        // Validate hairpin miRNA fasta file
        if (!params.hairpin) {
            error("Hairpin miRNA fasta file not found. Please specify using the `--hairpin` parameter.")
        }
    } else {
        // Validate MirGeneDB species
        if (!params.mirgenedb_species) {
            error("You specified to be using MirGeneDB, but the MirGeneDB species is not set. Please specify using the `--mirgenedb_species` parameter.")
        }
        // Validate MirGeneDB mature miRNA fasta file
        if (!params.mirgenedb_mature) {
            error("You specified to be using MirGeneDB, but the mature miRNA fasta file is not found. Please provide the file using the `--mirgenedb_mature` parameter.")
        }
        // Validate MirGeneDB hairpin miRNA fasta file
        if (!params.mirgenedb_hairpin) {
            error("You specified to be using MirGeneDB, but the hairpin miRNA fasta file is not found. Please provide the file using the `--mirgenedb_hairpin` parameter.")
        }
        // Validate MirGeneDB GFF file
        if (!params.mirgenedb_gff) {
            error("You specified to be using MirGeneDB, but the GFF file is not found. Please provide the file using the `--mirgenedb_gff` parameter.")
        }
    }

    // Validate UMI parameters
    def umi_error = umiParameterError(params.with_umi, params.umitools_bc_pattern)
    if (umi_error) {
        error(umi_error)
    }
    umiDeprecationWarnings(params.skip_umi_extract_before_dedup, params.umi_discard_read).each { message ->
        log.warn(message)
    }
    def umi_dedup_method_warning = umiDedupMethodWarning(params.with_umi, params.umi_dedup_method)
    if (umi_dedup_method_warning) {
        log.warn(umi_dedup_method_warning)
    }
}

//
// Warning for `--umi_dedup_method fgumi`: the method is experimental, and it has no effect without
// `--with_umi`. Returns a warning message, or null for the default method.
//
def umiDedupMethodWarning(with_umi, umi_dedup_method) {
    if (umi_dedup_method != 'fgumi') {
        return null
    }
    if (!with_umi) {
        return "'--umi_dedup_method fgumi' has no effect without '--with_umi': reads are not deduplicated."
    }
    return "'--umi_dedup_method fgumi' is experimental: UMI deduplication after alignment to a mature miRNA windows reference with fgumi (exact UMI matching). Reads that do not align to this reference are not deduplicated, '--umitools_method' is not used and no UMI saturation statistics are written. Its behaviour and outputs may change in future releases. See https://nf-co.re/smrnaseq/usage#experimental-alignment-based-deduplication-with-fgumi"
}

//
// Check UMI parameter combinations that cannot work. Returns an error message, or null if the
// combination is valid (or `--with_umi` is not set).
//
def umiParameterError(with_umi, bc_pattern) {
    if (with_umi && !bc_pattern) {
        return "'--with_umi' requires '--umitools_bc_pattern' to describe where the UMI is in the read (-profile qiaseq sets it for QIAseq miRNA libraries). See https://nf-co.re/smrnaseq/usage#umi-handling"
    }
    return null
}

//
// Warnings for deprecated UMI parameters that are set (they are ignored). Returns a list of messages.
//
def umiDeprecationWarnings(skip_umi_extract_before_dedup, umi_discard_read) {
    def messages = []
    if (skip_umi_extract_before_dedup != null) {
        messages << "'--skip_umi_extract_before_dedup' is deprecated and ignored: with '--with_umi', the UMI is always extracted before deduplication. It will be removed in a future release."
    }
    if (umi_discard_read != null) {
        messages << "'--umi_discard_read' is deprecated and ignored: the pipeline only uses read 1. It will be removed in a future release."
    }
    return messages
}

//
// Warning for UMI saturation, from the statistics TSV written by UMI_FROM_SEQUENCE (header line and
// one line per sample). One-mismatch UMI clustering (UMICollapse `dir`, `adj` or `cc`) merges distinct
// molecules once one insert sequence uses more than about 0.5% of all possible UMIs (4^UMI length),
// so its count is underestimated. Returns a warning message, or null.
//
def umiSaturationWarning(stats_text) {
    def threshold = 0.005
    def lines = stats_text.readLines().findAll { line -> line.trim() }
    if (lines.size() < 2) {
        return null
    }
    def header = lines[0].split('\t') as List
    def values = lines[1].split('\t') as List
    def row = [header, values].transpose().collectEntries()
    def fraction = row.top_sequence_umi_fraction as BigDecimal
    if (fraction <= threshold) {
        return null
    }
    // Not String.format: its decimal separator depends on the locale
    def percent = (fraction * 100).setScale(2, BigDecimal.ROUND_HALF_UP).toPlainString()
    return "UMI saturation in sample '${row.sample}': after deduplication, the most frequent sequence (${row.top_sequence}) has ${row.top_sequence_reads} reads, ${percent}% of all ${row.umi_space} possible ${row.umi_length} nt UMIs (warning threshold 0.5%). One-mismatch UMI clustering then merges distinct molecules, so the counts of the most abundant miRNAs are underestimated. Consider exact UMI matching (UMICollapse '-k 0'). See https://nf-co.re/smrnaseq/usage#choosing-a-umi-strategy"
}

// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    // Emit a warning if `single_end` is false
    if (metas[0].single_end == false) {
        log.warn "Sample ${metas[0].id} is detected as paired-end reads (fastq_1 and fastq_2). The pipeline only handles SE data. Samplesheets with fastq_1 and fastq_2 are supported but fastq_2 is removed."
        // Remove fastq_2 from the list and keep only fastq_1
        fastqs = fastqs.collect { fastq -> fastq.take(1) }
        metas[0].single_end = true
    }

    return [ metas[0], fastqs ]
}
//
// Validate the bowtie align arguments provided through params.bowtie_align_args
//
def validateBowtieGenomeArgs() {
    def rawArgs = params.bowtie_align_args
    if (rawArgs == null || rawArgs == '') {
        return
    }
    // A bare '--bowtie_align_args --best' reaches here as true or 'true'
    if (!(rawArgs instanceof CharSequence) || rawArgs.toString() == 'true') {
        error("[bowtie_align_args] Expected a string of Bowtie options but got '${rawArgs}'. A value holding a single flag must use the '=' form, e.g. --bowtie_align_args=--best, because '--bowtie_align_args --best' is parsed as the boolean true.")
    }

    // Allowlist: flag -> expected value type (null = flag only)
    def allowedArgs = [
        '-v'           : 'int_0_3',
        '-n'           : 'int_0_3',
        '--seedmms'    : 'int_0_3',
        '-e'           : 'int',
        '--maqerr'     : 'int',
        '-l'           : 'int',
        '--seedlen'    : 'int',
        '-k'           : 'int_min_2',
        '-m'           : 'int_min_1',
        '-M'           : 'int_min_1',
        '--maxbts'     : 'int',
        '--chunkmbs'   : 'int',
        '--seed'       : 'int',
        '--best'       : null,
        '--strata'     : null,
        '--nomaqround' : null,
        '--norc'       : null,
    ]

    def seen = []

    // The accumulator holds the flag still waiting for its value, if any
    def unfinishedFlag = rawArgs.toString().tokenize().inject(null) { pendingFlag, token ->
        if (pendingFlag != null) {
            validateBowtieArgValue(pendingFlag, allowedArgs[pendingFlag], token)
            return null
        }

        def separator = token.indexOf('=')
        def flag = separator >= 0 ? token.substring(0, separator) : token

        if (!allowedArgs.containsKey(flag)) {
            error("[bowtie_align_args] Argument '${flag}' is not allowed. Only these flags are permitted: ${allowedArgs.keySet().sort().join(', ')}")
        }
        if (separator >= 0 && !flag.startsWith('--')) {
            // Bowtie does not reject '-v=1': it silently misreads the value
            error("[bowtie_align_args] Short flag '${flag}' cannot use the '=' form, got '${token}'. Use '${flag} <value>' instead.")
        }

        seen.add(flag)
        def expectedType = allowedArgs[flag]

        if (separator >= 0) {
            if (expectedType == null) {
                error("[bowtie_align_args] Argument '${flag}' does not take a value, got '${token}'.")
            }
            validateBowtieArgValue(flag, expectedType, token.substring(separator + 1))
            return null
        }
        return expectedType != null ? flag : null
    }

    if (unfinishedFlag != null) {
        error("[bowtie_align_args] Argument '${unfinishedFlag}' requires a value but none was provided.")
    }

    if (seen.contains('-v') && (seen.contains('-n') || seen.contains('--seedmms'))) {
        error("[bowtie_align_args] '-v' and '-n'/'--seedmms' are mutually exclusive.")
    }
    // Bowtie accepts both, but '-M' then wins and the last value of either sets the ceiling
    if (seen.contains('-m') && seen.contains('-M')) {
        error("[bowtie_align_args] '-m' and '-M' are mutually exclusive.")
    }
}
//
// Check the value given to a bowtie flag that takes an integer
//
def validateBowtieArgValue(flag, expectedType, value) {
    if (!(value ==~ /\d+/)) {
        error("[bowtie_align_args] Argument '${flag}' requires an integer value, got '${value}'.")
    }
    if (expectedType == 'int_0_3' && !(value ==~ /0*[0-3]/)) {
        error("[bowtie_align_args] Argument '${flag}' must be between 0 and 3, got '${value}'.")
    }
    if (expectedType == 'int_min_1' && !(value ==~ /0*[1-9]\d*/)) {
        error("[bowtie_align_args] Argument '${flag}' must be at least 1, got '${value}'.")
    }
    if (expectedType == 'int_min_2' && !(value ==~ /0*([2-9]|[1-9]\d+)/)) {
        error("[bowtie_align_args] Argument '${flag}' must be at least 2, got '${value}'. The default genome mapping arguments include '--strata', which Bowtie only accepts with '-k' above 1. To report a single alignment per read, use '-M 1' instead, which reports one random best alignment for reads with more than one.")
    }
}
//
// Build one line of the miRTrace config file: <reads>,<display name>,<adapter>,<PHRED offset>
// miRTrace (1.0.1) cannot auto-detect adapters, and an empty adapter field means no adapter trimming.
// `--three_prime_adapter auto-detect` only applies to fastp, so it becomes an empty field here.
//
def mirtraceConfigLine(reads_name, sample_id, adapter, phred_offset) {
    def mirtrace_adapter = adapter == 'auto-detect' ? '' : adapter
    return "./${reads_name},${sample_id},${mirtrace_adapter},${phred_offset}\n".toString()
}

//
// Warning for miRTrace when the adapter is 'auto-detect' and miRTrace runs, i.e. a miRTrace
// species is set (null otherwise).
//
def mirtraceAdapterWarning(adapter, mirtrace_species) {
    if (adapter != 'auto-detect' || !mirtrace_species) {
        return null
    }
    return "miRTrace cannot auto-detect adapters, so with '--three_prime_adapter auto-detect' it runs without adapter trimming. Set '--three_prime_adapter' to the adapter sequence for adapter-trimmed miRTrace QC."
}

//
// Warn when there is no genome FASTA. Genome quantification and miRDeep2 only run with a FASTA
// (`--fasta`, or `--genome` with an iGenomes entry), so without one they are skipped. A given
// genome Bowtie index is then ignored, which gets its own message. miRDeep2 is only mentioned
// when it would otherwise run. Logs at most one warning and returns it, or null if there is none.
//
def warnMissingGenomeFasta(bowtie_index, fasta, skip_mirdeep) {
    if (fasta) {
        return null
    }
    def skipped = skip_mirdeep ? 'genome quantification is skipped' : 'genome quantification and miRDeep2 are skipped'
    def pronoun = skip_mirdeep ? 'it' : 'them'
    def how     = "use '--genome' with an iGenomes reference that provides a FASTA."
    def message = bowtie_index
        ? "A genome Bowtie index was given with '--bowtie_index', but no genome FASTA. The index is ignored, and ${skipped}. " +
            "To run ${pronoun}, also provide the FASTA the index was built from with '--fasta', or ${how}"
        : "No genome FASTA was given, so ${skipped}. To run ${pronoun}, provide a genome FASTA with '--fasta', or ${how}"
    log.warn(message.toString())
    return message.toString()
}

//
// First lines of a text file, transparently decompressing gzip (detected by its magic bytes)
//
def readFirstLines(file_path, n_lines) {
    return file_path.withInputStream { raw ->
        def input = new java.io.BufferedInputStream(raw)
        input.mark(2)
        def is_gzip = input.read() == 0x1f && input.read() == 0x8b
        input.reset()
        def stream = is_gzip ? new java.util.zip.GZIPInputStream(input) : input
        def reader = new java.io.BufferedReader(new java.io.InputStreamReader(stream, 'UTF-8'))
        reader.iterator().take(n_lines).toList()
    }
}

//
// Version of the miRNA reference database and its genome build, read from the comment header
// of the miRNA GFF3 (--mirna_gtf, --mirgenedb_gff or the GFF3 downloaded for --mirtrace_species).
// miRBase and MirGeneDB GFF3 files start with lines such as:
//   # microRNAs:               miRBase v22
//   # genome-build-id:         GRCh38
// Returns [ process, tool, version ] tuples for the software versions YAML. The database
// version is "unknown" when the header has no "microRNAs:" line.
//
// Also checks the file: a GFF that is an HTML page (e.g. the miRBase website serving
// <p>##gff-version 3<br>... instead of the plain file) stops the run, since mirtop cannot read
// it. A file that does not start with the "##gff-version" pragma only gives a warning, as
// custom GFF files may lack it. Gzipped files are decompressed before these checks.
//
def mirnaReferenceVersions(gff) {
    def lines = readFirstLines(gff, 200)
    def first_line = lines.find { line -> line.trim() }?.trim() ?: ''
    def first_line_lower = first_line.toLowerCase()
    def is_html = first_line.startsWith('<') ||
        ['<html', '<!doctype', '<p>', '<br'].any { tag -> first_line_lower.contains(tag) }
    if (is_html) {
        error("The miRNA GFF file '${gff.toUriString()}' is an HTML page, not a GFF3 file (it starts with '${first_line.take(40)}'). " +
            "If it was downloaded from the miRBase website, download the plain GFF3 from https://mirbase.org/download/<species>.gff3 " +
            "(e.g. https://mirbase.org/download/mmu.gff3) and pass the local file with '--mirna_gtf' (or '--mirgenedb_gff' with '--mirgenedb').")
    }
    if (!first_line.startsWith('##gff-version')) {
        log.warn "The miRNA GFF file '${gff.toUriString()}' does not start with a '##gff-version' line. Please check that it is a GFF3 file with miRNA coordinates."
    }
    def header = lines.takeWhile { line -> line.startsWith('#') }
    def database = header.findResult { line ->
        def match = line =~ /^#+\s*microRNAs:\s*(\S+)\s+(\S+)/
        match.find() ? [match.group(1), match.group(2)] : null
    }
    def genome_build = header.findResult { line ->
        def match = line =~ /^#+\s*genome-build-id:\s*(\S+)/
        match.find() ? match.group(1) : null
    }
    def versions = database
        ? [['MIRNA_REFERENCE', database[0], database[1]]]
        : [['MIRNA_REFERENCE', 'mirna_gff', "unknown (${gff.name})".toString()]]
    if (genome_build) {
        versions << ['MIRNA_REFERENCE', 'genome_build', genome_build]
    }
    return versions
}

//
// miRNA reference entries for the software versions YAML when no miRNA GFF3 is available.
// The mature/hairpin FASTA headers carry no release, so the version is only known if the
// directory of the mature or hairpin file is a release number, e.g.
// https://mirbase.org/download_version_files/22.1/mature.fa. The file names are reported as well.
//
def mirnaReferenceFallbackVersions(database, mature, hairpin) {
    def release = [mature, hairpin].findResult { fasta ->
        def match = fasta ? (fasta.toString() =~ /(?:^|\/)v?(\d+(?:\.\d+)*)\/[^\/]+$/) : null
        match?.find() ? match.group(1) : null
    }
    def version = release ? "v${release} (from file path)".toString() : 'unknown (no GFF)'
    def file_names = [mature, hairpin].collect { fasta -> fasta ? fasta.toString().tokenize('/').last() : 'none' }
    return [
        ['MIRNA_REFERENCE', database, version],
        ['MIRNA_REFERENCE', 'mature', file_names[0]],
        ['MIRNA_REFERENCE', 'hairpin', file_names[1]]
    ]
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
