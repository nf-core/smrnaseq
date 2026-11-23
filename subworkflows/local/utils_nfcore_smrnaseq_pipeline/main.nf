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
        '-m'           : 'int',
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
    if (expectedType == 'int_min_2' && !(value ==~ /0*([2-9]|[1-9]\d+)/)) {
        error("[bowtie_align_args] Argument '${flag}' must be at least 2, got '${value}'. The default genome mapping arguments include '--strata', which Bowtie only accepts with '-k' above 1. To report a single alignment per read, override ext.args for 'NFCORE_SMRNASEQ:GENOME_QUANT:BOWTIE_MAP_GENOME' instead: https://nf-co.re/smrnaseq/dev/docs/usage#custom-tool-arguments")
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
