/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSmrnaseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
def trimgalore_options    = modules['trimgalore']
// TODO if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }
if (params.mature) { reference_mature = file(params.mature, checkIfExists: true) } else { exit 1, "Mature miRNA fasta file not found: ${params.mature}" }
if (params.hairpin) { reference_hairpin = file(params.hairpin, checkIfExists: true) } else { exit 1, "Hairpin miRNA fasta file not found: ${params.hairpin}" }
mirna_gtf = params.mirtrace_species ? file("ftp://mirbase.org/pub/mirbase/CURRENT/genomes/${params.mirtrace_species}.gff3", checkIfExists: true) : false

include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )
include { FASTQC_TRIMGALORE } from '../subworkflows/nf-core/fastqc_trimgalore' addParams( fastqc_options: modules['fastqc'], trimgalore_options: trimgalore_options )
include { MIRNA_QUANT } from '../subworkflows/local/mirna_quant' addParams( samtools_options: modules['samtools_view'], map_options: modules['map_mirna'],
                                                                            samtools_sort_options: modules['samtools_sort'],
                                                                            samtools_index_options: modules['samtools_index'],
                                                                            samtools_stats_options: modules['samtools_index'] )



/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''
def cat_fastq_options          = modules['cat_fastq']
//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ } from '../modules/nf-core/modules/cat/fastq/main' addParams( options: cat_fastq_options   )
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SMRNASEQ {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    ).map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .dump(tag: 'map')
    .groupTuple(by: [0])
    .dump(tag: 'group')
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }

    // //
    // // MODULE: Run FastQC
    // //
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    //
    // MODULE: Run trimming
    // //
    // INPUT_TRIM (
    //     INPUT_CHECK.out.reads
    // )
    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .mix(ch_fastq.single)
    .dump(tag:'cat')
    .set { ch_cat_fastq }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))


    MIRNA_QUANT (
        reference_mature,
        reference_hairpin,
        mirna_gtf,
        FASTQC_TRIMGALORE.out.reads
    )

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSmrnaseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
