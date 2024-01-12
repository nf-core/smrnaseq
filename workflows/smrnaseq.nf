/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation


WorkflowSmrnaseq.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

WorkflowSmrnaseq.initialise(params, log)

// Check optional parameters
if (!params.mirtrace_species) {
    exit 1, "Reference species for miRTrace is not defined via the --mirtrace_species parameter."
}

// Genome options
def mirna_gtf_from_species = params.mirtrace_species ? "https://mirbase.org/download/CURRENT/genomes/${params.mirtrace_species}.gff3" : false
def mirna_gtf = params.mirna_gtf ?: mirna_gtf_from_species

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
ch_fastp_adapters          = Channel.fromPath(params.fastp_known_mirna_adapters, checkIfExists: true).collect() // collect to consume for all incoming samples to FASTP

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
if (!params.mirgenedb) {
    if (params.mature) { reference_mature = file(params.mature, checkIfExists: true) } else { exit 1, "Mature miRNA fasta file not found: ${params.mature}" }
    if (params.hairpin) { reference_hairpin = file(params.hairpin, checkIfExists: true) } else { exit 1, "Hairpin miRNA fasta file not found: ${params.hairpin}" }
} else {
    if (params.mirgenedb_mature) { reference_mature = file(params.mirgenedb_mature, checkIfExists: true) } else { exit 1, "Mature miRNA fasta file not found: ${params.mirgenedb_mature}" }
    if (params.mirgenedb_hairpin) { reference_hairpin = file(params.mirgenedb_hairpin, checkIfExists: true) } else { exit 1, "Hairpin miRNA fasta file not found: ${params.mirgenedb_hairpin}" }
    if (params.mirgenedb_gff) { mirna_gtf = file(params.mirgenedb_gff, checkIfExists: true) } else { exit 1, "MirGeneDB gff file not found: ${params.mirgenedb_gff}"}
}

include { INPUT_CHECK           } from '../subworkflows/local/input_check'
include { FASTQ_FASTQC_UMITOOLS_FASTP } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp'
include { DEDUPLICATE_UMIS      } from '../subworkflows/local/umi_dedup'
include { CONTAMINANT_FILTER    } from '../subworkflows/local/contaminant_filter'
include { MIRNA_QUANT           } from '../subworkflows/local/mirna_quant'
include { GENOME_QUANT          } from '../subworkflows/local/genome_quant'
include { MIRTRACE              } from '../subworkflows/local/mirtrace'
include { MIRDEEP2              } from '../subworkflows/local/mirdeep2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SMRNASEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    .reads
    .dump(tag: 'group')
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters & dedup UMIs if necessary / desired by the user
    //

    FASTQ_FASTQC_UMITOOLS_FASTP (
        ch_cat_fastq,
        params.skip_fastqc,
        params.with_umi,
        params.skip_umi_extract,
        params.umi_discard_read,
        params.skip_fastp,
        params.fastp_known_mirna_adapters,
        params.save_trimmed_fail,
        params.save_merged,
        params.min_trimmed_reads
    )
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)

    reads_for_mirna = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads

    //
    // SUBWORKFLOW: Deduplicate UMIs by mapping them to the genome
    //
    if (params.with_umi){
        if (fasta){
            fasta_ch = file(fasta)
            DEDUPLICATE_UMIS (
                fasta_ch,
                bt_index,
                FASTQC_UMITOOLS_FASTP.out.trim_reads
            )
            reads_for_mirna = DEDUPLICATE_UMIS.out.reads
            ch_versions = ch_versions.mix(DEDUPLICATE_UMIS.out.versions)
        }
    }


    //
    // SUBWORKFLOW: mirtrace QC
    //
    FASTQ_FASTQC_UMITOOLS_FASTP.out.adapter_seq
    .join( FASTQ_FASTQC_UMITOOLS_FASTP.out.reads )
    .map { meta, adapter_seq, reads -> [adapter_seq, meta.id, reads] }
    .groupTuple()
    .set { ch_mirtrace_inputs }

    MIRTRACE(ch_mirtrace_inputs)
    ch_versions = ch_versions.mix(MIRTRACE.out.versions.ifEmpty(null))


    //
    // SUBWORKFLOW: remove contaminants from reads
    //
    contamination_stats = Channel.empty()
    mirna_reads = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
    if (params.filter_contamination){
        CONTAMINANT_FILTER (
            reference_hairpin,
            params.rrna,
            params.trna,
            params.cdna,
            params.ncrna,
            params.pirna,
            params.other_contamination,
            FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        )

        contamination_stats = CONTAMINANT_FILTER.out.filter_stats
        ch_versions = ch_versions.mix(CONTAMINANT_FILTER.out.versions)
        mirna_reads = CONTAMINANT_FILTER.out.filtered_reads

    }

    MIRNA_QUANT (
        [ [:], reference_mature],
        [ [:], reference_hairpin],
        mirna_gtf,
        mirna_reads
    )
    ch_versions = ch_versions.mix(MIRNA_QUANT.out.versions.ifEmpty(null))

    //
    // GENOME
    //
    genome_stats = Channel.empty()
    if (params.fasta){
        GENOME_QUANT ( file(params.fasta), params.bowtie_index, MIRNA_QUANT.out.unmapped )
        genome_stats = GENOME_QUANT.out.stats
        ch_versions = ch_versions.mix(GENOME_QUANT.out.versions)

        hairpin_clean = MIRNA_QUANT.out.fasta_hairpin.map { it -> it[1] }
        mature_clean  = MIRNA_QUANT.out.fasta_mature.map { it -> it[1] }

        if (!params.skip_mirdeep) {
            MIRDEEP2 (
                FASTQ_FASTQC_UMITOOLS_FASTP.out.reads,
                GENOME_QUANT.out.fasta,
                GENOME_QUANT.out.index.collect(),
                hairpin_clean,
                mature_clean
            )
            ch_versions = ch_versions.mix(MIRDEEP2.out.versions)
        }
    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowSmrnaseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowSmrnaseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description = Channel.value(methods_description)
        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(contamination_stats.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(genome_stats.collect({it[1]}).ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.mature_stats.collect({it[1]}).ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.hairpin_stats.collect({it[1]}).ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.mirtop_logs.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRTRACE.out.results.collect().ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
