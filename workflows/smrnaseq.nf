/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSmrnaseq.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check optional parameters
if (!params.mirtrace_species){
    exit 1, "Reference species for miRTrace is not defined."
}
// Genome options
bt_index_from_species = params.genome ? params.genomes[ params.genome ].bowtie ?: false : false
bt_index              = params.bt_indices ?: bt_index_from_species
mirtrace_species_from_species = params.genome ? params.genomes[ params.genome ].mirtrace_species ?: false : false
mirtrace_species = params.mirtrace_species ?: mirtrace_species_from_species
fasta_from_species = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
fasta = params.fasta ?: fasta_from_species
mirna_gtf_from_species = params.mirtrace_species ? "https://mirbase.org/ftp/CURRENT/genomes/${params.mirtrace_species}.gff3" : false
mirna_gtf = params.mirna_gtf ? params.mirna_gtf : mirna_gtf_from_species

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
if (params.mature) { reference_mature = file(params.mature, checkIfExists: true) } else { exit 1, "Mature miRNA fasta file not found: ${params.mature}" }
if (params.hairpin) { reference_hairpin = file(params.hairpin, checkIfExists: true) } else { exit 1, "Hairpin miRNA fasta file not found: ${params.hairpin}" }

include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { FASTQC_TRIMGALORE } from '../subworkflows/nf-core/fastqc_trimgalore'
include { MIRNA_QUANT       } from '../subworkflows/local/mirna_quant'
include { GENOME_QUANT      } from '../subworkflows/local/genome_quant'
include { MIRTRACE          } from '../subworkflows/local/mirtrace'
include { MIRDEEP2          } from '../subworkflows/local/mirdeep2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

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
        ch_input
    )
    .reads
    .map {
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
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

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
    // SUBWORKFLOW: mirtrace QC
    //
    MIRTRACE (ch_cat_fastq)
    ch_versions = ch_versions.mix(MIRTRACE.out.versions.ifEmpty(null))


    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    reads_for_mirna = FASTQC_TRIMGALORE.out.reads
    MIRNA_QUANT (
        reference_mature,
        reference_hairpin,
        mirna_gtf,
        reads_for_mirna
    )
    ch_versions = ch_versions.mix(MIRNA_QUANT.out.versions.ifEmpty(null))

    //
    // GENOME
    //
    genome_stats = Channel.empty()
    if (fasta){
        fasta_ch = file(fasta)
        GENOME_QUANT ( fasta_ch, bt_index, MIRNA_QUANT.out.unmapped )
        GENOME_QUANT.out.stats
            .set { genome_stats }
        ch_versions = ch_versions.mix(GENOME_QUANT.out.versions)

        if (!params.skip_mirdeep) {
            MIRDEEP2 (FASTQC_TRIMGALORE.out.reads, GENOME_QUANT.out.fasta, GENOME_QUANT.out.indices, MIRNA_QUANT.out.fasta_hairpin, MIRNA_QUANT.out.fasta_mature)
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

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.mature_stats.collect({it[1]}).ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.hairpin_stats.collect({it[1]}).ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(genome_stats.collect({it[1]}).ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.mirtop_logs.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRTRACE.out.results.collect().ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect()
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
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
