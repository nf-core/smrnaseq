/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSmrnaseq.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check optional parameters
if( !params.mirtrace_species ){
    exit 1, "Reference species for miRTrace is not defined."
}
// Genome options
bt_index_from_species = params.genome ? params.genomes[ params.genome ].bowtie ?: false : false
bt_index              = params.bt_indices ?: bt_index_from_species
mirtrace_species_from_species = params.genome ? params.genomes[ params.genome ].mirtrace_species ?: false : false
mirtrace_species = params.mirtrace_species ?: mirtrace_species_from_species
fasta_from_species = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
fasta = params.fasta ?: fasta_from_species
mirna_gtf_from_species = params.mirtrace_species ? "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/${params.mirtrace_species}.gff3" : false
mirna_gtf = params.mirna_gtf ? params.mirna_gtf : mirna_gtf_from_species

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

include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )
include { FASTQC_TRIMGALORE } from '../subworkflows/nf-core/fastqc_trimgalore' addParams( fastqc_options: modules['fastqc'], trimgalore_options: trimgalore_options )
include { MIRNA_QUANT } from '../subworkflows/local/mirna_quant' addParams( samtools_options: modules['samtools_view'], map_options: modules['map_mirna'],
                                                                            samtools_sort_options: modules['samtools_sort'],
                                                                            samtools_index_options: modules['samtools_index'],
                                                                            samtools_stats_options: modules['samtools_index'],
                                                                            table_merge_options: modules['table_merge'] )
include { GENOME_QUANT } from '../subworkflows/local/genome_quant' addParams( samtools_options: modules['samtools_view'], map_options: modules['map_mirna'],
                                                                            samtools_sort_options: modules['samtools_sort'],
                                                                            samtools_index_options: modules['samtools_index'],
                                                                            samtools_stats_options: modules['samtools_index'] )
include { MIRTRACE } from '../subworkflows/local/mirtrace'
include { MIRDEEP2 } from '../subworkflows/local/mirdeep2'


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

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    CAT_FASTQ.out.reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }

    //
    // SUBWORKFLOW: mirtrace QC
    //
    MIRTRACE (ch_cat_fastq)
    ch_software_versions = ch_software_versions.mix(MIRTRACE.out.versions.ifEmpty(null))


    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.fastqc_versions.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.trimgalore_versions.first().ifEmpty(null))

    reads_for_mirna = FASTQC_TRIMGALORE.out.reads
    MIRNA_QUANT (
        reference_mature,
        reference_hairpin,
        mirna_gtf,
        reads_for_mirna
    )
    ch_software_versions = ch_software_versions.mix(MIRNA_QUANT.out.bowtie_versions.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(MIRNA_QUANT.out.samtools_versions.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(MIRNA_QUANT.out.seqcluster_versions.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(MIRNA_QUANT.out.mirtop_versions.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(MIRNA_QUANT.out.merge_versions.ifEmpty(null))

    //
    // GENOME
    //
    genome_stats = Channel.empty()
    if (fasta){
        fasta_ch = file(fasta)
        GENOME_QUANT ( fasta_ch, bt_index, MIRNA_QUANT.out.unmapped )
        GENOME_QUANT.out.stats
            .set { genome_stats }
        MIRDEEP2 (FASTQC_TRIMGALORE.out.reads, GENOME_QUANT.out.fasta , GENOME_QUANT.out.indices, MIRNA_QUANT.out.fasta_hairpin, MIRNA_QUANT.out.fasta_mature)
        ch_software_versions = ch_software_versions.mix(MIRDEEP2.out.versions.first().ifEmpty(null))

    }

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
    ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.mature_stats.collect({it[1]}).ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.hairpin_stats.collect({it[1]}).ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(genome_stats.collect({it[1]}).ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.mirtop_logs.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MIRTRACE.out.results.collect().ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions.ifEmpty(null))
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
