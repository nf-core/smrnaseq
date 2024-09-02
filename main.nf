#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/smrnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/smrnaseq
    Website: https://nf-co.re/smrnaseq
    Slack  : https://nfcore.slack.com/channels/smrnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NFCORE_SMRNASEQ         } from './workflows/smrnaseq'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_smrnaseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_smrnaseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_smrnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta               = getGenomeAttribute('fasta')
params.mirtrace_species    = getGenomeAttribute('mirtrace_species')
params.bowtie_index        = getGenomeAttribute('bowtie')
params.mirna_gtf           = getGenomeAttribute('mirna_gtf') //not in igenomes yet
params.rrna                = getGenomeAttribute('rrna') //not in igenomes yet
params.trna                = getGenomeAttribute('trna') //not in igenomes yet
params.cdna                = getGenomeAttribute('cdna') //not in igenomes yet
params.ncrna               = getGenomeAttribute('ncrna') //not in igenomes yet
params.pirna               = getGenomeAttribute('pirna') //not in igenomes yet
params.other_contamination = getGenomeAttribute('other_contamination') //not in igenomes yet


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    main:
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW : Prepare reference genome files
    //
    PREPARE_GENOME (
        params.fasta,
        params.bowtie_index,
        params.mirtrace_species,
        params.rrna,
        params.trna,
        params.cdna,
        params.ncrna,
        params.pirna,
        params.other_contamination,
        params.fastp_known_mirna_adapters,
        params.mirna_gtf
    )

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_SMRNASEQ (
        PREPARE_GENOME.out.has_fasta,
        PREPARE_GENOME.out.has_mirtrace_species,
        PREPARE_GENOME.out.mirna_adapters,
        PREPARE_GENOME.out.mirtrace_species,
        PREPARE_GENOME.out.reference_mature,
        PREPARE_GENOME.out.reference_hairpin,
        PREPARE_GENOME.out.mirna_gtf,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.bowtie_index,
        PREPARE_GENOME.out.rrna,
        PREPARE_GENOME.out.trna,
        PREPARE_GENOME.out.cdna,
        PREPARE_GENOME.out.ncrna,
        PREPARE_GENOME.out.pirna,
        PREPARE_GENOME.out.other_contamination,
        ch_versions,
        PIPELINE_INITIALISATION.out.samplesheet,
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_SMRNASEQ.out.multiqc_report
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
