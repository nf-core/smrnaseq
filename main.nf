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

params.fasta            = getGenomeAttribute('fasta')
params.mirtrace_species = getGenomeAttribute('mirtrace_species')
params.bowtie_index     = getGenomeAttribute('bowtie')


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
        params.input,
        params.mirtrace_species
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_SMRNASEQ (
        Channel.of(file(params.input, checkIfExists: true)),
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.mirtrace_species,
        PIPELINE_INITIALISATION.out.reference_mature,
        PIPELINE_INITIALISATION.out.reference_hairpin,
        PIPELINE_INITIALISATION.out.mirna_gtf,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.bowtie_index,
        ch_versions
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
