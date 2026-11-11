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
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.three_prime_adapter,
        params.phred_offset,
        params.help,
        params.help_full,
        params.show_hidden
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
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.three_prime_adapter,
        PIPELINE_INITIALISATION.out.phred_offset
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
        NFCORE_SMRNASEQ.out.multiqc_report
    )

    publish:
    bowtie_index_untar   = PREPARE_GENOME.out.bowtie_index_untar
    bowtie_index_genome  = PREPARE_GENOME.out.bowtie_index_built.mix(PREPARE_GENOME.out.fasta)
    fastqc_raw           = NFCORE_SMRNASEQ.out.fastqc_raw
    fastqc_trim          = NFCORE_SMRNASEQ.out.fastqc_trim
    fastp_reports        = NFCORE_SMRNASEQ.out.fastp_reports
    fastp_logs           = NFCORE_SMRNASEQ.out.fastp_logs
    fastp_reads_fail     = NFCORE_SMRNASEQ.out.fastp_reads_fail
    fastp_reads          = NFCORE_SMRNASEQ.out.fastp_reads
    umi_extract_logs     = NFCORE_SMRNASEQ.out.umi_extract_logs
    umi_extract_reads    = NFCORE_SMRNASEQ.out.umi_extract_reads
    umi_dedup            = NFCORE_SMRNASEQ.out.umi_dedup
    umi_fastp_reports    = NFCORE_SMRNASEQ.out.umi_fastp_reports
    umi_fastp_logs       = NFCORE_SMRNASEQ.out.umi_fastp_logs
    umi_fastp_reads_fail = NFCORE_SMRNASEQ.out.umi_fastp_reads_fail
    mirtrace             = NFCORE_SMRNASEQ.out.mirtrace
    contaminant_filter   = NFCORE_SMRNASEQ.out.contaminant_filter
    mirna_reference      = NFCORE_SMRNASEQ.out.mirna_reference
    mirna_index_mature   = NFCORE_SMRNASEQ.out.mirna_index_mature
    mirna_index_hairpin  = NFCORE_SMRNASEQ.out.mirna_index_hairpin
    mirna_bam_mature     = NFCORE_SMRNASEQ.out.mirna_bam_mature
    mirna_bam_hairpin    = NFCORE_SMRNASEQ.out.mirna_bam_hairpin
    mirna_bam_seqcluster = NFCORE_SMRNASEQ.out.mirna_bam_seqcluster
    mirna_seqcluster     = NFCORE_SMRNASEQ.out.mirna_seqcluster
    mirna_edger_qc       = NFCORE_SMRNASEQ.out.mirna_edger_qc
    mirna_mirtop         = NFCORE_SMRNASEQ.out.mirna_mirtop
    genome_quant_stats   = NFCORE_SMRNASEQ.out.genome_quant_stats
    genome_quant_bam     = NFCORE_SMRNASEQ.out.genome_quant_bam
    mirdeep2             = NFCORE_SMRNASEQ.out.mirdeep2
    multiqc              = NFCORE_SMRNASEQ.out.multiqc
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Published relative to `outputDir` (= params.outdir) using `workflow.output.mode`
    (= params.publish_dir_mode), see nextflow.config. The layout matches the one
    documented in docs/output.md.
----------------------------------------------------------------------------------------
*/

output {
    // Reference preparation
    bowtie_index_untar {
        path 'untar'
    }
    bowtie_index_genome {
        path 'bowtie_index/genome'
        enabled params.save_intermediates
    }

    // Read QC and trimming
    fastqc_raw {
        path 'fastqc/raw'
        enabled params.save_intermediates
    }
    fastqc_trim {
        path 'fastqc/trimmed'
        enabled params.save_intermediates
    }
    fastp_reports {
        path 'fastp/on_raw'
    }
    fastp_logs {
        path 'fastp/on_raw/log'
    }
    fastp_reads_fail {
        path 'fastp/on_raw'
        enabled params.save_trimmed_fail
    }
    fastp_reads {
        path 'fastp/fastq'
        enabled params.save_merged
    }

    // UMI handling
    umi_extract_logs {
        path 'umi_dedup/fastq_extracted_umi'
    }
    umi_extract_reads {
        path 'umi_dedup/fastq_extracted_umi'
        enabled params.save_umi_intermeds
    }
    umi_dedup {
        path 'umi_dedup/bam_deduplicated'
        enabled params.save_umi_intermeds
    }
    umi_fastp_reports {
        path 'fastp/on_umi_extracted'
    }
    umi_fastp_logs {
        path 'fastp/on_umi_extracted/log'
    }
    umi_fastp_reads_fail {
        path 'fastp/on_umi_extracted'
        enabled params.save_trimmed_fail
    }

    // miRNA QC
    mirtrace {
        path { meta, _files -> "mirtrace/${meta.id}" }
    }

    // Contamination filtering
    contaminant_filter {
        path 'contaminant_filter/filter'
        enabled params.save_intermediates
    }

    // miRNA quantification
    mirna_reference {
        path 'mirna_quant/reference'
        enabled params.save_intermediates
    }
    mirna_index_mature {
        path 'bowtie_index/mirna_mature'
        enabled params.save_intermediates
    }
    mirna_index_hairpin {
        path 'bowtie_index/mirna_hairpin'
        enabled params.save_intermediates
    }
    mirna_bam_mature {
        path 'mirna_quant/bam/mature'
        enabled params.save_intermediates
    }
    mirna_bam_hairpin {
        path 'mirna_quant/bam/hairpin'
        enabled params.save_intermediates
    }
    mirna_bam_seqcluster {
        path 'mirna_quant/bam/seqcluster'
        enabled params.save_intermediates
    }
    mirna_seqcluster {
        path 'mirna_quant/seqcluster'
        enabled params.save_intermediates
    }
    mirna_edger_qc {
        path 'mirna_quant/edger_qc'
    }
    mirna_mirtop {
        path 'mirna_quant/mirtop'
    }

    // Genome quantification
    genome_quant_stats {
        path 'genome_quant/bam'
    }
    genome_quant_bam {
        path 'genome_quant/bam'
        enabled params.save_intermediates
    }

    // Novel miRNA discovery
    mirdeep2 {
        path 'mirdeep2'
    }

    // Reports
    multiqc {
        path 'multiqc'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
