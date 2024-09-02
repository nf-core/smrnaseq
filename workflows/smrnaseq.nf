/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// nf-core modules
include { CAT_FASTQ                        } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                           } from '../modules/nf-core/fastqc/main'
include { FASTP as FASTP_LENGTH_FILTER     } from '../modules/nf-core/fastp'
include { MULTIQC                          } from '../modules/nf-core/multiqc/main'
include { UMICOLLAPSE as UMICOLLAPSE_FASTQ } from '../modules/nf-core/umicollapse/main'
include { UMITOOLS_EXTRACT                 } from '../modules/nf-core/umitools/extract/main'
// nf-core subworkflows
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp'
include { paramsSummaryMultiqc             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// local subworkflows
include { CONTAMINANT_FILTER               } from '../subworkflows/local/contaminant_filter'
include { GENOME_QUANT                     } from '../subworkflows/local/genome_quant'
include { MIRNA_QUANT                      } from '../subworkflows/local/mirna_quant'
include { MIRDEEP2                         } from '../subworkflows/local/mirdeep2'
include { MIRTRACE                         } from '../subworkflows/local/mirtrace'
// plugins
include { paramsSummaryMap                 } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_fastp_adapters                     = Channel.fromPath(params.fastp_known_mirna_adapters, checkIfExists: true).collect() // collect to consume for all incoming samples to FASTP

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCORE_SMRNASEQ {

    take:
    has_fasta            // boolean
    has_mirtrace_species // boolean
    ch_samplesheet       // channel: sample fastqs parsed from --input
    ch_mirna_adapters    // channel: [ val(string) ]
    ch_mirtrace_species  // channel: [ val(string) ]
    ch_reference_mature  // channel: [ val(meta), path(fasta) ]
    ch_reference_hairpin // channel: [ val(meta), path(fasta) ]
    ch_mirna_gtf         // channel: [ path(GTF) ]
    ch_fasta             // channel: [ val(meta), path(fasta) ]
    ch_bowtie_index      // channel: [genome.1.ebwt, genome.2.ebwt, genome.3.ebwt, genome.4.ebwt, genome.rev.1.ebwt, genome.rev.2.ebwt ]
    ch_versions          // channel: [ path(versions.yml) ]

    main:
    //
    // Create separate channels for samples that have single/multiple FastQ files to merge
    //
    ch_fastq = ch_samplesheet
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    ch_cat_fastq = CAT_FASTQ.out.reads.mix(ch_fastq.single)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters & dedup UMIs if necessary / desired by the user
    //
    if ( params.skip_fastp && params.skip_fastqc ) {
        exit 1, "At least one of skip_fastp or skip_fastqc must be false"
    }

    FASTQ_FASTQC_UMITOOLS_FASTP (
        ch_cat_fastq,
        params.skip_fastqc,
        params.with_umi,
        params.skip_umi_extract_before_dedup,
        Channel.value(params.umi_discard_read).ifEmpty([]),
        params.skip_fastp,
        ch_mirna_adapters,
        params.save_trimmed_fail,
        params.save_merged,
        Channel.value(params.min_trimmed_reads).ifEmpty([])
    )
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)

    ch_reads_for_mirna = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads

    // UMI Dedup for fastq input
    // This involves running on the sequencing adapter trimmed remnants of the entire reads
    // consisting of sequence + common sequence "miRNA adapter" + UMI
    // once collapsing happened, we will use umitools extract to get rid of the common miRNA sequence + the UMI to have only plain collapsed reads without any other clutter
    if (params.with_umi) {
        ch_fastq = Channel.value('fastq')
        ch_input_for_collapse = ch_reads_for_mirna.map{ meta, reads -> [meta, reads, []]} //Needs to be done to add a []
        UMICOLLAPSE_FASTQ(ch_input_for_collapse, ch_fastq)
        ch_versions = ch_versions.mix(UMICOLLAPSE_FASTQ.out.versions)
        UMITOOLS_EXTRACT(UMICOLLAPSE_FASTQ.out.fastq)

        // Filter out sequences smaller than params.fastp_min_length
        FASTP_LENGTH_FILTER (
            UMITOOLS_EXTRACT.out.reads,
            ch_mirna_adapters,
            false,
            params.save_trimmed_fail,
            params.save_merged
        )
        ch_versions = ch_versions.mix(FASTP_LENGTH_FILTER.out.versions)

        ch_reads_for_mirna = FASTP_LENGTH_FILTER.out.reads
    }

    //
    // MODULE: mirtrace QC
    //

    // Define the main adapter sequence channel
    ch_adapter_seq = FASTQ_FASTQC_UMITOOLS_FASTP.out.adapter_seq

    // Define a fallback channel with the default value "custom"
    ch_fallback_adapter_seq = ch_reads_for_mirna.map { meta, reads -> [meta, 'custom'] }

    // Change to fallback channel if ch_adapter_seq is empty
    ch_adapter_seq = ch_adapter_seq ? ch_fallback_adapter_seq : ch_adapter_seq

    // Now join the adapter sequence channel with the reads channel
    ch_mirtrace_inputs = ch_adapter_seq
        .join(ch_reads_for_mirna)
        .map { meta, adapter_seq, reads -> [adapter_seq, meta.id, reads] }
        .groupTuple()
        .map { adapter_seq, ids, reads_list -> [adapter_seq, ids, reads_list.flatten()] }

    //
    // SUBWORKFLOW: MIRTRACE
    //
    if (has_mirtrace_species){
        MIRTRACE(ch_mirtrace_inputs, ch_mirtrace_species)
        ch_versions = ch_versions.mix(MIRTRACE.out.versions)
    } else {
        log.warn "The parameter --mirtrace_species is absent. MIRTRACE quantification skipped."
    }

    //
    // SUBWORKFLOW: remove contaminants from reads
    //
    contamination_stats = Channel.empty()
    if (params.filter_contamination){
        CONTAMINANT_FILTER (
            ch_reference_hairpin.map{meta,file -> file},
            Channel.value(params.rrna).ifEmpty([]),
            Channel.value(params.trna).ifEmpty([]),
            Channel.value(params.cdna).ifEmpty([]),
            Channel.value(params.ncrna).ifEmpty([]),
            Channel.value(params.pirna).ifEmpty([]),
            Channel.value(params.other_contamination).ifEmpty([]),
            ch_reads_for_mirna
        )
        contamination_stats = CONTAMINANT_FILTER.out.filter_stats
        ch_versions = ch_versions.mix(CONTAMINANT_FILTER.out.versions)
        ch_reads_for_mirna = CONTAMINANT_FILTER.out.filtered_reads
    }

    //MIRNA_QUANT process should still run even if mirtrace_species is null, when mirgendb is true
    MIRNA_QUANT (
        ch_reference_mature,
        ch_reference_hairpin,
        ch_mirna_gtf,
        ch_reads_for_mirna,
        ch_mirtrace_species
    )
    ch_versions = ch_versions.mix(MIRNA_QUANT.out.versions)

    //
    // GENOME
    //
    genome_stats = Channel.empty()
    if (has_fasta){
        GENOME_QUANT (
            ch_bowtie_index,
            ch_fasta,
            MIRNA_QUANT.out.unmapped
        )
        genome_stats = GENOME_QUANT.out.stats
        ch_versions = ch_versions.mix(GENOME_QUANT.out.versions)

        hairpin_clean = MIRNA_QUANT.out.fasta_hairpin.map { it -> it[1] }
        mature_clean  = MIRNA_QUANT.out.fasta_mature.map { it -> it[1] }

        if (!params.skip_mirdeep) {
            MIRDEEP2 (
                ch_reads_for_mirna,
                GENOME_QUANT.out.fasta,
                GENOME_QUANT.out.index.collect(),
                hairpin_clean,
                mature_clean
            )
            ch_versions = ch_versions.mix(MIRDEEP2.out.versions)
        }
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
        // .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_smrnaseq_software_mqc_versions.yml', sort: true, newLine: true)
        // .set {ch_collated_versions}

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = Channel.empty()
    if (!params.skip_multiqc) {
        // summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        // ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_config        = Channel.fromPath(
            "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ?
            Channel.fromPath(params.multiqc_config, checkIfExists: true) :
            Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo ?
            Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
            Channel.empty()

        summary_params      = paramsSummaryMap(
            workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            Channel.fromPath(params.multiqc_methods_description, checkIfExists: true) :
            Channel.fromPath("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        // ch_methods_description                = Channel.value(
        //     methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        // ch_multiqc_files = ch_multiqc_files.mix(
        //     ch_methods_description.collectFile(
        //         name: 'methods_description_mqc.yaml',
        //         sort: true
        //     )
        // )

        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json.collect{it[1]}.ifEmpty([]))
        if(params.with_umi) {
            ch_multiqc_files = ch_multiqc_files.mix(UMICOLLAPSE_FASTQ.out.log.collect{it[1]}.ifEmpty([]))
        }
        ch_multiqc_files = ch_multiqc_files.mix(contamination_stats.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(genome_stats.collect({it[1]}).ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.mature_stats.collect({it[1]}).ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.hairpin_stats.collect({it[1]}).ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MIRNA_QUANT.out.mirtop_logs.collect().ifEmpty([]))
        if (has_mirtrace_species){
            ch_multiqc_files = ch_multiqc_files.mix(MIRTRACE.out.results.collect().ifEmpty([]))
        }

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
        ch_multiqc_report = MULTIQC.out.report

    }

    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
