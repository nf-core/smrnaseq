//
// Write one miRTrace config file per sample: <reads>,<display name>,<adapter>,<PHRED offset>
//

include { mirtraceAdapter } from '../utils_nfcore_smrnaseq_pipeline'
include { mirtraceConfigLine } from '../utils_nfcore_smrnaseq_pipeline'

workflow MIRTRACE_CONFIG {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]
    ch_fastp_adapter // channel: [ val(meta), val(adapter) ] adapter detected by fastp; empty channel if fastp was skipped
    ch_three_prime_adapter // channel: [ val(string) ]
    ch_phred_offset // channel: [ val(string) ]

    main:
    // Match samples by id: UMI processing can change the meta map before fastp. Samples without a fastp
    // result (e.g. --skip_fastp) are kept by the remainder join and get the fallback adapter.
    ch_reads_adapter = ch_reads
        .map { meta, reads -> [meta.id, meta, reads] }
        .join(ch_fastp_adapter.map { meta, adapter -> [meta.id, adapter] }, remainder: true)
        .filter { _id, meta, _reads, _fastp_adapter -> meta != null }
        .combine(ch_three_prime_adapter)
        .map { _id, meta, reads, fastp_adapter, three_prime_adapter ->
            [meta, reads, mirtraceAdapter(three_prime_adapter, fastp_adapter)]
        }

    // miRTrace cannot auto-detect adapters: warn once about the samples that run without adapter trimming
    ch_reads_adapter
        .filter { _meta, _reads, adapter -> adapter == 'auto-detect' }
        .map { meta, _reads, _adapter -> meta.id }
        .collect(sort: true)
        .subscribe { ids ->
            log.warn("miRTrace cannot auto-detect adapters, and fastp reported no adapter (or was skipped) for: ${ids.join(', ')}. miRTrace runs without adapter trimming for these samples. Set '--three_prime_adapter' to the adapter sequence for adapter-trimmed miRTrace QC.")
        }

    ch_config = ch_reads_adapter
        .transpose()
        .combine(ch_phred_offset)
        .collectFile { meta, reads, adapter, phred ->
            ["${meta.id}.data", mirtraceConfigLine(reads.getFileName().toString(), meta.id, adapter, phred)]
        }
        .map { config_file -> [[id: config_file.getBaseName()], config_file] }

    emit:
    config = ch_config // channel: [ val(meta), path(config) ]
}
