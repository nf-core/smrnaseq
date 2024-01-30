//
// Quantify mirna with bowtie and mirtop
//

include { MIRTRACE_RUN } from '../../modules/local/mirtrace'

workflow MIRTRACE {
    take:
    reads      // channel: [ val(adapterseq), [ val(ids) ], [ path(reads) ] ]

    main:
    ch_versions = Channel.empty()
    reads | MIRTRACE_RUN
    ch_versions.mix(MIRTRACE_RUN.out.versions)

    emit:
    results    = MIRTRACE_RUN.out.mirtrace
    versions   = ch_versions
}
