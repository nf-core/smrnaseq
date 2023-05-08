//
// Quantify mirna with bowtie and mirtop
//

include { MIRTRACE_RUN } from '../../modules/local/mirtrace'

workflow MIRTRACE {
    take:
    reads      // channel: [ val(adapterseq), [ val(ids) ], [ path(reads) ] ]

    main:
    reads
    | map { adapterseq, ids, read_collection -> [adapterseq, ids, read_collection*.first()]}
    | MIRTRACE_RUN

    emit:
    results    = MIRTRACE_RUN.out.mirtrace
    versions   = MIRTRACE_RUN.out.versions
}
