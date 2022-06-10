//
// Quantify mirna with bowtie and mirtop
//

include { MIRTRACE_RUN } from '../../modules/local/mirtrace'

workflow MIRTRACE {
    take:
    reads      // channel: [ val(meta), [ reads ] ]

    main:
    reads
        .map { it[1] }
        .flatten()
        .dump(tag:'mirtrace')
        .set { all_reads }

    MIRTRACE_RUN ( all_reads.collect() )

    emit:
    results    = MIRTRACE_RUN.out.mirtrace
    versions   = MIRTRACE_RUN.out.versions
}
