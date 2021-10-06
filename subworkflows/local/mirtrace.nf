//
// Quantify mirna with bowtie and mirtop
//

params.options = [:]

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
    version   = MIRTRACE_RUN.out.version


}
