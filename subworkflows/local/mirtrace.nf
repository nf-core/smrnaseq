//
// Quantify mirna with bowtie and mirtop
//

include { MIRTRACE_RUN } from '../../modules/local/mirtrace'

workflow MIRTRACE {
    take:
    reads      // channel: [ val(adapterseq), [ val(ids) ], [ path(reads) ] ]

    main:

    //Staging the files as path() but then getting the filenames for the config file that mirtrace needs
    //Directly using val(reads) as in previous versions is not reliable as staging between work directories is not 100% reliable if not explicitly defined via nextflow itself
    //mirtrace is a bit peculiar in parsing these config files, so looked it up in the source how its done. this way should work
    ch_mirtrace_config =
    reads.map { adapter, ids, reads -> [adapter, ids,reads]}
    .transpose()
    .collectFile { adapter, id, path -> "./${path.getFileName().toString()},${id},${adapter},${params.phred_offset}\n" } // operations need a channel, so, should be outside the module

    MIRTRACE_RUN (
        reads,
        ch_mirtrace_config
    )

    emit:
    results    = MIRTRACE_RUN.out.mirtrace
    versions   = MIRTRACE_RUN.out.versions
}
