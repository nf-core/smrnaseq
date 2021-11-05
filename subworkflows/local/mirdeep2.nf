//
// Quantify mirna with bowtie and mirtop
//

params.options = [:]

include { MIRDEEP2_PIGZ } from '../../modules/local/mirdeep2_prepare'
include { MIRDEEP2_MAPPER } from '../../modules/local/mirdeep2_mapper'
include { MIRDEEP2_RUN } from '../../modules/local/mirdeep2_run'

workflow MIRDEEP2 {
    take:
    reads        // channel: [ val(meta), [ reads ] ]
    fasta
    indices
    hairpin
    mature

    main:

    MIRDEEP2_PIGZ ( reads )

    MIRDEEP2_MAPPER ( MIRDEEP2_PIGZ.out.reads, indices )

    MIRDEEP2_RUN ( fasta, MIRDEEP2_MAPPER.out.mirdeep2_inputs, hairpin, mature )

    emit:
    versions_prepare = MIRDEEP2_PIGZ.out.versions
    versions_mapper  = MIRDEEP2_MAPPER.out.versions
    versions_run     = MIRDEEP2_RUN.out.versions

}
