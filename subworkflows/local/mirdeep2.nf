//
// Quantify mirna with bowtie and mirtop
//

params.options = [:]

include { MIRDEEP2_PIGZ } from '../../modules/local/mirdeep2'
include { MIRDEEP2_MAPPER } from '../../modules/local/mirdeep2'
include { MIRDEEP2_RUN } from '../../modules/local/mirdeep2'

workflow MIRDEEP2 {
    take:
    reads        // channel: [ val(meta), [ reads ] ]
    fasta
    indeces
    hairpin
    mature

    main:

    MIRDEEP2_PIGZ ( reads )

    MIRDEEP2_MAPPER ( MIRDEEP2_PIGZ.out.reads, indeces )

    MIRDEEP2_RUN ( fasta, MIRDEEP2_MAPPER.out.mirdeep2_inputs, hairpin, mature )

    emit:
        version   = MIRDEEP2_MAPPER.out.version


}