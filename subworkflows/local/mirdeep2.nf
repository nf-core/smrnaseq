//
// Quantify mirna with bowtie and mirtop
//

include { PIGZ_UNCOMPRESS } from '../../modules/nf-core/pigz/uncompress/main'
include { MIRDEEP2_MAPPER } from '../../modules/local/mirdeep2_mapper'
include { MIRDEEP2_RUN    } from '../../modules/local/mirdeep2_run'

workflow MIRDEEP2 {
    take:
    ch_reads_for_mirna // channel: [ val(meta), [ reads ] ]
    ch_fasta           // channel: [ val(meta), path(fasta) ]
    ch_bowtie_index    // channel: [ genome.1.ebwt, genome.2.ebwt, genome.3.ebwt, genome.4.ebwt, genome.rev.1.ebwt, genome.rev.2.ebwt ]
    ch_hairpin_clean   // channel: [ path(hairpin.fa) ]
    ch_mature_clean    // channel: [ path(mature.fa)  ]

    main:
    ch_versions = Channel.empty()

    PIGZ_UNCOMPRESS ( ch_reads_for_mirna )
    ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS.out.versions.first())

    MIRDEEP2_MAPPER ( PIGZ_UNCOMPRESS.out.file, ch_bowtie_index )
    ch_versions = ch_versions.mix(MIRDEEP2_MAPPER.out.versions.first())

    MIRDEEP2_RUN ( ch_fasta.map{meta,file->file}, MIRDEEP2_MAPPER.out.mirdeep2_inputs, ch_hairpin_clean, ch_mature_clean )
    ch_versions = ch_versions.mix(MIRDEEP2_RUN.out.versions.first())

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
