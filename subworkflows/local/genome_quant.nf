//
// Quantify mirna with bowtie and mirtop
//

params.options = [:]

include { INDEX_GENOME } from '../../modules/local/bowtie_genome'


workflow GENOME_QUANT {
    take:
    fasta
    bt_index
    mirtrace_species
    gtf        // channle: GTF file
    mature     // channel: fasta file
    hairpin    // channel: fasta file
    reads      // channel: [ val(meta), [ reads ] ]

    main:
    if (!bt_index && !params.skip_mirdeep){
        bt_indices = INDEX_GENOME( fasta ).out.bt_indeces
    } else if (bt_index) {
        bt_indices = Channel.fromPath("${bt_index}**ebwt", checkIfExists: true).ifEmpty { exit 1, "Bowtie1 index directory not found: ${bt_dir}" }

    } else {
        bt_indeces = Channel.empty()
    }
    //emit:
    // fasta_mirna   = mirna_formatted


}
