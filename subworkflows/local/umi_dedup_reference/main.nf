//
// Reference for alignment-based UMI deduplication (`--umi_dedup_method fgumi`, experimental):
// mature miRNA "windows" (mature sequence plus flanking hairpin nucleotides, see UMI_DEDUP_WINDOWS)
// of the selected species, and their Bowtie 1 index.
//
// UMI deduplication runs before MIRNA_QUANT, so the mature and hairpin FASTA files are parsed and
// formatted here in the same way as in MIRNA_QUANT (same species filter, same MirGeneDB star FASTA).
//

include { PARSE_FASTA_MIRNA  as PARSE_MATURE_UMI_DEDUP
          PARSE_FASTA_MIRNA  as PARSE_HAIRPIN_UMI_DEDUP  } from '../../../modules/local/parse_fasta_mirna/main'
include { FORMAT_FASTA_MIRNA as FORMAT_MATURE_UMI_DEDUP
          FORMAT_FASTA_MIRNA as FORMAT_HAIRPIN_UMI_DEDUP } from '../../../modules/local/format_fasta_mirna/main'
include { UMI_DEDUP_WINDOWS                              } from '../../../modules/local/umi_dedup_windows/main'
include { BOWTIE_BUILD       as INDEX_UMI_DEDUP          } from '../../../modules/nf-core/bowtie/build/main'

workflow UMI_DEDUP_REFERENCE {
    take:
    ch_reference_mature  // channel: [ val(meta), path(fasta) ]
    ch_reference_hairpin // channel: [ val(meta), path(fasta) ]
    ch_mirtrace_species  // channel: [ val(string) ]

    main:
    // Same species and star FASTA as in MIRNA_QUANT
    ch_parse_species_input = params.mirgenedb ? channel.value(params.mirgenedb_species) : (params.mirna_species ? channel.value(params.mirna_species) : ch_mirtrace_species)
    val_mature_star = params.mirgenedb && params.mirgenedb_star ? file(params.mirgenedb_star, checkIfExists: true) : null

    PARSE_MATURE_UMI_DEDUP ( ch_reference_mature, ch_parse_species_input, val_mature_star )
    FORMAT_MATURE_UMI_DEDUP ( PARSE_MATURE_UMI_DEDUP.out.parsed_fasta )

    PARSE_HAIRPIN_UMI_DEDUP ( ch_reference_hairpin, ch_parse_species_input, null )
    FORMAT_HAIRPIN_UMI_DEDUP ( PARSE_HAIRPIN_UMI_DEDUP.out.parsed_fasta )

    UMI_DEDUP_WINDOWS (
        FORMAT_MATURE_UMI_DEDUP.out.formatted_fasta,
        FORMAT_HAIRPIN_UMI_DEDUP.out.formatted_fasta
    )

    INDEX_UMI_DEDUP ( UMI_DEDUP_WINDOWS.out.fasta.map { _meta, fasta -> [ [ id: 'umi_dedup_windows' ], fasta ] } )

    emit:
    index           = INDEX_UMI_DEDUP.out.index           // channel: [ val(meta), path(bowtie) ]
    fasta           = UMI_DEDUP_WINDOWS.out.fasta         // channel: [ val(meta), path(fasta) ]
    windows         = UMI_DEDUP_WINDOWS.out.windows       // channel: [ val(meta), path(tsv) ]
    versions_tuples = PARSE_MATURE_UMI_DEDUP.out.versions.mix( // channel: [ process, tool, version ]
        PARSE_HAIRPIN_UMI_DEDUP.out.versions,
        FORMAT_MATURE_UMI_DEDUP.out.versions,
        FORMAT_HAIRPIN_UMI_DEDUP.out.versions
    )
}
