// Synthetic alignment-only harness: retain the production process names so
// conf/modules.config, rather than test-specific Bowtie arguments, is exercised.
include { BOWTIE_BUILD } from '../../modules/nf-core/bowtie/build/main'
include {
    BOWTIE_ALIGN as BOWTIE_MAP_MATURE
    BOWTIE_ALIGN as BOWTIE_MAP_HAIRPIN
    BOWTIE_ALIGN as BOWTIE_MAP_SEQCLUSTER
    BOWTIE_ALIGN as BOWTIE_MAP_GENOME
} from '../../modules/nf-core/bowtie/align/main'

process ALIGNMENT_SUMMARY {
    container 'community.wave.seqera.io/library/bowtie_htslib_samtools:e1e242368ffcb5d3'
    conda '../../modules/nf-core/bowtie/align/environment.yml'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    path "${meta.id}.tsv"

    script:
    """
    samtools view -F 4 $bam | cut -f 1-3 > ${meta.id}.tsv
    """
}

workflow MIRNA_QUANT {
    take:
    reads
    index

    main:
    BOWTIE_MAP_MATURE(reads.map { meta, fastq -> [meta + [id: 'mature'], fastq] }, index, false)
    BOWTIE_MAP_HAIRPIN(reads.map { meta, fastq -> [meta + [id: 'hairpin'], fastq] }, index, false)
    BOWTIE_MAP_SEQCLUSTER(reads.map { meta, fastq -> [meta + [id: 'seqcluster'], fastq] }, index, false)

    emit:
    bam = BOWTIE_MAP_MATURE.out.bam.mix(BOWTIE_MAP_HAIRPIN.out.bam, BOWTIE_MAP_SEQCLUSTER.out.bam)
}

workflow GENOME_QUANT {
    take:
    reads
    index

    main:
    BOWTIE_MAP_GENOME(reads.map { meta, fastq -> [meta + [id: 'genome'], fastq] }, index, false)

    emit:
    bam = BOWTIE_MAP_GENOME.out.bam
}

workflow NFCORE_SMRNASEQ {
    main:
    reads = channel.value([[id: 'synthetic', single_end: true], file("${projectDir}/reads.fastq")])
    BOWTIE_BUILD(channel.value([[id: 'synthetic'], file("${projectDir}/reference.fa")]))
    MIRNA_QUANT(reads, BOWTIE_BUILD.out.index)
    GENOME_QUANT(reads, BOWTIE_BUILD.out.index)
    ALIGNMENT_SUMMARY(MIRNA_QUANT.out.bam.mix(GENOME_QUANT.out.bam))
}

workflow {
    NFCORE_SMRNASEQ()
}
