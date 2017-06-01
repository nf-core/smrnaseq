#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
          N G I    S M A L L    R N A - S E Q    B E S T    P R A C T I C E
========================================================================================
 Small-RNA-Seq Best Practice Analysis Pipeline. Started May 2016.
 #### Homepage / Documentation
 https://github.com/SciLifeLab/NGI-smRNAseq
 #### Authors
 Phil Ewels <phil.ewels@scilifelab.se>
 Chuan Wang <chuan.wang@scilifelab.se>
 Rickard Hammarén <rickard.hammaren@scilifelab.se>
----------------------------------------------------------------------------------------
*/

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 1.4

// Configurable variables
params.genome = false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bt2index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.mature = params.genome ? params.genomes[ params.genome ].mature ?: false : false
params.hairpin = params.genome ? params.genomes[ params.genome ].hairpin ?: false : false
params.name = "miRNA-Seq Best practice"
params.saveReference = false
params.reads = "data/*.fastq.gz"
params.outdir = './results'

// Check that we have a mature / hairpin reference
if( !params.mature || !params.hairpin ){
    exit 1, "Missing mature / hairpin reference indexes! Is --genome specified?"
}

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs=file(params.rlocation)
nxtflow_libs.mkdirs()

log.info "==========================================="
log.info " RNAbp : Small RNA-Seq Best Practice v${version}"
log.info "==========================================="
log.info "Reads                : ${params.reads}"
if(params.genome)   log.info "Genome               : ${params.genome}"
if(params.bt2index) log.info "Bowtie2 Index        : ${params.bt2index}"
if(params.gtf)      log.info "Annotation           : ${params.gtf}"
log.info "miRBase mature       : ${params.mature}"
log.info "miRBase hairpin      : ${params.hairpin}"
log.info "Current user         : $USER"
log.info "Current path         : $PWD"
log.info "Script dir           : $baseDir"
log.info "R libraries          : ${params.rlocation}"
log.info "Output dir           : ${params.outdir}"
log.info "Config Profile       : ${workflow.profile}"
log.info "==========================================="

// Set up nextflow objects

// Validate inputs
if( params.mature ){
    mature = file(params.mature)
    if( !mature.exists() ) exit 1, "Mature file not found: ${params.mature}"
}
if( params.hairpin ){
    hairpin = file(params.hairpin)
    if( !hairpin.exists() ) exit 1, "Hairpin file not found: ${params.hairpin}"
}
if( params.gtf ){
    gtf = file(params.gtf)
    if( !gtf.exists() ) exit 1, "GTF file not found: ${params.gtf}"
}
if( params.bt2index ){
    bt2_index = file("${params.bt2index}.1.bt2")
    bt2_indices = Channel.fromPath( "${params.bt2index}*" ).toList()
    if( !bt2_index.exists() ) exit 1, "Reference genome Bowtie 2 not found: ${params.bt2index}"
}
if( !params.gtf || !params.bt2index) {
    log.info "No GTF / Bowtie2 index supplied - host reference genome analysis will be skipped."
}
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"


/*
 * Create a channel for input read files
 */
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { raw_reads_fastqc; raw_reads_trimgalore }


/*
 * PREPROCESSING - Build Bowtie index for mature and hairpin
 */
process makeBowtieIndex {

    publishDir path: { params.saveReference ? "${params.outdir}/bowtie/reference" : params.outdir },
               saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file mature from mature
    file hairpin from hairpin

    output:
    file 'mature_idx.*' into mature_index
    file 'hairpin_idx.*' into hairpin_index

    script:
    """
    fasta_formatter -w 0 -i $mature -o mature_igenome.fa
    fasta_nucleotide_changer -d -i mature_igenome.fa -o mature_idx.fa
    bowtie-build mature_idx.fa mature_idx
    fasta_formatter -w 0 -i $hairpin -o hairpin_igenome.fa
    fasta_nucleotide_changer -d -i hairpin_igenome.fa -o hairpin_idx.fa
    bowtie-build hairpin_idx.fa hairpin_idx
    """
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$reads"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    file reads from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    tag "$reads"
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    file reads from raw_reads_trimgalore

    output:
    file '*.gz' into trimmed_reads_bowtie, trimmed_reads_bowtie2
    file '*trimming_report.txt' into trimgalore_results

    script:
    """
    trim_galore --small_rna --gzip $reads
    """
}


/*
 * STEP 3 - Bowtie miRBase mature miRNA
 */
process bowtie_miRBase_mature {
    tag "$reads"
    publishDir "${params.outdir}/bowtie/miRBase_mature", mode: 'copy', pattern: '*.mature_unmapped.fq.gz'

    input:
    file reads from trimmed_reads_bowtie
    file index from mature_index

    output:
    file '*.mature.bam' into miRBase_mature_bam
    file '*.mature_unmapped.fq.gz' into mature_unmapped_reads

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    bowtie \\
        $index_base \\
        -q <(zcat $reads) \\
        -p 2 \\
        -t \\
        -n 0 \\
        -l 15 \\
        -e 99999 \\
        -k 10 \\
        --best \\
        --chunkmbs 2048 \\
        --un ${prefix}.mature_unmapped.fq \\
        -S \\
        | samtools view -bS - > ${prefix}.mature.bam

    gzip ${prefix}.mature_unmapped.fq
    """
}

/*
 * STEP 4 - Bowtie against miRBase hairpin
 */
process bowtie_miRBase_hairpin {
    tag "$reads"
    publishDir "${params.outdir}/bowtie/miRBase_hairpin", mode: 'copy', pattern: '*.hairpin_unmapped.fq.gz'

    input:
    file reads from mature_unmapped_reads
    file index from hairpin_index

    output:
    file '*.hairpin.bam' into miRBase_hairpin_bam
    file '*.hairpin_unmapped.fq.gz' into hairpin_unmapped_reads

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() - '.mature_unmapped.fq.gz'
    """
    bowtie \\
        $index_base \\
        -p 2 \\
        -t \\
        -n 1 \\
        -l 15 \\
        -e 99999 \\
        -k 10 \\
        --best \\
        --chunkmbs 2048 \\
        -q <(zcat $reads) \\
        --un ${prefix}.hairpin_unmapped.fq \\
        -S \\
        | samtools view -bS - > ${prefix}.hairpin.bam

    gzip ${prefix}.hairpin_unmapped.fq
    """
}


/*
 * STEP 5 - Post-alignment processing for miRBase mature and hairpin
 */
def wrap_mature_and_hairpin = { file ->
    if ( file.contains("mature") ) return "miRBase_mature/$file"
    if ( file.contains("hairpin") ) return "miRBase_hairpin/$file"
}

process miRBasePostAlignment {
    publishDir "${params.outdir}/bowtie", mode: 'copy', saveAs: wrap_mature_and_hairpin

    input:
    file input from miRBase_mature_bam.mix(miRBase_hairpin_bam)

    output:
    file "${input.baseName}.count" into miRBase_counts
    file "${input.baseName}.sorted.bam" into miRBase_bam
    file "${input.baseName}.sorted.bam.bai" into miRBase_bai

    script:
    """
    samtools sort ${input.baseName}.bam -o ${input.baseName}.sorted.bam
    samtools index ${input.baseName}.sorted.bam
    samtools idxstats ${input.baseName}.sorted.bam > ${input.baseName}.count
    """
}


/*
 * STEP 6 - edgeR miRBase feature counts processing
 */
process edgeR_miRBase {
    publishDir "${params.outdir}/edgeR", mode: 'copy', saveAs: wrap_mature_and_hairpin

    input:
    file input_files from miRBase_counts.toSortedList()

    output:
    file '*.{txt,pdf}' into edgeR_miRBase_results

    script:
    """
    edgeR_miRBase.r $params.rlocation $input_files
    """
}


/*
 * STEP 7.1 and 7.2 IF A GENOME SPECIFIED ONLY!
 */
if( params.gtf && params.bt2index) {

    /*
     * STEP 7.1 - Bowtie 2 against reference genome
     */
    process bowtie2 {
        tag "$reads"
        publishDir "${params.outdir}/bowtie2", mode: 'copy'

        input:
        file reads from trimmed_reads_bowtie2
        file index from bt2_index
        file bt2_indices

        output:
        file '*.bowtie2.bam' into bowtie2_bam
        stdout into bowtie2_log

        script:
        index_base = index.toString() - '.1.bt2'
        prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        bowtie2 \\
            -x $index_base \\
            -U $reads \\
            -k 10 \\
            --very-sensitive \\
            -p 8 \\
            -t \\
            | samtools view -bT $index_base - > ${prefix}.bowtie2.bam
        """
    }


    /*
     * STEP 7.2 - NGI-Visualizations of Bowtie 2 alignment statistics
     */
    process ngi_visualizations {
        tag "$bowtie2_bam"
        publishDir "${params.outdir}/bowtie2/ngi_visualizations", mode: 'copy'

        input:
        file gtf from gtf
        file bowtie2_bam

        output:
        file '*.{png,pdf}' into bowtie2_ngi_visualizations

        script:
        // Note! ngi_visualizations needs to be installed!
        // See https://github.com/NationalGenomicsInfrastructure/ngi_visualizations
        """
        #!/usr/bin/env python
        from ngi_visualizations.biotypes import count_biotypes
        count_biotypes.main('$gtf','$bowtie2_bam')
        """
    }

}


/*
 * STEP 8 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file ('fastqc/*') from fastqc_results.flatten().toList()
    file ('trim_galore/*') from trimgalore_results.flatten().toList()
    file ('edgeR/*') from edgeR_miRBase_results.flatten().toList()
    // if( params.gtf && params.bt2index) file ('bowtie2_log/*') from bowtie2_log.flatten().toList()

    output:
    file '*multiqc_report.html' into multiqc_html
    file '*multiqc_data' into multiqc_data

    script:
    """
    multiqc -f .
    """
}
