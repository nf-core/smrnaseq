#!/usr/bin/env nextflow

/*
========================================================================================
                    R N A - S E Q    T W O    P O I N T    Z E R O
========================================================================================
 New RNA-Seq Best Practice Analysis Pipeline. Started March 2016.
 @Authors
 Phil Ewels <phil.ewels@scilifelab.se>
 Rickard Hammarén <rickard.hammaren@scilifelab.se>
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow main.nf
 
 Pipeline variables can be configured with the following command line options:
 --genome [GRCh37 | GRCm38]
 --index [path to STAR index]
 --gtf [path to GTF file]
 --files [path to input files]
 
 For example:
 $ nextflow rnaseq.nf --files path/to/data/*fq.gz
----------------------------------------------------------------------------------------
 Pipeline overview:
 - MultiQC
----------------------------------------------------------------------------------------
 GA project GA_14_20 RNA-Seq Pipeline. See planning document:
 https://docs.google.com/document/d/1_I4r-yYLl_nA5SzMKtABjDKxQxHSb5N9FMWyomVSWVU/edit#heading=h.uc2543wvne80
----------------------------------------------------------------------------------------
*/



/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1

// Reference genome index
params.genome = 'GRCh37'
params.index = params.genomes[ params.genome ].star
params.gtf   = params.genomes[ params.genome ].gtf
params.bed12 = params.genomes[ params.genome ].bed12

params.name = "miRNA-Seq Best practice"

// Input files
params.reads = "data/*.fastq.gz"

// Output path
params.out = "$PWD"

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs=file(params.rlocation)

log.info "===================================="
log.info " RNAbp : miRNA-Seq Best Practice v${version}"
log.info "===================================="
log.info "Reads        : ${params.reads}"
log.info "Genome       : ${params.genome}"
log.info "Index        : ${params.index}"
log.info "Annotation   : ${params.gtf}"
log.info "Current home : $HOME"
log.info "Current user : $USER"
log.info "Current path : $PWD"
log.info "R libraries  : ${params.rlocation}"
log.info "Script dir   : $baseDir"
log.info "Working dir  : $workDir"
log.info "Output dir   : ${params.out}"
log.info "===================================="

// Create R library directories if not already existing
nxtflow_libs.mkdirs()

// Set up nextflow objects
index = file(params.index)
gtf   = file(params.gtf)
bed12 = file(params.bed12)

// Validate inputs
if( !index.exists() ) exit 1, "Missing STAR index: ${index}"
if( !gtf.exists() )   exit 2, "Missing GTF annotation: ${gtf}"
if( !bed12.exists() ) exit 2, "Missing BED12 annotation: ${bed12}"

//Setting up a directory to save results to 
results_path = './results'

/*
 * Create a channel for read files 
 */
 
Channel
     .fromPath( params.reads )
     .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
     .map { path ->  
        def prefix = readPrefix(path, params.reads)
        tuple(prefix, path)
     }
     .groupTuple(sort: true)
     .set { read_files }
 
read_files.into  { read_files_fastqc; read_files_trimming;name_for_star }
/*
 * STEP 1 - FastQC
 */

process fastqc {
    tag "reads: $name"

    module 'bioinfo-tools'
    module 'FastQC'

    memory '2 GB'
    time '1h'

    publishDir "$results_path/fastqc"

    input:
    set val(name), file(reads:'*') from read_files_fastqc

    output:
    file '*_fastqc.html' into fastqc_html
    file '*_fastqc.zip' into fastqc_zip

    """
    fastqc -q ${reads}
    """
}


/*
 * STEP 2 - Trim Galore!
 */

process trim_galore {
    tag "reads: $name"

    module 'bioinfo-tools'
    module 'FastQC'
    module 'cutadapt'
    module 'TrimGalore'

    cpus 3
    memory '3 GB'
    time '8h'

    publishDir "$results_path/trim_galore"

    input:
    set val(name), file(reads:'*') from read_files_trimming
    

    output:
    file '*fq.gz' into trimmed_reads
    file '*trimming_report.txt' into results
    script:

    """
    trim_galore --gzip --fastqc_args "-q" $reads
     """
}


/*
 * STEP 3
 */
process bowtie{

    module 'Bioinfo-tools'
    module 'bowtie'

    cpus 4
    memory '16 GB'
    input:
    file (reads:'*') from trimmed_reads
    output:
    file 'out.bam' into bow_out
    """
    bowtie -p 4 -S -q -n 1 -e 80 -l 30 -a -m 5 --best --strata ${reads} out.sam
    samtools view -b -S out.sam | out.bam
    """ 
}


/*
 * STEP  Feature counts
 */


process featureCounts {
    
    module 'bioinfo-tools'
    module 'subread'
    
    memory '4 GB'
    time '2h'
    
    publishDir "$results_path/featureCounts"
    input:
    file bam
    file gtf from gtf
    
    output:
    file '*_gene.featureCounts.txt' into results
    file '*_biotype.featureCounts.txt' into results
    file '*_rRNA_counts.txt' into results
    file '*.summary' into results
    file 'featureCounts.done' into featureCounts_done    
    """
    featureCounts -a $gtf -g gene_id -o ${bam_featurecounts}_gene.featureCounts.txt -p -s 2 $bam_featurecounts
    featureCounts -a $gtf -g gene_biotype -o ${bam_featurecounts}_biotype.featureCounts.txt -p -s 2 $bam_featurecounts
    cut -f 1,7 ${bam_featurecounts}_biotype.featureCounts.txt | sed '1,2d' | grep 'rRNA' > ${bam_featurecounts}_rRNA_counts.txt
    echo done >featureCounts.done
    """
}
