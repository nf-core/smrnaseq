#!/usr/bin/env nextflow

/*
========================================================================================
                                         miRNA 
========================================================================================
 miRNA  Analysis Pipeline. Started May 2016.
 @Authors
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
params.bowtie = params.genomes[ params.genome ].bowtie
params.name = "miRNA-Seq Best practice"

// Input files
params.input = "data/*.fastq.gz"

// Output path
params.out = "$PWD"

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs=file(params.rlocation)

log.info "===================================="
log.info " RNAbp : miRNA-Seq Best Practice v${version}"
log.info "===================================="
log.info "Reads        : ${params.input}"
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
bowtie=file(params.bowtie)

// Validate inputs
if( !index.exists() ) exit 1, "Missing STAR index: ${index}"
if( !gtf.exists() )   exit 2, "Missing GTF annotation: ${gtf}"
if( !bed12.exists() ) exit 2, "Missing BED12 annotation: ${bed12}"
//if( !bowtie.exists() ) exit 2, "Missing  annotation: ${bowtie}"

//Setting up a directory to save results to 
results_path = './results'

dataset = Channel.fromPath(params.input)
dataset.into { fastQC_input; trimgalore_input }  

/*
 * STEP 1 - FastQC
 */

process fastqc {
    tag "reads: $reads"

    module 'bioinfo-tools'
    module 'FastQC'

    memory '2 GB'
    time '1h'

    publishDir "$results_path/fastqc"

    input:
    file reads from fastQC_input
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
    tag "reads: $reads"

    module 'bioinfo-tools'
    module 'FastQC'
    module 'cutadapt'
    module 'TrimGalore'

    cpus 3
    memory '3 GB'
    time '8h'

    publishDir "$results_path/trim_galore"

    input:
    file(reads) from trimgalore_input

    output:
    file '*fq.gz' into trimmed_reads, trimmed_reads_miRdeep2
    file '*trimming_report.txt' 
    script:

    """
    trim_galore --gzip --fastqc_args "-q" ${reads}
    """
}

/*
 * STEP 3 - Bowtie
*/

process bowtie{
    
    module 'bioinfo-tools'
    module 'bowtie'
    module 'samtools'
    
    memory '2 GB'
    time '6h'

    input:
    file(reads:'*') from trimmed_reads 
    
    output:
    file '*.bam' into bam 
   
    """
    f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1}
    prefix=\$f
    bowtie -p 2 -t -n 0 -l 15 -e 99999 -k 200 --best -S --chunkmbs 2048 $bowtie $reads | samtools view -bS > \$prefix.Aligned.bam 
    """
}

/*
 * STEP 4
 */
 
//Under development
process miRdeep2{

    module 'Bioinfo-tools'
    module 'mirdeep2'

    cpus 4
    memory '16 GB'
    input:
    //file (reads:'*') from trimmed_reads_miRdeep2
    output:
    file '*.reads' into mapper_output
    """
    mapper.pl -e $reads -p params.genome -s ${reads}.reads -t ${reads}.mapping  
    """ 
}
*/

/*
 * STEP 5 -  Feature counts
 */


process featureCounts {
    
    module 'bioinfo-tools'
    module 'subread'
    
    memory '4 GB'
    time '2h'
    
    publishDir "$results_path/featureCounts"
    input:
    file '*.bam' from bam
    file gtf from gtf
    
    output:
    file '*_gene.featureCounts.txt' 
    file '*_biotype.featureCounts.txt' 
    file '*_rRNA_counts.txt' 
    file '*.summary' 
    file 'featureCounts.done' into featureCounts_done    
    """
    featureCounts -a $gtf -g gene_id -o ${bam}_gene.featureCounts.txt -p -s 2 $bam_featurecounts
    featureCounts -a $gtf -g gene_biotype -o ${bam}_biotype.featureCounts.txt -p -s 2 $bam_featurecounts
    cut -f 1,7 ${bam}_biotype.featureCounts.txt | sed '1,2d' | grep 'rRNA' > ${bam}_rRNA_counts.txt
    echo done >featureCounts.done
    """
}

