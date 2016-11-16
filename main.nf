#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
          N G I    S M A L L    R N A - S E Q    B E S T    P R A C T I C E
========================================================================================
 small RNA-seq  Analysis Pipeline. Started May 2016.
 @Authors
 Phil Ewels <phil.ewels@scilifelab.se>
 Chuan Wang <chuan.wang@scilifelab.se>
 Rickard Hammarén <rickard.hammaren@scilifelab.se>
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow main.nf

 Pipeline variables can be configured with the following command line options:
 --genome [ID]
 --index [path to bowtie2 index] (alternative  to --genome)
 --gtf [path to GTF file]
 --reads [path to input files]

 For example:
 $ nextflow main.nf --reads '*.fastq.gz'

 NOTE! Paired-end data is not supported by this pipeline.
 For paired-end data, use Read 1 only:
 $ nextflow main.nf --reads '*.R1.fastq.gz'

 NOTE! With the option --genome 'ALL', the entire dataset of mature miRNAs and hairpins
 in miRBase will be used as reference regardless of species. Meanwhile the alignment against
 host reference genome will be skipped.

----------------------------------------------------------------------------------------
 Pipeline overview:
 - 1:   FastQC for raw sequencing reads quility control
 - 2:   Trim Galore! for adapter trimming
 - 3:   Bowtie 1 alignment against miRBase mature miRNA
 - 4:   Bowtie 1 alignment against miRBase hairpin for the unaligned reads in step 3
 - 5:   Post-alignment processing for miRBase mature and hairpin
 - 6:   edgeR analysis on miRBase counts for both mature miRNA and miRNA precursors (hairpins)
          - TMM normalization and a table of top expression features
          - MDS plot clustering samples
          - Heatmap of sample similarities
 - 7.1: Bowtie 2 alignment against host reference genome (for specified genome only)
 - 7.2: NGI-Visualization of Bowtie 2 alignment statistics (for specified genome only)
 - 8:   MultiQC
----------------------------------------------------------------------------------------
*/


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 1.3

// Reference genome index

if (params.genome != "ALL") {
  params.genome = 'GRCh37'
  params.gtf   = params.genomes[ params.genome ].gtf
  params.bed12 = params.genomes[ params.genome ].bed12
  params.index = params.genomes[ params.genome ].bowtie2
  params.mature  = params.miRBase_mature[ params.genome ].bowtie
  params.hairpin = params.miRBase_hairpin[ params.genome ].bowtie
}
else{
  params.mature  = params.miRBase_mature[ params.genome ].bowtie
  params.hairpin = params.miRBase_hairpin[ params.genome ].bowtie
}

params.name = "miRNA-Seq Best practice"
params.reads = "data/*.fastq.gz"
params.outdir = './results'

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs=file(params.rlocation)
nxtflow_libs.mkdirs()

log.info "===================================="
log.info " RNAbp : Small RNA-Seq Best Practice v${version}"
log.info "===================================="
log.info "Reads                : ${params.reads}"
log.info "Genome               : ${params.genome}"
log.info "Bowtie2 Index        : ${params.index}"
log.info "Annotation           : ${params.gtf}"
log.info "miRBase mature       : ${params.mature}"
log.info "miRBase hairpin      : ${params.hairpin}"
log.info "Current home         : $HOME"
log.info "Current user         : $USER"
log.info "Current path         : $PWD"
log.info "R libraries          : ${params.rlocation}"
log.info "Script dir           : $baseDir"
log.info "Working dir          : $workDir"
log.info "Output dir           : ${params.outdir}"
log.info "Config Profile : ${workflow.profile}"
log.info "===================================="

// Set up nextflow objects

if (params.genome != "ALL") {
  gtf = file(params.gtf)
  bed12 = file(params.bed12)
  index = file(params.index)
  mature = file(params.mature)
  hairpin = file(params.hairpin)
}
else{
  mature = file(params.mature)
  hairpin = file(params.hairpin)
}

// Validate inputs
if( !mature.exists() ) exit 1, "Missing Bowtie 1 miRBase_mature index: ${mature}"
if( !hairpin.exists() ) exit 1, "Missing Bowtie 1 miRBase_hairpin index: ${hairpin}"
if( !index.exists() ) exit 1, "Missing Bowtie 2 index: ${index}"
if( !gtf.exists() )   exit 2, "Missing GTF annotation: ${gtf}"
if( !bed12.exists() ) exit 2, "Missing BED12 annotation: ${bed12}"
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"


/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: -1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { raw_reads_fastqc; raw_reads_trimgalore }


/*
 * STEP 1 - FastQC
 */

process fastqc {

    tag "$reads"

    memory { 2.GB * task.attempt }
    time { 4.h * task.attempt }
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    file reads from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q ${reads}
    """
}


/*
 * STEP 2 - Trim Galore!
 */

process trim_galore {

    tag "$reads"

    cpus 2
    memory { 4.GB * task.attempt }
    time { 8.h * task.attempt }
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

    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    publishDir "${params.outdir}/bowtie/miRBase_mature", mode: 'copy', pattern: '*.mature_unmapped.fq.gz'

    input:
    file reads from trimmed_reads_bowtie

    output:
    file '*.mature.bam' into miRBase_mature_bam
    file '*.mature_unmapped.fq.gz' into mature_unmapped_reads

    script:
    /* Note! Only read 1 is used for alignment. */

    """
    f='$reads';f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_trimmed};f=\${f%_1};f=\${f%.R1};f=\${f%_R1}

    bowtie \\
      $mature \\
      -q <(zcat $reads) \\
      -p 2 \\
      -t \\
      -n 0 \\
      -l 15 \\
      -e 99999 \\
      -k 10 \\
      --best \\
      --chunkmbs 2048 \\
      --un \${f}.mature_unmapped.fq \\
      -S \\
      | samtools view -bS - > \${f}.mature.bam

    gzip \${f}.mature_unmapped.fq
    """
}

/*
 * STEP 4 - Bowtie against miRBase hairpin
 */

process bowtie_miRBase_hairpin {

    tag "$reads"

    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    publishDir "${params.outdir}/bowtie/miRBase_hairpin", mode: 'copy', pattern: '*.hairpin_unmapped.fq.gz'

    input:
    file reads from mature_unmapped_reads

    output:
    file '*.hairpin.bam' into miRBase_hairpin_bam
    file '*.hairpin_unmapped.fq.gz' into hairpin_unmapped_reads

    script:
    """
    f='$reads';f=\${f%.mature_unmapped.fq.gz}

    bowtie \\
      $hairpin \\
      -p 2 \\
      -t \\
      -n 1 \\
      -l 15 \\
      -e 99999 \\
      -k 10 \\
      --best \\
      --chunkmbs 2048 \\
      -q <(zcat $reads) \\
      --un \${f}.hairpin_unmapped.fq \\
      -S \\
      | samtools view -bS - > \${f}.hairpin.bam

    gzip \${f}.hairpin_unmapped.fq
    """
}


def wrap_mature_and_hairpin(file){

      if ( file.contains("mature") )
        return "miRBase_mature/$file"

      if ( file.contains("hairpin") )
        return "miRBase_hairpin/$file"
}


/*
 * STEP 5 - Post-alignment processing for miRBase mature and hairpin
 */

process miRBasePostAlignment {

    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    publishDir "${params.outdir}/bowtie", mode: 'copy', saveAs: this.&wrap_mature_and_hairpin

    input:
    file input from miRBase_mature_bam .mix(miRBase_hairpin_bam)

    output:
    file '*.count' into miRBase_counts
    file '*.sorted.bam' into miRBase_bam
    file '*.sorted.bam.bai' into miRBase_bai

    script:
    """
    f='$input'; f=\${f%.bam}

    samtools sort \${f}.bam \${f}.sorted
    samtools index \${f}.sorted.bam
    samtools idxstats \${f}.sorted.bam > \${f}.count
    """
}


/*
 * STEP 6 - edgeR miRBase feature counts processing
 */

process edgeR_miRBase {

    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    publishDir "${params.outdir}/edgeR", mode: 'copy', saveAs: this.&wrap_mature_and_hairpin

    input:
    file input_files from miRBase_counts.toSortedList()

    output:
    file '*.{txt,pdf}' into edgeR_miRBase_results

    script:
    """
    #!/usr/bin/env Rscript

    # Load / install required packages
    .libPaths( c( "${params.rlocation}", .libPaths() ) )
    if (!require("limma")){
        source("http://bioconductor.org/biocLite.R")
        biocLite("limma", suppressUpdates=TRUE, lib="${params.rlocation}")
        library("limma")
    }

    if (!require("edgeR")){
        source("http://bioconductor.org/biocLite.R")
        biocLite("edgeR", suppressUpdates=TRUE, lib="${params.rlocation}")
        library("edgeR")
    }

    if (!require("statmod")){
        install.packages("statmod", dependencies=TRUE, repos='http://cloud.r-project.org/', lib="${params.rlocation}")
        library("statmod")
    }

    if (!require("data.table")){
        install.packages("data.table", dependencies=TRUE, repos='http://cloud.r-project.org/', lib="${params.rlocation}")
        library("data.table")
    }

    if (!require("gplots")) {
        install.packages("gplots", dependencies=TRUE, repos='http://cloud.r-project.org/', lib="${params.rlocation}")
        library("gplots")
    }

    if (!require("methods")) {
        install.packages("methods", dependencies=TRUE, repos='http://cloud.r-project.org/', lib="${params.rlocation}")
        library("methods")
    }

    # Read in files
    datafiles = c( "${(input_files as List).join('", "')}" )

    # Put mature and hairpin count files in separated file lists
    filelist<-list()
    filelist[[1]]<-datafiles[grep(".mature.count",datafiles)]
    filelist[[2]]<-datafiles[grep(".hairpin.count",datafiles)]
    names(filelist)<-c("mature","hairpin")

    for (i in 1:2) {

    header<-names(filelist)[i]

    # Prepare the combined data frame with gene ID as rownames and sample ID as colname
    data<-do.call("cbind", lapply(filelist[[i]], fread, header=FALSE, select=c(3)))
    data<-as.data.frame(data)

    temp <- fread(filelist[[i]][1],header=FALSE, select=c(1))
    rownames(data)<-temp\$V1
    colnames(data)<-gsub(".count","",basename(filelist[[i]]))

    data<-data[rownames(data)!="*",]

    # Remove genes with 0 reads in all samples
    row_sub = apply(data, 1, function(row) all(row ==0 ))
    data<-data[!row_sub,]

    # Normalization
    dataDGE<-DGEList(counts=data,genes=rownames(data))
    o <- order(rowSums(dataDGE\$counts), decreasing=TRUE)
    dataDGE <- dataDGE[o,]
    dataNorm <- calcNormFactors(dataDGE)

    # Print normalized read counts to file
    dataNorm_df<-as.data.frame(cpm(dataNorm))
    write.table(dataNorm_df,file=paste(header,"_normalized_CPM.txt",sep=""),sep='\\t',quote=FALSE)

    # Print heatmap based on normalized read counts
    pdf(paste(header,"_CPM_heatmap.pdf",sep=""))
    heatmap.2(cpm(dataNorm),col=redgreen(100),key=TRUE,scale="row",density.info="none",trace="none")
    dev.off()

    # Make MDS plot (only perform with 3 or more samples)
    if (ncol(data)>2){
      pdf(paste(header,"_edgeR_MDS_plot.pdf",sep=""))
      MDSdata <- plotMDS(dataNorm)
      dev.off()
    }

    # Print distance matrix to file
    write.table(MDSdata\$distance.matrix, paste(header,"_edgeR_MDS_distance_matrix.txt",sep=""), quote=FALSE, sep="\\t")

    # Print plot x,y co-ordinates to file
    MDSxy = MDSdata\$cmdscale.out
    colnames(MDSxy) = c(paste(MDSdata\$axislabel, '1'), paste(MDSdata\$axislabel, '2'))

    write.table(MDSxy, paste(header,"_edgeR_MDS_plot_coordinates.txt",sep=""), quote=FALSE, sep="\\t")

    # Get the log counts per million values
    logcpm <- cpm(dataNorm, prior.count=2, log=TRUE)

    # Calculate the euclidean distances between samples
    dists = dist(t(logcpm))

    # Plot a heatmap of correlations
    pdf(paste(header,"_log2CPM_sample_distances_heatmap.pdf",sep=""))
    hmap <- heatmap.2(as.matrix(dists),main="Sample Correlations", key.title="Distance", trace="none",dendrogram="row", margin=c(9, 9))
    dev.off()

    # Plot the heatmap dendrogram
    pdf(paste(header,"_log2CPM_sample_distances_dendrogram.pdf",sep=""))
    plot(hmap\$rowDendrogram, main="Sample Dendrogram")
    dev.off()

    # Write clustered distance values to file
    write.table(hmap\$carpet, paste(header,"_log2CPM_sample_distances.txt",sep=""), quote=FALSE, sep="\\t")
    }

    file.create("corr.done")
    #Printing sessioninfo to standard out
    print("Sample correlation info:")
    sessionInfo()
    """
}


/*
 * STEP 7.1 and 7.2 FOR SPECIFIED GENOME ONLY!
 */

if (params.genome != "ALL") {

/*
 * STEP 7.1 - Bowtie 2 against reference genome
 */

process bowtie2 {

    tag "$reads"

    cpus 8
    memory { 64.GB * task.attempt }
    time { 48.h * task.attempt }
    publishDir "${params.outdir}/bowtie2", mode: 'copy'

    input:
    file reads from trimmed_reads_bowtie2

    output:
    file '*.bowtie2.bam' into bowtie2_bam

    script:
    """
    f='$reads';f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_trimmed};f=\${f%_1};f=\${f%.R1};f=\${f%_R1}

    bowtie2 \\
      -x $index \\
      -U $reads \\
      -k 10 \\
      --very-sensitive \\
      -p 8 \\
      -t \\
      | samtools view -bT $index - > \${f}.bowtie2.bam
    """
}


/*
 * STEP 7.2 - NGI-Visualizations of Bowtie 2 alignment statistics
 */

process ngi_visualizations {

    tag "$bowtie2_bam"

    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    publishDir "${params.outdir}/bowtie2/ngi_visualizations", mode: 'copy'

    input:
    file bowtie2_bam

    output:
    file '*.{png,pdf}' into bowtie2_ngi_visualizations

    script:
    /* Note! Ngi-visualizations (https://github.com/NationalGenomicsInfrastructure/ngi_visualizations) needs to be installed! */
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

    cpus 2
    memory { 16.GB * task.attempt }
    time { 4.h * task.attempt }
    publishDir "${params.outdir}/MultiQC", mode: 'copy'


    input:
    file ('fastqc/*') from fastqc_results.toSortedList()
    file ('trim_galore/*') from trimgalore_results.toSortedList()
    file ('edgeR/*') from edgeR_miRBase_results.toSortedList()

    output:
    file '*multiqc_report.html' into multiqc_html
    file '*multiqc_data' into multiqc_data

    script:
    """
    # Load MultiQC with environment module if not already in PATH
    type multiqc >/dev/null 2>&1 || { module load MultiQC; };
    multiqc -f .
    """
}
