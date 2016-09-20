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
 $ nextflow rnaseq.nf --reads path/to/data/*fq.gz
----------------------------------------------------------------------------------------
 Pipeline overview:
 - 1:   FastQC for raw sequencing reads quility control
 - 2:   Trim Galore! for adapter trimming
 - 3:   Bowtie 1 alignment against miRBase mature miRNA
 - 4:   Bowtie 1 alignment against miRBase hairpin for the unaligned reads in step 3
 - 5:   Post-alignment processing for miRBase mature and hairpin
 - 6:   edgeR analysis on miRBase mature miRNA counts
          - TMM normalization and a table of top expression mature miRNA
          - MDS plot clustering samples
          - Heatmap of sample similarities
 - 7.1: Bowtie 2 alignment against host reference genome
 - 7.2: Post-alignment processing for Bowtie 2
 - 7.3: NGI-Visualization of Bowtie 2 alignment statistics
 - 8:   MultiQC
----------------------------------------------------------------------------------------
*/


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 1.1

// Reference genome index
params.genome = 'GRCh37'
params.gtf   = params.genomes[ params.genome ].gtf
params.bed12 = params.genomes[ params.genome ].bed12
params.index = params.genomes[ params.genome ].bowtie2
params.mature  = params.miRBase_mature[ params.genome ].bowtie
params.hairpin = params.miRBase_hairpin[ params.genome ].bowtie
params.name = "miRNA-Seq Best practice"
params.outdir = './results'

// Input files
params.reads = "data/*{1,2}*.fastq.gz"

// Output path
params.out = "$PWD"

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs=file(params.rlocation)

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
log.info "===================================="

// Create R library directories if not already existing
nxtflow_libs.mkdirs()

// Set up nextflow objects
gtf = file(params.gtf)
bed12 = file(params.bed12)
index = file(params.index)
mature = file(params.mature)
hairpin = file(params.hairpin)

// Validate inputs(NOT WORKING!!)
// if( !mature.exists() ) exit 1, "Missing Bowtie 1 miRBase_mature index: ${mature}"
// if( !hairpin.exists() ) exit 1, "Missing Bowtie 1 miRBase_hairpin index: ${hairpin}"
// if( !index.exists() ) exit 1, "Missing Bowtie 2 index: ${index}"
// if( !gtf.exists() )   exit 2, "Missing GTF annotation: ${gtf}"
// if( !bed12.exists() ) exit 2, "Missing BED12 annotation: ${bed12}"

// Setting up a directory to save results to
results_path = './results'

/*
 * Create a channel for input read files
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

read_files.into { fastQC_input; trimgalore_input; name_for_bowtie_miRBase_mature; name_for_bowtie_miRBase_hairpin; name_for_bowtie2; name_for_samtools }


/*
 * STEP 1 - FastQC
 */

process fastqc {

    tag "$name"

    module 'bioinfo-tools'
    module 'FastQC'

    memory { 2.GB * task.attempt }
    time { 4.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'warning' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    set val(name), file(reads:'*') from fastQC_input

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

    tag "$name"

    module 'bioinfo-tools'
    module 'FastQC'
    module 'cutadapt'
    module 'TrimGalore'

    cpus 2
    memory { 3.GB * task.attempt }
    time { 8.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    set val(name), file(reads:'*') from trimgalore_input

    output:
    file '*.gz' into trimmed_reads_bowtie, trimmed_reads_bowtie2
    file '*trimming_report.txt' into trimgalore_results

    script:
    /* Note! Use Trim_Galore Version 0.4.2 (Released 2016-09-07) or newer versions! */

    single = reads instanceof Path

    if(single) {
      """
      trim_galore --small_rna --gzip $reads
      """
    } else {
      """
      trim_galore --paired --small_rna --gzip $reads
      """
    }
}


/*
 * STEP 3 - Bowtie miRBase mature miRNA
 */

process bowtie_miRBase_mature {

    tag "$name"

    module 'bioinfo-tools'
    module 'bowtie'
    module 'samtools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie_miRBase_mature", mode: 'copy'

    input:
    file(reads:'*') from trimmed_reads_bowtie
    set val(name) from name_for_bowtie_miRBase_mature

    output:
    file '*.mature.bam' into miRBase_mature_bam
    file '*.mature_unmapped.fq.gz' into mature_unmapped_reads

    script:
    /* Note! Only read 1 is used for alignment. */

    single = reads instanceof Path

    if(single) {
      """
      f='$reads';f=(\$f);input=\${f[0]};f=\${input%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.R1_val_1*};f=\${f%.R2_val_2*};f=\${f%_val_1};f=\${f%_val_2};f=\${f%_trimmed};f=\${f%_1}
      prefix=\$f

      bowtie \\
        $mature \\
        -q <(zcat \$input) \\
        -p 2 \\
        -t \\
        -n 0 \\
        -l 15 \\
        -e 99999 \\
        -k 10 \\
        --best \\
        --chunkmbs 2048 \\
        --un \${prefix}.mature_unmapped.fq \\
        -S \\
        | samtools view -bS - > \${prefix}.mature.bam

      gzip \${prefix}.mature_unmapped.fq
      """
    } else {
      """
      f='$reads';f=(\$f);input=\${f[0]};f=\${input%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.R1_val_1*};f=\${f%.R2_val_2*};f=\${f%_val_1};f=\${f%_val_2};f=\${f%_trimmed};f=\${f%_1}
      prefix=\$f

      bowtie \\
        $mature \\
        -q <(zcat \$input) \\
        -p 2 \\
        -t \\
        -n 0 \\
        -l 15 \\
        -e 99999 \\
        -k 10 \\
        --best \\
        --chunkmbs 2048 \\
        --un \${prefix}.mature_unmapped.fq \\
        -S \\
        | samtools view -bS - > \${prefix}.mature.bam

      gzip \${prefix}.mature_unmapped.fq
      """
    }
}

/*
 * STEP 4 - Bowtie against miRBase hairpin
 */

process bowtie_miRBase_hairpin {

    tag "$name"

    module 'bioinfo-tools'
    module 'bowtie'
    module 'samtools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie_miRBase_hairpin", mode: 'copy'

    input:
    file(reads:'*') from mature_unmapped_reads
    set val(name) from name_for_bowtie_miRBase_hairpin

    output:
    file '*.hairpin.bam' into miRBase_hairpin_bam
    file '*.hairpin_unmapped.fq.gz' into hairpin_unmapped_reads

    script:
    /* Note! Only read 1 is used for alignment. */

    single = reads instanceof Path

    if(single) {
      """
      f='$reads';f=(\$f);input=\${f[0]};f=\${input%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.mature_unmapped};f=\${f%.R1_val_1*};f=\${f%.R2_val_2*};f=\${f%_val_1};f=\${f%_val_2};f=\${f%_trimmed};f=\${f%_1}
      prefix=\$f

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
        -q <(zcat \$input) \\
        --un \${prefix}.hairpin_unmapped.fq \\
        -S \\
        | samtools view -bS - > \${prefix}.hairpin.bam

      gzip \${prefix}.hairpin_unmapped.fq
      """
    } else {
      """
      f='$reads';f=(\$f);input=\${f[0]};f=\${input%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.mature_unmapped};f=\${f%.R1_val_1*};f=\${f%.R2_val_2*};f=\${f%_val_1};f=\${f%_val_2};f=\${f%_trimmed};f=\${f%_1}
      prefix=\$f

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
        -q <(zcat \$input) \\
        --un \${prefix}.hairpin_unmapped.fq \\
        -S \\
        | samtools view -bS - > \${prefix}.hairpin.bam

      gzip \${prefix}.hairpin_unmapped.fq
      """
    }
}


/*
 * STEP 5 - Post-alignment processing for miRBase mature and hairpin
 */

process miRBasePostAlignment {

    module 'bioinfo-tools'
    module 'samtools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/miRBaseCounts", mode: 'copy'

    input:
    file(input:'*') from miRBase_mature_bam .mix(miRBase_hairpin_bam)

    output:
    file '*.count' into miRBase_counts

    script:
    """
    f='$input';f=(\$f);f=\${f[0]};f=\${f%.bam}
    prefix=\$f

    samtools sort \${prefix}.bam \${prefix}.sorted
    samtools index \${prefix}.sorted.bam
    samtools idxstats \${prefix}.sorted.bam > \${prefix}.count
    """
}


/*
 * STEP 6 - edgeR miRBase mature miRNA counts processing
 */

process edgeR_miRBase_mature {

    module 'bioinfo-tools'
    module 'R/3.2.3'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/edgeR", mode: 'copy'

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

    # Make MDS plot
    pdf(paste(header,"_edgeR_MDS_plot.pdf",sep=""))
    MDSdata <- plotMDS(dataNorm)
    dev.off()

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
 * STEP 7.1 - Bowtie 2 against reference genome
 */
process bowtie2 {

    tag "$name"

    module 'bioinfo-tools'
    module 'bowtie2'
    module 'samtools'

    cpus 4
    memory { 32.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie2", mode: 'copy'

    input:
    file(reads:'*') from trimmed_reads_bowtie2
    set val(name) from name_for_bowtie2

    output:
    file '*.bowtie2.bam' into bowtie2_bam

    script:
    /* Note! Only read 1 is used for alignment. */

    single = reads instanceof Path

    if(single) {
      """
      f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.R1_val_1*};f=\${f%.R2_val_2*};f=\${f%_val_1};f=\${f%_val_2};f=\${f%_trimmed};f=\${f%_1}
      prefix=\$f

      bowtie2 \\
        -x $index \\
        -U $reads \\
        -k 10 \\
        --very-sensitive \\
        -p 1 \\
        -t \\
        | samtools view -bT $index - > \${prefix}.bowtie2.bam
      """
    } else {
      """
      f='$reads';f=(\$f);input=\${f[0]};f=\${input%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.R1_val_1*};f=\${f%.R2_val_2*};f=\${f%_val_1};f=\${f%_val_2};f=\${f%_trimmed};f=\${f%_1}
      prefix=\$f

      bowtie2 \\
        -x $index \\
        -U \$input \\
        -k 10 \\
        --very-sensitive \\
        -p 1 \\
        -t \\
        | samtools view -bT $index - > \${prefix}.bowtie2.bam
      """
    }
}


/*
 * STEP 7.2 - Bowtie 2 post-alignment processing
 */

process bowtie2_postalignment {

    tag "$name"

    module 'bioinfo-tools'
    module 'samtools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie2", mode: 'copy'

    input:
    file bowtie2_bam
    set val(name) from name_for_samtools

    output:
    file '*.sorted.bam' into bowtie2_bam_sorted
    file '*.sorted.bam.bai' into bowtie2_bai
    stdout into bowtie2_logs

    script:
    """
    f='$bowtie2_bam';f=(\$f);f=\${f[0]};f=\${f%.bam}
    prefix=\$f

    samtools sort \${prefix}.bam \${prefix}.sorted
    samtools index \${prefix}.sorted.bam
    rm \${prefix}.bam
    """
}



/*
 * STEP 7.3 - NGI-Visualizations of Bowtie 2 alignment statistics
 */

process ngi_visualizations {

    tag "$bowtie2_bam_sorted"

    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie2/ngi_visualizations", mode: 'copy'
    errorStrategy 'ignore'

    input:
    file bowtie2_bam_sorted

    output:
    file '*.{png,pdf}' into bowtie2_ngi_visualizations

    script:
    /* Note! Ngi-visualizations (https://github.com/NationalGenomicsInfrastructure/ngi_visualizations) needs to be installed! */
    """
    #!/usr/bin/env python
    from ngi_visualizations.biotypes import count_biotypes
    count_biotypes.main('$gtf','$bowtie2_bam_sorted')
    """
}


/*
 * STEP 8 - MultiQC
 */

process multiqc {

    module 'bioinfo-tools'
    module 'MultiQC'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 4.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    errorStrategy 'ignore'

    input:
    file ('fastqc/*') from fastqc_results.toSortedList()
    file ('trim_galore/*') from trimgalore_results.toSortedList()
    file ('edgeR/*') from edgeR_miRBase_results.toSortedList()

    output:
    file '*multiqc_report.html'
    file '*multiqc_data'

    script:
    """
    # Load MultiQC with environment module if not already in PATH
    type multiqc >/dev/null 2>&1 || { module load multiqc; };
    multiqc -f -t ngi .
    """
}


/*
 * Helper function, given a file Path
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 *
 * For example:
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 *
 * Returns:
 *   'file_alpha'
 */
def readPrefix( Path actual, template ) {

    final fileName = actual.getFileName().toString()

    def filePattern = template.toString()
    int p = filePattern.lastIndexOf('/')
    if( p != -1 ) filePattern = filePattern.substring(p+1)
    if( !filePattern.contains('*') && !filePattern.contains('?') )
        filePattern = '*' + filePattern

    def regex = filePattern
        .replace('.','\\.')
        .replace('*','(.*)')
        .replace('?','(.?)')
        .replace('{','(?:')
        .replace('}',')')
        .replace(',','|')

    def matcher = (fileName =~ /$regex/)
    if( matcher.matches() ) {
        def end = matcher.end(matcher.groupCount() )
        def prefix = fileName.substring(0,end)
        while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') )
            prefix=prefix[0..-2]
        return prefix
    }
    return fileName
}
