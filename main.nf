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
 - 3.1: Bowtie 1 alignment against miRBase mature miRNA
 - 3.2: Post-alignment processing of miRBase mature miRNA counts
 - 3.3: edgeR analysis on miRBase mature miRNA counts
        - TMM normalization and a table of top expression mature miRNA
        - MDS plot clustering samples
        - Heatmap of sample similarities
 - 4.1: Bowtie 1 alignment against miRBase hairpin for the unaligned reads in step 3
 - 4.2: Post-alignment processing of miRBase hairpin counts
 - 4.3: edgeR analysis on miRBase hairpin counts
        - TMM normalization and a table of top expression hairpin
        - MDS plot clustering samples
        - Heatmap of sample similarities
 - 5.1: Bowtie 2 alignment against host reference genome
 - 5.2: Post-alignment processing of Bowtie 2
 - 6:   NGI-Visualization of Bowtie 2 alignment statistics
 - 7:   MultiQC
----------------------------------------------------------------------------------------
*/


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 1.0

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
    /* I'm using the latest version of trimgalore downloaded from github. The version on uppmax has a bug when handling smallRNA-seq data. */

    single = reads instanceof Path

    if(single) {
      """
      perl /home/chuanw/trimgalore/trim_galore.pl --small_rna --gzip $reads
      """
    } else {
      """
      perl /home/chuanw/trimgalore/trim_galore.pl --paired --small_rna --gzip $reads
      """
    }
}


/*
 * STEP 3.1 - Bowtie miRBase mature miRNA
 */

process bowtie_miRBase_mature {

    tag "$name"

    module 'bioinfo-tools'
    module 'bowtie'
    module 'samtools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie_miRBase_mature", mode: 'copy'

    input:
    file(reads:'*') from trimmed_reads_bowtie
    set val(name) from name_for_bowtie_miRBase_mature

    output:
    file '*.mature.sam' into miRBase_mature_sam
    file '*.mature_unmapped.fq.gz' into mature_unmapped_reads

    script:
    /* Note! Only read 1 is used for alignment. */

    single = reads instanceof Path

    if(single) {
      """
      f='$reads';f=(\$f);input=\${f[0]};f=\${input%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.R1_val_1*};f=\${f%.R2_val_2*};f=\${f%_val_1};f=\${f%_val_2};f=\${f%_trimmed};f=\${f%_1}
      prefix=\$f

      bowtie $mature -q <(zcat \$input) -v 0 -k 10 -S \${prefix}.mature.sam --un \${prefix}.mature_unmapped.fq

      gzip \${prefix}.mature_unmapped.fq
      """
    } else {
      """
      f='$reads';f=(\$f);input=\${f[0]};f=\${input%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.R1_val_1*};f=\${f%.R2_val_2*};f=\${f%_val_1};f=\${f%_val_2};f=\${f%_trimmed};f=\${f%_1}
      prefix=\$f

      bowtie $mature -q <(zcat \$input) -v 0 -k 10 -S \${prefix}.mature.sam --un \${prefix}.mature_unmapped.fq

      gzip \${prefix}.mature_unmapped.fq
      """
    }
}

/*
 * STEP 3.2 - Bowtie against miRBase mature miRNA post-alignment processing
 */

process bowtie_miRBase_mature_postalignment {

    tag "$miRBase_mature_sam"

    module 'bioinfo-tools'
    module 'samtools'
    module 'BEDTools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie_miRBase_mature", mode: 'copy'

    input:
    file miRBase_mature_sam

    output:
    file '*.sorted.bam' into miRBase_mature_bam
    file '*.sorted.bam.bai' into miRBase_mature_bai
    file '*.count' into miRBase_mature_count

    script:
    """
    f='$miRBase_mature_sam';f=(\$f);f=\${f[0]};f=\${f%.sam}
    prefix=\$f

    samtools view -bS \${prefix}.sam > \${prefix}.bam
    samtools sort \${prefix}.bam \${prefix}.sorted
    samtools index \${prefix}.sorted.bam
    samtools idxstats \${prefix}.sorted.bam > \${prefix}.count

    rm \${prefix}.sam
    rm \${prefix}.bam
    """
}


/*
 * STEP 3.3 - edgeR miRBase mature miRNA counts processing
 */

process edgeR_miRBase_mature {

    module 'bioinfo-tools'
    module 'R/3.2.3'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie_miRBase_mature/edgeR", mode: 'copy'

    input:
    file input_files from miRBase_mature_count.toSortedList()

    output:
    file '*.{txt,pdf}' into edgeR_miRBase_mature_results

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

    # Prepare the combined data frame with gene ID as rownames and sample ID as colname
    data<-do.call("cbind", lapply(datafiles, fread, header=FALSE, select=c(3)))
    data<-as.data.frame(data)

    temp <- fread(datafiles[1],header=FALSE, select=c(1))
    rownames(data)<-temp\$V1
    colnames(data)<-gsub(".count","",basename(datafiles))

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
    pdf('edgeR_MDS_plot.pdf')
    MDSdata <- plotMDS(dataNorm)
    dev.off()

    # Print distance matrix to file
    write.table(MDSdata\$distance.matrix, 'edgeR_MDS_distance_matrix.txt', quote=FALSE, sep="\\t")

    # Print plot x,y co-ordinates to file
    MDSxy = MDSdata\$cmdscale.out
    colnames(MDSxy) = c(paste(MDSdata\$axislabel, '1'), paste(MDSdata\$axislabel, '2'))
    write.table(MDSxy, 'edgeR_MDS_plot_coordinates.txt', quote=FALSE, sep="\\t")

    # Get the log counts per million values
    logcpm <- cpm(dataNorm, prior.count=2, log=TRUE)

    # Calculate the euclidean distances between samples
    dists = dist(t(logcpm))

    # Plot a heatmap of correlations
    pdf('log2CPM_sample_distances_heatmap.pdf')
    hmap <- heatmap.2(as.matrix(dists),main="Sample Correlations", key.title="Distance", trace="none",dendrogram="row", margin=c(9, 9))
    dev.off()

    # Plot the heatmap dendrogram
    pdf('log2CPM_sample_distances_dendrogram.pdf')
    plot(hmap\$rowDendrogram, main="Sample Dendrogram")
    dev.off()

    # Write clustered distance values to file
    write.table(hmap\$carpet, 'log2CPM_sample_distances.txt', quote=FALSE, sep="\\t")

    file.create("corr.done")
    #Printing sessioninfo to standard out
    print("Sample correlation info:")
    sessionInfo()
    """
}

/*
 * STEP 4.1 - Bowtie against miRBase hairpin
 */

process bowtie_miRBase_hairpin {

    tag "$name"

    module 'bioinfo-tools'
    module 'bowtie'
    module 'samtools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie_miRBase_hairpin", mode: 'copy'

    input:
    file(reads:'*') from mature_unmapped_reads
    set val(name) from name_for_bowtie_miRBase_hairpin

    output:
    file '*.hairpin.sam' into miRBase_hairpin_sam
    file '*.hairpin_unmapped.fq.gz' into hairpin_unmapped_reads

    script:
    /* Note! Only read 1 is used for alignment. */

    single = reads instanceof Path

    if(single) {
      """
      f='$reads';f=(\$f);input=\${f[0]};f=\${input%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.mature_unmapped};f=\${f%.R1_val_1*};f=\${f%.R2_val_2*};f=\${f%_val_1};f=\${f%_val_2};f=\${f%_trimmed};f=\${f%_1}
      prefix=\$f

      bowtie $hairpin -q <(zcat \$input) -v 0 -k 10 -S \${prefix}.hairpin.sam --un \${prefix}.hairpin_unmapped.fq

      gzip \${prefix}.hairpin_unmapped.fq
      """
    } else {
      """
      f='$reads';f=(\$f);input=\${f[0]};f=\${input%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%.mature_unmapped};f=\${f%.R1_val_1*};f=\${f%.R2_val_2*};f=\${f%_val_1};f=\${f%_val_2};f=\${f%_trimmed};f=\${f%_1}
      prefix=\$f

      bowtie $hairpin -q <(zcat \$input) -v 0 -k 10 -S \${prefix}.hairpin.sam --un \${prefix}.hairpin_unmapped.fq

      gzip \${prefix}.hairpin_unmapped.fq
      """
    }
}

/*
 * STEP 4.2 - Bowtie against miRBase hairpin post-alignment processing
 */

process bowtie_miRBase_mature_postalignment {

    tag "$miRBase_hairpin_sam"

    module 'bioinfo-tools'
    module 'samtools'
    module 'BEDTools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie_miRBase_hairpin", mode: 'copy'

    input:
    file miRBase_hairpin_sam

    output:
    file '*.sorted.bam' into miRBase_hairpin_bam
    file '*.sorted.bam.bai' into miRBase_hairpin_bai
    file '*.count' into miRBase_hairpin_count

    script:
    """
    f='$miRBase_hairpin_sam';f=(\$f);f=\${f[0]};f=\${f%.sam}
    prefix=\$f

    samtools view -bS \${prefix}.sam > \${prefix}.bam
    samtools sort \${prefix}.bam \${prefix}.sorted
    samtools index \${prefix}.sorted.bam
    samtools idxstats \${prefix}.sorted.bam > \${prefix}.count

    rm \${prefix}.sam
    rm \${prefix}.bam
    """
}


/*
 * STEP 4.3 - edgeR miRBase hairpin counts processing
 */

process edgeR_miRBase_hairpin {

    module 'bioinfo-tools'
    module 'R/3.2.3'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie_miRBase_hairpin/edgeR", mode: 'copy'

    input:
    file input_files from miRBase_hairpin_count.toSortedList()

    output:
    file '*.{txt,pdf}' into edgeR_miRBase_hairpin_results

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

    # Prepare the combined data frame with gene ID as rownames and sample ID as colname
    data<-do.call("cbind", lapply(datafiles, fread, header=FALSE, select=c(3)))
    data<-as.data.frame(data)

    temp <- fread(datafiles[1],header=FALSE, select=c(1))
    rownames(data)<-temp\$V1
    colnames(data)<-gsub(".count","",basename(datafiles))

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
    pdf('edgeR_MDS_plot.pdf')
    MDSdata <- plotMDS(dataNorm)
    dev.off()

    # Print distance matrix to file
    write.table(MDSdata\$distance.matrix, 'edgeR_MDS_distance_matrix.txt', quote=FALSE, sep="\\t")

    # Print plot x,y co-ordinates to file
    MDSxy = MDSdata\$cmdscale.out
    colnames(MDSxy) = c(paste(MDSdata\$axislabel, '1'), paste(MDSdata\$axislabel, '2'))
    write.table(MDSxy, 'edgeR_MDS_plot_coordinates.txt', quote=FALSE, sep="\\t")

    # Get the log counts per million values
    logcpm <- cpm(dataNorm, prior.count=2, log=TRUE)

    # Calculate the euclidean distances between samples
    dists = dist(t(logcpm))

    # Plot a heatmap of correlations
    pdf('log2CPM_sample_distances_heatmap.pdf')
    hmap <- heatmap.2(as.matrix(dists),main="Sample Correlations", key.title="Distance", trace="none",dendrogram="row", margin=c(9, 9))
    dev.off()

    # Plot the heatmap dendrogram
    pdf('log2CPM_sample_distances_dendrogram.pdf')
    plot(hmap\$rowDendrogram, main="Sample Dendrogram")
    dev.off()

    # Write clustered distance values to file
    write.table(hmap\$carpet, 'log2CPM_sample_distances.txt', quote=FALSE, sep="\\t")

    file.create("corr.done")
    #Printing sessioninfo to standard out
    print("Sample correlation info:")
    sessionInfo()
    """
}


/*
 * STEP 5.1 - Bowtie 2 against reference genome
 */
process bowtie2 {

    tag "$name"

    module 'bioinfo-tools'
    module 'bowtie2'
    module 'samtools'

    cpus 4
    memory { 32.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie2", mode: 'copy'

    input:
    file(reads:'*') from trimmed_reads_bowtie2
    set val(name) from name_for_bowtie2

    output:
    file '*.sam' into bowtie2_sam

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
        -S \${prefix}.sam
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
        -S \${prefix}.sam
      """
    }
}


/*
 * STEP 5.2 - Bowtie 2 post-alignment processing
 */

process bowtie2_postalignment {

    tag "$name"

    module 'bioinfo-tools'
    module 'samtools'
    module 'BEDTools'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie2", mode: 'copy'

    input:
    file bowtie2_sam
    set val(name) from name_for_samtools

    output:
    file '*.sorted.bam' into bowtie2_bam
    file '*.sorted.bam.bai' into bowtie2_bai
    file '*.sorted.bed' into bowtie2_bed
    stdout into bowtie2_logs

    script:
    """
    f='$bowtie2_sam';f=(\$f);f=\${f[0]};f=\${f%.sam}
    prefix=\$f

    samtools view -bT $index \${prefix}.sam > \${prefix}.bam
    samtools sort \${prefix}.bam \${prefix}.sorted
    samtools index \${prefix}.sorted.bam
    rm \${prefix}.sam
    rm \${prefix}.bam
    bedtools bamtobed -i \${prefix}.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > \${prefix}.sorted.bed
    """
}


/*
 * STEP 6 - NGI-Visualizations of Bowtie 2 alignment statistics
 */

process ngi_visualizations {

    tag "$bowtie2_bam"

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/bowtie2/ngi_visualizations", mode: 'copy'
    errorStrategy 'ignore'

    input:
    file bowtie2_bam

    output:
    file '*.{png,pdf}' into bowtie2_ngi_visualizations

    script:
    """
    python /home/chuanw/ngi_visualizations/count_biotypes.py -g $gtf $bowtie2_bam
    """
}


/*
 * STEP 7 - MultiQC
 */

process multiqc {

    module 'bioinfo-tools'
    module 'MultiQC'

    cpus 2
    memory { 16.GB * task.attempt }
    time { 120.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    errorStrategy 'ignore'

    input:
    file ('fastqc/*') from fastqc_results.toSortedList()
    file ('trim_galore/*') from trimgalore_results.toSortedList()
    file ('bowtie_miRBase_mature/edgeR/*') from edgeR_miRBase_mature_results.toSortedList()
    file ('bowtie_miRBase_hairpin/edgeR/*') from edgeR_miRBase_hairpin_results.toSortedList()

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
