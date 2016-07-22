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
params.gtf   = params.genomes[ params.genome ].gtf
params.bed12 = params.genomes[ params.genome ].bed12
params.bowtie2 = params.genomes[ params.genome ].bowtie2
params.bowtie_miRBase_mature = params.genomes.bowtie_miRBase_mature
params.bowtie_miRBase_hairpin = params.genomes.bowtie_miRBase_hairpin
params.bowtie_rfam = params.genomes.bowtie_rfam
params.name = "miRNA-Seq Best practice"

// Input files
params.input = "data/*.fastq.gz"

// Output path
params.out = "$PWD"

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs=file(params.rlocation)

log.info "===================================="
log.info " RNAbp : Small RNA-Seq Best Practice v${version}"
log.info "===================================="
log.info "Reads        : ${params.input}"
log.info "Genome       : ${params.genome}"
log.info "Index        : ${params.index}"
log.info "Annotation   : ${params.gtf}"
log.info "miRBase mature  : ${params.bowtie_miRBase_mature}"
log.info "miRBase hairpin : ${params.bowtie_miRBase_hairpin}"
log.info "Rfam         : ${params.bowtie_rfam}"
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
if( !bowtie_miRBase_mature.exists() ) exit 1, "Missing Bowtie 1 miRBase_mature index: ${bowtie_miRBase_mature}"
if( !bowtie_miRBase_hairpin.exists() ) exit 1, "Missing Bowtie 1 miRBase_hairpin index: ${bowtie_miRBase_hairpin}"
if( !bowtie_rfam.exists() ) exit 1, "Missing Bowtie 1 rfam index: ${bowtie_rfam}"
if( !bowtie2.exists() ) exit 1, "Missing Bowtie 2 index: ${bowtie2}"
if( !gtf.exists() )   exit 2, "Missing GTF annotation: ${gtf}"
if( !bed12.exists() ) exit 2, "Missing BED12 annotation: ${bed12}"

// Setting up a directory to save results to 
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
    /* Trim Galore should automatically detect the small RNA adapter sequence */
    """
    trim_galore --gzip --fastqc_args "-q" ${reads}
    """
}


/*
 * STEP 3 - Bowtie against miRBase mature RNA
 */
process bowtie_miRBase_mature {
    
    module 'bioinfo-tools'
    module 'bowtie'
    module 'samtools'
    
    memory '2 GB'
    time '6h'

    input:
    file(reads:'*') from trimmed_reads 
    
    output:
    file '*.bam' into miRBase_mature_bam 
   
    """
    f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1}
    bowtie \
        -p 2 \
        -t \
        -n 0 \
        -l 15 \
        -e 99999 \
        -k 200 \
        --best \
        -S \
        --chunkmbs 2048 \
        $bowtie_miRBase_mature \
        $reads \
        | samtools view -bS > \${f}.aligned_miRBase_mature.bam 
    """
}


/*
 * STEP 3 - Bowtie against miRBase hairpin RNA
 */
process bowtie_miRBase_mature {
    
    module 'bioinfo-tools'
    module 'bowtie'
    module 'samtools'
    
    memory '2 GB'
    time '6h'

    input:
    file(reads:'*') from trimmed_reads 
    
    output:
    file '*.bam' into miRBase_hairpin_bam
   
    """
    f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1}
    bowtie \
        -p 2 \
        -t \
        -n 0 \
        -l 15 \
        -e 99999 \
        -k 200 \
        --best \
        -S \
        --chunkmbs 2048 \
        $bowtie_miRBase_hairpin \
        $reads \
        | samtools view -bS > \${f}.aligned_miRBase_hairpin.bam 
    """
}



/*
 * STEP 3 - Bowtie against rfam
 */
process bowtie_rfam {
    
    module 'bioinfo-tools'
    module 'bowtie'
    module 'samtools'
    
    memory '2 GB'
    time '6h'

    input:
    file(reads:'*') from trimmed_reads 
    
    output:
    file '*.bam' into rfam_bam
   
    """
    f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1}
    bowtie \
        -p 2 \
        -t \
        -n 0 \
        -l 15 \
        -e 99999 \
        -k 200 \
        --best \
        -S \
        --chunkmbs 2048 \
        $bowtie_rfam \
        $reads \
        | samtools view -bS > \${f}.aligned_rfam.bam 
    """
}



/*
 * STEP 3 - Samtools idxstats to count bowtie 1 alignments
 */
process samtools_idxstats {
    
    module 'bioinfo-tools'
    module 'samtools'
    
    memory '8 GB'
    time '6h'

    input:
    file '*.bam' from miRBase_mature_bam
    // TODO: WORK OUT HOW TO GET ALL THREE CHANNELS TO COME HERE
    // ALSO: miRBase_hairpin_bam, rfam_bam 
    
    output:
    file '*.counts.txt' into samtools_counts
    file '*.sorted.bam'
    file '*.sorted.bam.bai'
   
    """
    samtools sort $miRBase_mature_bam > ${miRBase_mature_bam}.sorted.bam
    samtools index ${miRBase_mature_bam}.sorted.bam
    samtools idxstats ${miRBase_mature_bam}.sorted.bam > ${miRBase_mature_bam}_counts.txt
    """
    
}



/*
 * STEP 3 - Bowtie 2 against host genome
 */
process bowtie2 {
    
    module 'bioinfo-tools'
    module 'bowtie2'
    module 'samtools'
    
    memory '8 GB'
    time '6h'

    input:
    file(reads:'*') from trimmed_reads 
    
    output:
    file '*.bam' into bowtie2_bam
   
    """
    f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1}
    bowtie2 \
        -p 1 \
        -t \
        -x $bowtie2 \
        -U $reads \
        | samtools view -bS > \${f}.aligned_${params.genome}.bam 
    """
}



/*
 * Filter aligned reads based on alignment length
 */
process bowtie2 {
    
    module 'bioinfo-tools'
    module 'NGSUtils'
    
    memory '8 GB'
    time '6h'

    input:
    file '*.bam' from bowtie2_bam 
    
    output:
    file '*_filtered.bam' into bowtie2_bam_16_30nt
   
    """
    bamutils filter \
        $bowtie2_bam \
        ${bowtie2_bam}_filtered.bam \
        -minlen 16 \
        -maxlen 30 \
        -mapped
    """
}

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
    file '*.bam' from bowtie2_bam
    file gtf from gtf
    
    output:
    file '*_gene.featureCounts.txt' 
    file '*_biotype.featureCounts.txt' 
    file '*_rRNA_counts.txt' 
    file '*.summary' 
    file 'featureCounts.done' into featureCounts_done    
    """
    featureCounts -a $gtf -g gene_id -o ${bam}_gene.featureCounts.txt -p -s 2 $bowtie2_bam
    featureCounts -a $gtf -g gene_biotype -o ${bam}_biotype.featureCounts.txt -p -s 2 $bowtie2_bam
    cut -f 1,7 ${bam}_biotype.featureCounts.txt | sed '1,2d' | grep 'rRNA' > ${bam}_rRNA_counts.txt
    echo done >featureCounts.done
    """
}





/*
 * STEP 10 - edgeR MDS and heatmap
 */
/*
///////////////////////////////////////
// UNDER DEVELOPMENT
// Should still work the same as the mRNA pipeline though..?
///////////////////////////////////////


process sample_correlation {
    module 'bioinfo-tools'
    module 'R/3.2.3'
    
    memory { 16.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/sample_correlation", mode: 'copy'
    
    input:
    file input_files from geneCounts.toList()
    bam_count
    
    output:
    file '*.{txt,pdf}' into sample_correlation_results
    
    when:
    num_bams > 2 && (!params.sampleLevel)

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
    
    if (!require("data.table")){
        install.packages("data.table", dependencies=TRUE, repos='http://cloud.r-project.org/', lib="${params.rlocation}")
        library("data.table")
    }
    
    if (!require("gplots")) {
        install.packages("gplots", dependencies=TRUE, repos='http://cloud.r-project.org/', lib="${params.rlocation}")
        library("gplots")
    }
    
    # Load input counts data
    datafiles = c( "${(input_files as List).join('", "')}" )
    
    # Load count column from all files into a list of data frames
    # Use data.tables fread as much much faster than read.table
    # Row names are GeneIDs
    temp <- lapply(datafiles, fread, skip="Geneid", header=TRUE, colClasses=c(NA, rep("NULL", 5), NA))
    
    # Merge into a single data frame
    merge.all <- function(x, y) {
        merge(x, y, all=TRUE, by="Geneid")
    }
    data <- data.frame(Reduce(merge.all, temp))
    
    # Clean sample name headers
    colnames(data) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(data))
    
    # Set GeneID as row name
    rownames(data) <- data[,1]
    data[,1] <- NULL
    
    # Convert data frame to edgeR DGE object
    dataDGE <- DGEList( counts=data.matrix(data) )
    
    # Normalise counts
    dataNorm <- calcNormFactors(dataDGE)
    
    # Make MDS plot
    pdf('edgeR_MDS_plot.pdf')
    MDSdata <- plotMDS(dataNorm)
    dev.off()
    
    # Print distance matrix to file
    write.table(MDSdata\$distance.matrix, 'edgeR_MDS_distance_matrix.txt', quote=FALSE, sep="\t")
    
    # Print plot x,y co-ordinates to file
    MDSxy = MDSdata\$cmdscale.out
    colnames(MDSxy) = c(paste(MDSdata\$axislabel, '1'), paste(MDSdata\$axislabel, '2'))
    write.table(MDSxy, 'edgeR_MDS_plot_coordinates.txt', quote=FALSE, sep="\t")
    
    # Get the log counts per million values
    logcpm <- cpm(dataNorm, prior.count=2, log=TRUE)
    
    # Calculate the euclidean distances between samples
    dists = dist(t(logcpm))
    
    # Plot a heatmap of correlations
    pdf('log2CPM_sample_distances_heatmap.pdf')
    hmap <- heatmap.2(as.matrix(dists),
      main="Sample Correlations", key.title="Distance", trace="none",
      dendrogram="row", margin=c(9, 9)
    )
    dev.off()
    
    # Plot the heatmap dendrogram
    pdf('log2CPM_sample_distances_dendrogram.pdf')
    plot(hmap\$rowDendrogram, main="Sample Dendrogram")
    dev.off()
    
    # Write clustered distance values to file
    write.table(hmap\$carpet, 'log2CPM_sample_distances.txt', quote=FALSE, sep="\t")
    
    file.create("corr.done")
    """
}


*/





/*
 * STEP 4
 */
/*
//////////////////////////////////
// Under development
//////////////////////////////////
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




