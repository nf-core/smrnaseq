#!/usr/bin/env Rscript

# Originally written by Phil Ewels and Chuan Wang and released under the MIT license.
# Contributions by Alexander Peltzer, Anabella Trigila, James Fellows Yates, Sarah Djebali, Kevin Menden, Konrad Stawinski and Lorena Pantano also released under the MIT license. See LICENSE https://github.com/nf-core/smrnaseq/blob/master/LICENSE for details.

# Command line arguments
args = commandArgs(trailingOnly=TRUE)

input <- as.character(args[1:length(args)])
library("limma")
library("edgeR")
library("statmod")
library("data.table")
library("gplots")
library("methods")

# Put mature and hairpin count files in separated file lists
filelist<-list()
filelist[[1]]<-input[grep(".mature.sorted",input)]
filelist[[2]]<-input[grep(".hairpin.sorted",input)]
names(filelist)<-c("mature","hairpin")
print(filelist)

for (i in 1:2) {
    header<-names(filelist)[i]

    # Prepare the combined data frame with gene ID as rownames and sample ID as colname
    data<-do.call("cbind", lapply(filelist[[i]], fread, header=FALSE, select=c(3)))
    unmapped<-do.call("cbind", lapply(filelist[[i]], fread, header=FALSE, select=c(4)))
    data<-as.data.frame(data)
    unmapped<-as.data.frame(unmapped)
    temp <- fread(filelist[[i]][1],header=FALSE, select=c(1))
    rownames(data)<-temp$V1
    rownames(unmapped)<-temp$V1
    colnames(data)<-gsub("_mature.*","",basename(filelist[[i]]))
    colnames(unmapped)<-gsub("_mature.*","",basename(filelist[[i]]))

    data<-data[rownames(data)!="*",,drop=FALSE]
    unmapped<-unmapped[rownames(unmapped)=="*",,drop=FALSE]

    # Write the summary table of unmapped reads
    write.table(unmapped,file=paste(header,"_unmapped_read_counts.txt",sep=""),sep='\t',quote=FALSE)

    # Remove genes with 0 reads in all samples
    row_sub = apply(data, 1, function(row) all(row ==0 ))
    # Only subset if at least one sample is remaining
    nr_keep <- sum(row_sub)
    if (nr_keep > 0){
        data<-data[!row_sub,, drop=FALSE]
    }
    #Also check for colSums > 0, otherwise DGEList will fail if samples have entirely colSum == 0 #Fixes #134
    drop_colsum_zero <- (colSums(data, na.rm=T) != 0) # T if colSum is not 0, F otherwise
    data <- data[, drop_colsum_zero] # all the non-zero columns

    write.csv(t(data),file=paste(header,"_counts.csv",sep=""))

    # Normalization
    dataDGE<-DGEList(counts=data,genes=rownames(data))
    o <- order(rowSums(dataDGE$counts), decreasing=TRUE)
    dataDGE <- dataDGE[o,]

    # Save log10(TPM)
    tpm = cpm(dataDGE, normalized.lib.sizes=F, log = F, prior.count = 0.001)
    tpm = tpm + 0.001
    tpm = log10(tpm)
    ttpm = t(tpm)
    write.table(ttpm,file=paste(header,"_logtpm.txt",sep=""),sep='\t',quote=FALSE)
    write.csv(ttpm,file=paste(header,"_logtpm.csv",sep=""))

    # TMM
    dataNorm <- calcNormFactors(dataDGE)

    # Print normalized read counts to file
    dataNorm_df<-as.data.frame(cpm(dataNorm))
    write.table(dataNorm_df,file=paste(header,"_normalized_CPM.txt",sep=""),sep='\t',quote=FALSE)

    if (length(filelist[[1]]) > 1){ # with more than 1 sample
        # Print heatmap based on normalized read counts
        pdf(paste(header,"_CPM_heatmap.pdf",sep=""))
        heatmap.2(cpm(dataNorm),col=redgreen(100),key=TRUE,scale="row",density.info="none",trace="none")
        dev.off()
    }

    # Make MDS plot (only perform with 3 or more samples)
    if (ncol(dataNorm$counts) > 2){
        pdf(paste(header,"_edgeR_MDS_plot.pdf",sep=""))
        MDSdata <- plotMDS(dataNorm)
        dev.off()

        # Print distance matrix to file
        write.table(MDSdata$distance.matrix, paste(header,"_edgeR_MDS_distance_matrix.txt",sep=""), quote=FALSE, sep="\t")

        # Print plot x,y co-ordinates to file
        MDSxy = data.frame(x=MDSdata$x, y=MDSdata$y)
        colnames(MDSxy) = c(paste(MDSdata$axislabel, '1'), paste(MDSdata$axislabel, '2'))

        write.table(MDSxy, paste(header,"_edgeR_MDS_plot_coordinates.txt",sep=""), quote=FALSE, sep="\t")

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
        plot(hmap$rowDendrogram, main="Sample Dendrogram")
        dev.off()

        # Write clustered distance values to file
        write.table(hmap$carpet, paste(header,"_log2CPM_sample_distances.txt",sep=""), quote=FALSE, sep="\t")
    } else {
    warning("Not enough samples to create an MDS plot. At least 3 samples are required.")
    }
}

file.create("corr.done")

# Print sessioninfo to standard out
print("Sample correlation info:")
sessionInfo()
