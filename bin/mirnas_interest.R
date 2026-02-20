#!/usr/bin/env Rscript

# -------------------------------
# Written by Karla Ruiz and released under the MIT license. 
# Human Technopole. National Facility for Data Handling and Analysis - IU2 OMICS
#
# This script performs the heatmap and targets identification of a given list of miRNAs of interest.
# It requires the DGE_miRNAs.RData file and the miRNAs/genes of interest file as input.
# It generates heatmaps for the miRNAs of interest and identifies their targets using the multiMiR package.
# -------------------------------



#----- Load libraries ----
suppressMessages(library(multiMiR))
suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(optparse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gplots))

#----- Input files -----
# Parse the command-line arguments
option_list <- list(make_option(c("-h", "--help"), action = "store_true", default = FALSE, help = "Show this help message and exit"))
parser <- OptionParser(usage = "%prog rdata mirnaslist organism", 
                       prog = "mirnas_interest_targets",
                       description = "-----------------------------------------------------------
            Description: This script performs the heatmap and targets identification of a given list of miRNAs of interest.
  Required input:
     1. DGE_miRNAs.RData: RData object that contains the results of the differential expression analysis.
     2. miRNAs/genes of interest file: Path to miRNAs of interest file (txt). If there are miRNAs of interest, they should be indicated in the first column. 
     If there are genes of interest, they should be indicated in the second column. File is in a txt format and should have column names.
     3. organism: organism of study ('human' or 'mouse').
                                  -----------------------------------------------------------")
arguments <- parse_args(parser, args <- commandArgs(trailingOnly = TRUE), positional_arguments = TRUE)

# Check if required arguments are provided
if (length(arguments$args) < 3) {
  stop("Error: DE_miRNAs.RData, miRNAs/genes of interest file and organism are required.")
}

# Use the arguments
rdata_file <- arguments$args[1]

mir_interest_file <- arguments$args[2]
paste("miRNAs of interest file is:", mir_interest_file)

organism <- arguments$args[3]
print(paste("Organism is:", organism))

#----- ANALYSIS ------
load(rdata_file)
mir_interest.df <- read.table(mir_interest_file, header = TRUE, sep = "\t")

# Extract and format the miRNAs of interest
mir_interest <- as.character(mir_interest.df[[1]])
mir_interest <- gsub("\\.", "-", mir_interest)

#----- Heatmap overall mirnas interest -----
dge<- dge_res[[1]]$dge
logcounts <- dge_res[[1]]$logcounts
rownames(logcounts) <- tolower(rownames(logcounts))

mir_log_both <- intersect(tolower(rownames(logcounts)), tolower(mir_interest))
mir_log_both_lcpm <- logcounts[mir_log_both , ]

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- RColorBrewer::brewer.pal(8, "Set1")[group]

png("heatmap_mirnas_interest.png", height = 7, width = 7, res = 600, units = "in")
par(mfrow=c(4, 4))
heatmap.2(mir_log_both_lcpm, col=rev(morecols(50)), 
  trace="none", main=paste("miRNAs of interest"),
  ColSideColors=col.cell,
  scale="row", margins =c(10,10),
  cexCol = 1, cexRow = 1
)
dev.off()
      
#pdf
pdf("heatmap_mirnas_interest.pdf", height = 7, width = 7)
par(mfrow=c(4, 4))
heatmap.2(mir_log_both_lcpm, col=rev(morecols(50)), trace="none", 
  main=paste("miRNAs of interest"),
  ColSideColors=col.cell,
  scale="row", margins =c(10,10),
  cexCol = 0.8, cexRow = 1
)
dev.off()

#-----Heatmap of miRNAs of interest -----
for (contrast_name in names(dge_res)) {
  # Subset log counts for miRNAs of interest
  plot_title <- gsub("_", " ", contrast_name)
  plot_title<- tools::toTitleCase(plot_title)
  
  logcounts <- dge_res[[contrast_name]]$logcounts
  num<- intersect(tolower(rownames(logcounts)), tolower(mir_interest))
  print(paste("miRNAs of interest found for", contrast_name, ":", length(num)))
  
  if (length(num) > 0) {
    mir_interest_lcpm <- logcounts[tolower(rownames(logcounts)) %in% tolower(mir_interest), , drop = FALSE]
    
    #Subset for the contrasts:
    contr_heatmap <- unlist(strsplit(contrast_name, "_vs_"))
    samples_to_keep <-dge$samples %>%
      filter(group %in% contr_heatmap)
      
    mir_interest_lcpm <- mir_interest_lcpm[, colnames(mir_interest_lcpm) %in% samples_to_keep$samples]
    all(colnames(mir_interest_lcpm) %in% samples_to_keep)      
    
    mypalette <- brewer.pal(11,"RdYlBu")
    morecols <- colorRampPalette(mypalette)
    col.cell <- RColorBrewer::brewer.pal(8, "Set1")[samples_to_keep$group]
    
    # Create PNG heatmap
    png(paste0("heatmap_mirnas_interest_", contrast_name, ".png"), height = 7, width = 10, res = 600, units = "in")
    par(mfrow = c(4, 4))
    heatmap.2(mir_interest_lcpm,
              col=rev(morecols(50)), 
              trace = "none", 
              main= paste("miRNAs of interest ", plot_title),
              ColSideColors = col.cell, 
              scale = "row",
              margins = c(10, 10),
              cexCol = 1, cexRow = 1
    )
    dev.off()
    
    # Create PDF heatmap
    pdf(paste0("heatmap_mirnas_interest_", contrast_name, ".pdf"), height = 7, width = 10)
    par(mfrow = c(4, 4))
    heatmap.2(mir_interest_lcpm,
              col=rev(morecols(50)), 
              trace = "none", 
              main= paste("miRNAs of interest ", plot_title),
              ColSideColors = col.cell, 
              scale = "row",
              margins = c(10, 10),
              cexCol = 1, cexRow = 1
    )
    dev.off()

  } else {
    
    print(paste("No miRNAs of interest were found for", contrast_name))
  }
}

#----- Targets identification -----
# miRNAs of interest
if (!is.null(mir_interest)) {
  mir_interest_res <- get_multimir(mirna = mir_interest, summary = TRUE)
  mir_interest_res.df <- as.data.frame(unique(mir_interest_res@data))
  mir_interest_res_summary <- as.data.frame(mir_interest_res@summary)
  if (!nrow(mir_interest_res.df) == 0 && !ncol(mir_interest_res.df) == 0) {
      write.xlsx(file = 'mirnas_interest_targets.xlsx', mir_interest_res.df, rowNames = FALSE)
      write.xlsx(file = 'mirnas_interest_targets_summary.xlsx', mir_interest_res_summary, rowNames = FALSE)
    } else {
      print('No results available in database.')
    } 
} else {
  print("No miRNAs of interest were provided.")
}

# Genes of interest
if (ncol(mir_interest.df) >= 2) {
  genes_interest <- as.character(mir_interest.df[[2]])
  mir_targ_genes <- get_multimir(org = organism, target = genes_interest, table = 'predicted',
                                 summary = TRUE, predicted.cutoff.type = 'n', predicted.cutoff = 500000)
  
  mir_targ_genes.df <- as.data.frame(unique(mir_targ_genes@data))
  mir_targ_genes_summary <- as.data.frame(mir_targ_genes@summary)
  
  if (!nrow(mir_targ_genes.df) == 0 && !ncol(mir_targ_genes.df) == 0) {
      write.xlsx(file = 'genes_interest_targets.xlsx', mir_targ_genes.df, rowNames = FALSE)
      write.xlsx(file = 'genes_interest_targets_summary.xlsx', mir_targ_genes_summary, rowNames = FALSE)
    } else {
      print('No results available in database.')
    }   
} else {
  print("No genes of interest were provided.")
}

