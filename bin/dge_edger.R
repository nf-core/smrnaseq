#!/usr/bin/env Rscript

# -------------------------------
# Written by Karla Ruiz and released under the MIT license.
# Human Technopole. National Facility for Data Handling and Analysis - IU2 OMICS
#
# This script performs Differential Expression analysis on miRNA data using edgeR.
# It requires the miRNA counts file, metadata file, and contrast as input.
# It performs normalization, filtering, and statistical testing using edgeR.
# It saves the results in a text file and generates plots for the analysis.
# -------------------------------


#----- Load libraries ----
suppressMessages(library(tidyverse))
suppressMessages(library(tibble))
suppressMessages(library(edgeR))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gplots))
suppressMessages(library(optparse))
suppressMessages(library(openxlsx))


#----- Input files -----
option_list <- list(
  make_option(c("--LFC"), action="store", type = "double", default = 0.5, help = "Log Fold Change cutoff (default: 0.5)"),
  make_option(c("--FDR"), action="store", type = "double", default = 0.05, help = "Adjusted p-value cutoff (default: 0.05)")
)

# Parse the command-line arguments
parser <- OptionParser(usage = "%prog [options] counts metadata contrast",
                       option_list = option_list, prog = "dge_edger",
                       description = "-----------------------------------------------------------
  Description: This script performs Differential Expression analysis on miRNA data using edgeR.
  Required input:
     1. miRNA counts file: Path to the miRNA counts data file (CSV format).
     2. metadata file: Path to the metadata file (CSV format).
     3. contrast: is a string of format test:reference, with words separated by colon (:), where:
           - 'test' is the name of the factor level to be tested (e.g. treat). If multiple levels should be aggregated, they must be separated by '/' with no blank spaces (e.g. treat1/treat2/treat3)
           - 'reference' is the name of the factor level to be used as reference (e.g. ctrl). If multiple levels should be aggregated, they must be separated by '/' with no blank spaces (e.g. ctrl1/ctrl2).
                        -----------------------------------------------------------")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly = TRUE), positional_arguments = TRUE)
opt <- arguments$options

# Check if required arguments are provided
if (length(arguments$args) < 4) {
  stop("Error: At least four arguments (miRNA counts file, metadata file, and contrast) are required.")
}

# Use the arguments
mircounts_file <- arguments$args[1]
print(paste("miRNA counts file is:", mircounts_file))

metadata_file <- arguments$args[2]
print(paste("Metadata file is:", metadata_file))

contrast <- arguments$args[3]
print(paste("Contrast is:", contrast))

# Optional arguments
LFC.cutoff <- opt$LFC
print(paste("LFC cutoff is:", LFC.cutoff))

padj.cutoff <- opt$FDR
print(paste("Adjusted p-value cutoff (FDR) is:", padj.cutoff))

######

#----FUNCTIONS----
#Import miRNA counts and metadata file
import_data <- function(mircounts_file, metadata_file) {
  # Read the raw count table
  rawCountTable <- read.delim(mircounts_file)
  metadata <- read.delim(metadata_file, sep= ',', row.names = 1)

  # Transform the data
  rawCountTable <- rawCountTable %>%
    column_to_rownames(var = "miRNA") %>%
    as.data.frame()

  rownames(rawCountTable) <- gsub("\\.", "-", rownames(rawCountTable))
  colnames(rawCountTable) <- gsub("\\.", "-", colnames(rawCountTable))

  column_order <- rownames(metadata)
  rawCountTable <- rawCountTable[, column_order]

  return(list(rawCountTable = rawCountTable, metadata= metadata))
  #print("miRNA data has been imported")
}

# DGE
perform_dge_analysis <- function(rawCountTable, metadata) {
  # Create design matrix
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)

  parts <- strsplit(contrast, ":")[[1]]
  tests <- unlist(strsplit(parts[1], "/"))
  controls <- unlist(strsplit(parts[2], "/"))
  x <- sapply(1:length(tests), function(i) {
    if (i <= length(controls)) {
      paste0(tests[i], "-", controls[i])
    } else {
      paste0(tests[i], "-")
    }
  })

  x <- as.character(x)
  print(paste("Contrasts:", x))

  # Create DGEList object
  dge <- DGEList(counts = rawCountTable, group = group, samples = samples)

  # Remove low expressed genes
  keep <- filterByExpr(dge, min.count = 5)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  # Normalization
  dge <- calcNormFactors(dge, method="TMM")

  # Estimate dispersions
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)

  #GENERAL PLOTS
  ##MDS/PCA
  #png
  if (ncol(dge$counts) > 2){
    suppressWarnings({
      png("MDS.png", height = 7, width = 7, res = 600, units = "in")
      plotMDS(dge, col=as.numeric(dge$samples$group), cex = 0.8)
      legend("bottomright",legend=levels(dge$samples$group), col=1:length(unique(dge$samples$group)), pch=20, horiz = T, cex = 0.5)
      dev.off()

      #pdf
      pdf("MDS.pdf", height = 7, width = 7)
      plotMDS(dge, col=as.numeric(dge$samples$group), cex = 0.8)
      legend("bottomright",legend=levels(dge$samples$group), col=1:length(unique(dge$samples$group)), pch=20, horiz = T, cex = 0.5)
      dev.off()

      ## Heatmap
      logcounts <- cpm(dge, log = TRUE)
      cpm_counts<-as.data.frame(cpm(dge))

      # Calculate variance for each gene
      var_genes <- apply(logcounts, 1, var)

      # Select top variable genes
      select_var <- names(sort(var_genes, decreasing = TRUE))[1:10]

      # Subset log counts for highly variable genes
      highly_variable_lcpm <- logcounts[select_var, ]
      #write.table(cpm_counts,file=file.path(dge_ctr_folder, paste0("normalized_CPM_", contrast_name, ".txt")), row.names = T)

      mypalette <- brewer.pal(11,"RdYlBu")
      morecols <- colorRampPalette(mypalette)
      #col.cell <- c("purple","orange")[group]
      col.cell <- RColorBrewer::brewer.pal(8, "Set1")[group]

      #png
      png("heatmap.png", height = 7, width = 7, res = 600, units = "in")
      par(mfrow=c(4, 4))
      heatmap.2(highly_variable_lcpm,
                col=rev(morecols(50)),
                trace="none",
                main=paste("Top variable genes"),
                ColSideColors=col.cell,
                scale="row",
                margins =c(10,10),
                cexCol = 1, cexRow = 1
      )
      dev.off()

      #pdf
      pdf("heatmap.pdf", height = 7, width = 7)
      par(mfrow=c(4, 4))
      heatmap.2(highly_variable_lcpm,
                col=rev(morecols(50)),
                trace="none",
                main=paste("Top variable genes"),
                ColSideColors=col.cell,
                scale="row",
                margins =c(10,10),
                cexCol = 0.8, cexRow = 1
      )
      dev.off()

    })
  } else {
    warning("Not enough samples.")

  }

  #CONTRASTS ANALYSES -----
  my.contrasts<- makeContrasts(contrasts =x, levels=design)
  colnames(my.contrasts)<- gsub("-", "_vs_", colnames(my.contrasts))

  results_list <- list()
  for (i in 1:ncol(my.contrasts)) {
    contrast_name <- colnames(my.contrasts)[i]
    print(paste("Processing contrast:", contrast_name))
    contr <- my.contrasts[, i]  # Get the contrast for the current iteration

    dge_ctr_folder <- paste0("./DE_", contrast_name)

    if (dir.exists(dge_ctr_folder)) {
      #cat("Folder created successfully:", dge_ctr_folder, "\n")
    } else {
      dir.create(dge_ctr_folder)
    }

    # Fit the model for the current contrast
    fit <- glmQLFit(dge, design)
    glmtest <- glmQLFTest(fit, contrast = contr)

    # Perform the likelihood ratio test
    glmfit <- glmFit(dge, design)
    lrt <- glmLRT(glmfit, contrast = contr)

    # Summarize the results
    sum_res_glm <- summary(decideTests(glmtest))
    sum_res_lrt <- summary(decideTests(lrt))

    # Return top tags
    res <- topTags(lrt, adjust.method = 'fdr', n = nrow(lrt$table))

    # Convert the DGE results to a data frame
    res_DGE_lab <- as.data.frame(res, row.names = rownames(res))

    # Initialize the 'diffexpressed' column
    res_DGE_lab$diffexpressed <- "NO"
    res_DGE_lab$diffexpressed[res_DGE_lab$logFC >= LFC.cutoff & res_DGE_lab$FDR <= padj.cutoff] <- "UP"
    res_DGE_lab$diffexpressed[res_DGE_lab$logFC <= -LFC.cutoff & res_DGE_lab$FDR <= padj.cutoff] <- "DOWN"

    # Add miRNA names and reorder columns
    res_DGE_lab <- res_DGE_lab %>%
      mutate(miRNA = rownames(res_DGE_lab)) %>%
      select(miRNA, everything()) %>%
      as.data.frame() %>%
      { rownames(.) <- NULL; . }

    # Getting top variable miRNAs
    # Set topgenes based on the number of rows - config file
    if (dim(res_DGE_lab)[1] < 100) {
      topgenes <- dim(res_DGE_lab)[1]  # Use the maximum number of rows if less than 50
    } else {
      topgenes <- 100  # Otherwise, use 50
    }

    print(paste("Top variable genes in", contrast_name, topgenes))

    #Export a summary
    output_file <-  file.path(dge_ctr_folder, paste0("DGE_analysis_summary_", contrast_name, ".txt"))
    file_conn <- file(output_file, open = "wt")
    cat(paste("Differential expression analysis performed with miRNAs data:", contrast_name,"\n",
              "\n",
              "Thresholds used: LFC.cutoff", LFC.cutoff, " and padj.cutoff:", padj.cutoff,"\n",
              "Total significant miRNAs identified:", length(res_DGE_lab$miRNA), ", from which: \n"), file = file_conn)

    for (i in 1:nrow(sum_res_lrt)) {
      cat(paste(rownames(sum_res_lrt)[i], sum_res_lrt[i, ], sep = "\t"), "\n", file = file_conn)
    }
    close(file_conn)

    # Write results to a CSV file, including the contrast name
    write.table(res_DGE_lab, file = file.path(dge_ctr_folder, paste0("DGE_", contrast_name, ".txt")), row.names = FALSE)
    write.xlsx(res_DGE_lab, file = file.path(dge_ctr_folder, paste0("DGE_", contrast_name, ".xlsx")), rowNames = FALSE)

    #----- PLOTS -----
    #folder
    dge_png_plots_folder <- paste0(dge_ctr_folder, "/plots/png")
    dge_pdf_plots_folder <- paste0(dge_ctr_folder, "/plots/pdf")

    if (dir.exists(dge_png_plots_folder)) {
      #cat("Folder created successfully:", dge_png_plots_folder, "\n")
    } else {
      dir.create(dge_png_plots_folder, recursive = TRUE)
    }

    if (dir.exists(dge_pdf_plots_folder)) {
      #cat("Folder created successfully:", dge_pdf_plots_folder, "\n")
    } else {
      dir.create(dge_pdf_plots_folder, recursive = TRUE)
    }

    plot_title <- gsub("_", " ", contrast_name)
    #plot_title<- gsub("([a-z])([0-9])", "\\1 \\2", plot_title)
    plot_title<- tools::toTitleCase(plot_title)

    ## Volcanoplot
    vplot<- EnhancedVolcano(res$table,
                            lab = rownames(res),
                            x = "logFC",
                            y = "FDR",
                            pCutoffCol = 'FDR',
                            labSize = 3.0,
                            pCutoff = padj.cutoff,
                            #ylim = c(4, 19),
                            FCcutoff = LFC.cutoff,
                            subtitle = "",
                            title = plot_title,
                            legendPosition = "right",
                            legendLabSize = 10,
                            legendIconSize = 3)
    ggsave(file.path(dge_png_plots_folder, paste0("volcanoplot_", contrast_name, ".png")), vplot, height = 7, width = 10, dpi = 600)
    ggsave(file.path(dge_pdf_plots_folder, paste0("volcanoplot_", contrast_name, ".pdf")), vplot, height = 7, width = 10, dpi = 600)

    ## Heatmap
    # Select top variable genes
    select_var_mir <- names(sort(var_genes, decreasing = TRUE))[1:topgenes]

    # Subset log counts for highly variable genes
    highly_variable_mir_lcpm <- logcounts[select_var_mir, ]

    #Subset for the contrasts:
    contr_heatmap <- unlist(strsplit(contrast_name, "_vs_"))
    samples_to_keep <-dge$samples %>%
      filter(group %in% contr_heatmap)

    highly_variable_mir_lcpm <- highly_variable_mir_lcpm[, colnames(highly_variable_mir_lcpm) %in% samples_to_keep$samples]
    all(colnames(highly_variable_mir_lcpm) %in% samples_to_keep)

    #write.table(cpm_counts,file=file.path(dge_ctr_folder, paste0("normalized_CPM_", contrast_name, ".txt")), row.names = T)

    mypalette <- brewer.pal(11,"RdYlBu")
    morecols <- colorRampPalette(mypalette)
    col.cell <- RColorBrewer::brewer.pal(8, "Set1")[samples_to_keep$group]

    #png
    png(file.path(dge_png_plots_folder, paste0("heatmap_", contrast_name, ".png")), height = 10, width = 10, res = 600, units = "in")
    par(mfrow=c(4, 4))
    heatmap.2(highly_variable_mir_lcpm,
              col=rev(morecols(50)),
              trace="none",
              main=paste("Top", topgenes, "variable genes ", plot_title),
              ColSideColors=col.cell,
              scale="row",
              margins =c(10,10),
              cexCol = 1, colRow = 1
    )
    dev.off()

    #pdf
    pdf(file.path(dge_pdf_plots_folder, paste0("heatmap_", contrast_name, ".pdf")), height = 10, width = 10)
    par(mfrow=c(4, 4))
    heatmap.2(highly_variable_mir_lcpm,
              col=rev(morecols(50)),
              trace="none",
              main=paste("Top", topgenes, "variable genes ", plot_title),
              ColSideColors=col.cell,
              scale="row",
              margins =c(10,10),
              cexCol = 1, colRow = 1
    )
    dev.off()


    # Store the results in the results list
    results_list[[contrast_name]] <- list(
      res = res,
      dge = dge,
      lrt = lrt,
      res_DGE_lab = res_DGE_lab,
      sum_res_glm = sum_res_glm,
      sum_res_lrt = sum_res_lrt,
      highly_variable_lcpm = highly_variable_lcpm,
      selected_genes = select_var,
      logcounts = logcounts,
      cpm_counts = cpm_counts,
      dge_png_plots_folder = dge_png_plots_folder,
      dge_pdf_plots_folder = dge_pdf_plots_folder
    )
  }

  return(results_list = results_list)
}


# ---- ANALYSIS -----

data_mirnas <- import_data(mircounts_file, metadata_file)
mircounts <- data_mirnas$rawCountTable
metadata <- data_mirnas$metadata

# Create factors for group and samples
contrast<- tolower(contrast)
conditions <-unlist(strsplit(contrast, "[:/]"))

#Buscar el contraste en la tabla de metadatos
col_name <- c()
for (condition in conditions) {
  found <- sapply(metadata, function(col) any(grepl(condition, tolower(col))))
  if (any(found)) {
    col_name <- names(metadata)[which(found)]
  }
}

# Crear el factor
group <- factor(tolower(metadata[[col_name]]))
samples <- factor(rownames(metadata))

## DGE
dge_res <- perform_dge_analysis(mircounts, metadata)

save(mircounts, metadata, conditions, dge_res, group, samples, file="DGE_miRNAs.RData")
