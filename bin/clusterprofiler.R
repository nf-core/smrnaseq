#!/usr/bin/env Rscript

# -------------------------------
# Written by Karla Ruiz and released under the MIT license. 
# Human Technopole. National Facility for Data Handling and Analysis - IU2 OMICS
#
# This script performs the enrichment analysis of the miRNAs' targets from the DGE analysis on miRNA data.
# It requires the miRNAs_targets.RData file and the organism name as input.
# It performs Gene Ontology (GO) and KEGG enrichment analysis using the clusterProfiler package.
# It saves the results in an Excel file and generates plots for the enrichment analysis.
# -------------------------------

###------ENRICHMENT-----
#BiocManager::install(c("clusterProfiler", "AnnotationHub", "ggarchery"))
suppressMessages(library(enrichplot))
suppressMessages(library(ggplot2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(AnnotationHub))
suppressMessages(library(ggarchery))
suppressMessages(library(openxlsx))
suppressMessages(library(optparse))

#----- Input files -----
# Parse the command-line arguments
parser <- OptionParser(usage = "%prog [options] rdata organism",
                       prog = "mir_targets",
                       description = "-----------------------------------------------------------
  Description: This script performs the enrichment analysis of the miRNAs' targets from the DGE analysis on miRNA data.
  Required input:
     1. miRNAs_targets.RData: RData object that contains the results of targets identification.
     2. organism: organism of study
                        -----------------------------------------------------------")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly = TRUE), positional_arguments = TRUE)

# Check if at least the required arguments are provided
if (length(args) < 2) {
  stop("Error: miRNAs_targets.RData and organism are required.")
}

# Use the arguments
rdata_file <- arguments$args[1]

organism <- arguments$args[2]
print(paste("Organism is:", organism))

#----- FUNCTIONS------
#Search for the organism database
species_codes <- c("mmu", "hsa", "rno", "mml", "dre", "dme","ath","cfa","bta","cel","gga","ptr","xtr")
  
orgdb<- c( "org.Mmu.eg.db.sqlite", #mmu
   "org.Hs.eg.db.sqlite", #hsa
   "org.Rn.eg.db.sqlite", #rno
   "org.Mm.eg.db.sqlite", #mml  
   "org.Dr.eg.db.sqlite", #dre
   "org.Dm.eg.db.sqlite", #dme
   "org.At.tair.db.sqlite", #ath
   "org.Cf.eg.db.sqlite", #cfa
   "org.Bt.eg.db.sqlite", #bta
   "org.Ce.eg.db.sqlite", #cel
   "org.Gg.eg.db.sqlite", #gga
   "org.Pt.eg.db.sqlite", #ptr
   "org.Xl.eg.db.sqlite" #xtr
   )

hub <- AnnotationHub()

orgdb_ann <- function(organism) {
  index <- match(organism, species_codes) # Encuentra la posición del organismo en el vector species_codes
  if (!is.na(index)) {
    db_name <- orgdb[index] # Obtener el nombre de la base de datos
    query_result <- query(hub, db_name) # Buscar en AnnotationHub
    
    if (length(query_result) > 0) {
      return(query_result[[length(query_result)]]) # Devolver el último resultado encontrado
    } else {
      stop("Database not found.")
    }
  } else {
    stop("Organism not found.")
  }
}

orgdb <- orgdb_ann(organism)

#Save database info
orgdb <- orgdb_ann(organism)
info_db <- metadata(orgdb)

write.csv(info_db, "./enrichment/orgdb_metadata.csv", row.names = FALSE)


# Enrich function
perform_enrichment <- function(input, orgdb, ontology, contrast_name, datatype, enrichment_dir) {
  ego <- enrichGO(
    gene          = input$target_symbol,
    keyType       = "SYMBOL",
    OrgDb         = orgdb,
    ont           = ontology,
    pAdjustMethod = "bonferroni",
    pvalueCutoff  = 0.1,
    readable      = TRUE)
  
  if (nrow(ego@result) > 0) {
    ego2 <- pairwise_termsim(ego)
    
    # Plots dir
    plots_dir <- file.path(enrichment_dir, "plots")
    dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Save results
    output_file <- file.path(enrichment_dir, paste0("GO_", ontology, "_", datatype, ".xlsx"))
    write.xlsx(ego@result, file = output_file)
    
    # Create plots
    go_barplot <- ego %>%
      mutate(qscore = -log10(p.adjust)) %>%
      barplot(x = "qscore") + ggtitle(paste("GO", ontology, contrast_name, datatype))
    
    go_dotplot <- dotplot(ego, showCategory = 20) + ggtitle(paste("GO", ontology, contrast_name, datatype))
    go_treeplot <- treeplot(ego2, hclust_method = "average") + ggtitle(paste("GO", ontology, contrast_name, datatype))
    go_emmapplot <- emapplot(ego2) + ggtitle(paste("GO", ontology, contrast_name, datatype))
    
    # Lista de gráficos y nombres
    plots <- list(
      barplot = go_barplot,
      dotplot = go_dotplot,
      treeplot = go_treeplot,
      emapplot = go_emmapplot
    )
    
    # Export plots
    for (plot_name in names(plots)) {
      ggsave(filename = file.path(plots_dir, paste0(ontology, "_", plot_name, "_", datatype, ".png")), 
             plot = plots[[plot_name]], width = 8, height = 6, dpi = 300)
    }
    
    return(list(ego = ego, ego2 = ego2, plots = plots))
  } else {
    message(paste("No significant data", ontology, datatype))
    return(NULL)
  }
}

# ----- Analysis ----Procesar cada contraste en DGE
ontologies <- c("BP", "CC", "MF")

results_enrich_targets <- list()
for (contrast_name in names(dge_res)) {
  print(paste("Processing contrast:", contrast_name))
  
  # Folder results
  dge_ctr_folder <- file.path("./enrichment", contrast_name)
  dir.create(dge_ctr_folder, showWarnings = FALSE, recursive = TRUE)
  
  enrichr_input <- results_mirna_targets[[contrast_name]]$enrichr_input
  for (datatype in names(enrichr_input)) {
    input <- enrichr_input[[datatype]]
    
    enrichment_results <- file.path(dge_ctr_folder, datatype)
    dir.create(enrichment_results, showWarnings = FALSE, recursive = TRUE)
    
    # GO analysis
    print(paste("Processing:", ont, "for", contrast_name, datatype))
    enrichment_results_list <- list()
    for (ont in ontologies) {
      enrichment_results_list[[ont]] <- perform_enrichment(input, orgdb, ont, contrast_name, datatype, enrichment_results)
    }
    
    #KEGG analysis
    print(paste("Processing KEGG for", contrast_name, datatype))
    kegg <- enrichKEGG(gene         = input$target_entrez,
                     organism     = organism,
                     pvalueCutoff = 0.1)
    
    if (nrow(kegg@result) > 0) {
      write.xlsx(kegg@result, file = file.path(enrichment_results, paste0("KEGG_", datatype, ".xlsx")))
      
      # Plots
      barplot_kegg<- barplot(kegg) + ggtitle(paste0("KEGG ", contrast_name, " ", datatype))
      dotplot_kegg<- dotplot(kegg) + ggtitle(paste0("KEGG ", contrast_name, " ", datatype))
      
      ggsave(file.path(enrichment_results, paste0("KEGG_dotplot_", datatype, ".png")), dotplot_kegg, height = 7, width = 7, dpi = 600)
      ggsave(file.path(enrichment_results, paste0("KEGG_barplot_", datatype, ".png")), barplot_kegg, height = 7, width = 7, dpi = 600)
    }
    
    results_enrich_targets[[contrast_name]][[datatype]] <- list(GO = enrichment_results_list, KEGG = kegg)
  }
}

save(results_enrich_targets, file = "enrichr.RData")
