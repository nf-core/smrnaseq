#!/usr/bin/env Rscript

###------ENRICHMENT-----
suppressMessages(library(enrichR))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(optparse))

#----- Input files -----
# Parse the command-line arguments
parser <- OptionParser(usage = "%prog [options] rdata organism",
                       prog = "mir_targets",
                       description = "-----------------------------------------------------------
  Description: This script performs the identification of miRNAs targets from the DGE analysis on miRNA data.
  Required input:
     1. miRNAs_targets.RData: RData object that contains the results of targets identification.
     2. organism: organism of study ('human' or 'mouse')
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

#Connect to enrichR site
setEnrichrSite("Enrichr")

# Conditional statement to choose the appropriate databases
if (organism %in% c("hsa", "human")) {
  dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022", "KEGG_2021_Human")
  plot_info <- list(
    "GO_Biological_Process_2023" = list(title = "Biological Process", suffix = "GO_Biol_Proc"),
    "GO_Molecular_Function_2023" = list(title = "Molecular Function", suffix = "GO_Mol_Fun"),
    "GO_Cellular_Component_2023" = list(title = "Cellular Component", suffix = "GO_Cel_Comp"),
    "KEGG_2021_Human" = list(title = "KEGG", suffix = "KEGG"),
    "Reactome_2022" = list(title = "Reactome", suffix = "Reactome")
  )

} else if (organism %in% c("mmu", "mouse")) {
  dbs <- c("Mouse_Gene_Atlas", "WikiPathways_2024_Mouse", "KEGG_2019_Mouse")
  plot_info <- list(
    "Mouse_Gene_Atlas" = list(title = "Mouse Gene Atlas", suffix = "Mouse_Gen_Atl"),
    "WikiPathways_2024_Mouse" = list(title = "WikiPathways", suffix = "WikiPathways"),
    "KEGG_2019_Mouse" = list(title = "KEGG", suffix = "KEGG")
  )

} else {
  stop("Invalid organism specified. Please choose either 'human' or 'mouse'.")
}

#-----ANALYSIS-----
#Verificar que el input coincida incluyendo los contrastes
load (rdata_file)

results_enrich_targets <- list()
for (contrast_name in names(dge_res)) {
  print(paste("Processing contrast:", contrast_name))
  dge_ctr_folder <- paste0("./enrich_", contrast_name)

  if (dir.exists(dge_ctr_folder)) {
    #cat("Folder created successfully:", dge_ctr_folder, "\n")
  } else {
    dir.create(dge_ctr_folder)
  }

  enrichr_input <- results_mirna_targets[[contrast_name]]$enrichr_input

  for (datatype in names(enrichr_input)) {
    input <- enrichr_input[[datatype]]
    enrichment_results <- paste0(dge_ctr_folder, '/', datatype)

  if (!dir.exists(enrichment_results)) {
      dir.create(enrichment_results)
  }

  # Perform enrichment analysis
  enriched <- enrichr(input$target_symbol, dbs)
  write.xlsx(enriched, file = file.path(enrichment_results, paste0("enrichment_", datatype, ".xlsx")))

  #Plot
  for (name in names(plot_info)) {
    if (name %in% names(enriched)) {
      goplot <- plotEnrich(enriched[[name]], showTerms = 20, numChar = 40,
                           y = "Count", orderBy = "P.value", title = plot_info[[name]]$title)
      ggsave(file.path(enrichment_results, paste0(plot_info[[name]]$suffix, "_", datatype, ".png")),
             goplot, height = 7, width = 7, dpi = 600)

      print(paste(name, "plot created"))
    } else {
      warning(paste(name, "not found in the results."))
    }
    }
  }
}

save(results_enrich_targets, file = "enrichr.RData")
