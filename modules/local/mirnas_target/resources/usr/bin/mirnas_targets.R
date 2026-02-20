#!/usr/bin/env Rscript

#----- Load libraries ----
suppressMessages(library(multiMiR))
suppressMessages(library(tidyverse))
suppressMessages(library(openxlsx))
suppressMessages(library(optparse))

#----- Input files -----
# Parse the command-line arguments
option_list <- list(
  make_option(c("--drug"), action="store", type = "character", default = NULL, help = "Drug of interest (e.g: 'cisplatin'). If multiple drugs should be aggregated, they must be separated by '/' with no blank spaces (e.g. cisplatin/oxaliplatin)."),
  make_option(c("--disease"), action="store", type = "character", default = NULL, help = "Disease of interest (e.g: 'lung cancer'). If multiple diseases should be aggregated, they must be separated by '/' with no blank spaces (e.g. lung cancer/bladder cancer).")
)

parser <- OptionParser(usage = "%prog [options] rdata mirtrace_species",
                       option_list = option_list, prog = "mir_targets",
                       description = "-----------------------------------------------------------
  Description: This script performs the identification of miRNAs targets from the DGE analysis on miRNA data.
  Required input:
     1. miRNAs RData: RData object that contains the results of the differential expression analysis.
     2. organism: organism of study ('human' or 'mouse').
                        -----------------------------------------------------------")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly = TRUE), positional_arguments = TRUE)
opt <- arguments$options

# Check if at least the required arguments are provided
if (length(args) < 2) {
  stop("Error: RData and organism are required.")
}

# Use the arguments
rdata_file <- arguments$args[1]

organism <- arguments$args[2]
print(paste("Organism is:", organism))

# Optional arguments
drug <- opt$drug
disease <- opt$disease

#----- miRNA-Target -----
load(rdata_file)

#export databases info
db.info<- multimir_dbInfo()
write.table(db.info, file="miRNA_targetDB_info.txt" , sep = "\t", row.names = FALSE, quote = FALSE)

#------ANALYSIS------
results_mirna_targets <- list()
db_list_val <- c("mirecords", "mirtarbase")
db_list_pred <- c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb", "pictar", "pita")
suppressWarnings(
for (contrast_name in names(dge_res)) {
  print(paste("Processing contrast:", contrast_name))
  dge_ctr_folder <- paste0("./targets_", contrast_name)

  if (dir.exists(dge_ctr_folder)) {
    #cat("Folder created successfully:", dge_ctr_folder, "\n")
  } else {
    dir.create(dge_ctr_folder)
  }

  res_DGE_lab <- dge_res[[contrast_name]]$res_DGE_lab
  data_res_val <- list()
  summary_res_val<- list()
  print(paste("Processing validated started", Sys.time()))
  for (val_dbs in db_list_val){
  multimir_validated <- get_multimir(org = organism,
                                     mirna = res_DGE_lab$miRNA,
                                     table = val_dbs,
                                     summary = TRUE)

  data_res_val[[val_dbs]] <- multimir_validated@data
  summary_res_val[[val_dbs]] <- multimir_validated@summary
  }
  print(paste("Processing validated db ended", Sys.time()))

  data_res_pred <- list()
  summary_res_pred<- list()
  print(paste("Processing predicted db started", Sys.time()))
  for (pred_dbs in db_list_pred){
  multimir_predicted <- get_multimir(org = organism,
                                     mirna = res_DGE_lab$miRNA,
                                     table = pred_dbs,
                                     summary = TRUE)
  data_res_pred[[pred_dbs]] <- multimir_predicted@data
  summary_res_pred[[pred_dbs]] <- multimir_predicted@summary
  }
  print(paste("Processing predicted db ended", Sys.time()))

  # Convert to data frames
  multimir_validated.df <- do.call(rbind, data_res_val)
  multimir_validated.df <- unique(multimir_validated.df)
  multimir_validated_summary <- bind_rows(summary_res_val)
  multimir_validated_summary[is.na(multimir_validated_summary)] <- 0
  multimir_validated_summary <- multimir_validated_summary %>%
    group_by(mature_mirna_acc, mature_mirna_id, target_symbol, target_entrez, target_ensembl) %>%
    summarise(
      mirecords = sum(mirecords),
      mirtarbase = sum(mirtarbase),
      .groups = 'drop'  # This prevents grouping in the final output
    ) %>%
    filter(!is.na(target_ensembl) & target_ensembl != "") %>%
    mutate(sum = (mirecords != 0) + (mirtarbase != 0))

  #multimir_validated.df <- as.data.frame(unique(multimir_validated@data))
  #multimir_validated_summary <- as.data.frame(multimir_validated@summary)

  multimir_predicted.df <- do.call(rbind, data_res_pred)
  multimir_predicted.df <- unique(multimir_predicted.df)
  multimir_predicted_summary <- bind_rows(summary_res_pred)
  multimir_predicted_summary[is.na(multimir_predicted_summary)] <- 0
  multimir_predicted_summary <- multimir_predicted_summary %>%
    group_by(mature_mirna_acc, mature_mirna_id, target_symbol, target_entrez, target_ensembl) %>%
    summarise(
      diana_microt = sum(diana_microt),
      elmmo = sum(elmmo),
      microcosm = sum(microcosm),
      miranda = sum(miranda),
      mirdb = sum(mirdb),
      pictar = sum(pictar),
      pita = sum(pita),
      .groups = 'drop'  # This prevents grouping in the final output
    ) %>%
    filter(!is.na(target_ensembl) & target_ensembl != "")  %>%
    mutate(sum = (diana_microt != 0) + (elmmo != 0) + (microcosm != 0) + (miranda != 0) +
             (mirdb != 0) + (pictar != 0) + (pita != 0))

  #multimir_predicted.df <- as.data.frame(unique(multimir_predicted@data))
  #multimir_predicted_summary <- as.data.frame(multimir_predicted@summary)

  #multimir_all.df <- full_join(multimir_validated.df, multimir_predicted.df,
  #          by = c("database", "mature_mirna_acc", "mature_mirna_id",
  #                 "target_symbol", "target_entrez", "target_ensembl", "type"))
  #multimir_all_summary <- merge(multimir_validated_summary, multimir_predicted_summary, by = "mature_mirna_acc", all = T)

  #multimir_all.df <- as.data.frame(unique(multimir_all@data))
  #multimir_all_summary <- as.data.frame(multimir_all@summary)

  # Classify by up/down regulated
  # UP
  upreg_mirnas <- res_DGE_lab %>%
    filter(diffexpressed == "UP")

  if (nrow(upreg_mirnas) > 0) {
    upreg_mir_targets <- multimir_validated_summary %>%
      filter(tolower(mature_mirna_id) %in% tolower(upreg_mirnas$miRNA))

    upreg_df <- upreg_mir_targets %>%
      group_by(mature_mirna_id) %>%
      summarise(target_genes = n(), .groups = 'drop') %>%
      mutate(diffexpressed = "UP")

    uplist <- list("Upregulated_miRNA_targets" = upreg_mir_targets, "Number of targets per miRNA" = upreg_df)

    write.xlsx(uplist, file = file.path(dge_ctr_folder, paste0('validated_UpRegulated_mirnas_targets_', contrast_name, ".xlsx")), rowNames = FALSE)
  } else {
    message("No UPREGULATED miRNAs are present.")
  }

  # DOWN
  downreg_mirnas <- res_DGE_lab %>%
    filter(diffexpressed == "DOWN")

  if (nrow(downreg_mirnas) > 0) {
    downreg_mir_targets <- multimir_validated_summary %>%
      filter(tolower(mature_mirna_id) %in% tolower(downreg_mirnas$miRNA))

    downreg_df <- downreg_mir_targets %>%
      group_by(mature_mirna_id) %>%
      summarise(target_genes = n(), .groups = 'drop') %>%
      mutate(diffexpressed = "DOWN")

    downlist <- list("Downregulated_miRNA_targets" = downreg_mir_targets, "Number of targets per miRNA" = downreg_df)

    write.xlsx(file = file.path(dge_ctr_folder, paste0('validated_DownRegulated_mirnas_targets_', contrast_name, ".xlsx")), downlist, rowNames = FALSE)

  } else {
    message("No DOWNREGULATED miRNAs are present.")
  }

  enrichr_input <- list("All_targets" = multimir_validated_summary)
  if (nrow(upreg_mirnas) > 0) {
    enrichr_input[["Upregulated"]] <- upreg_mir_targets
  }
  if (nrow(downreg_mirnas) > 0) {
    enrichr_input[["Downregulated"]] <- downreg_mir_targets
  }

  # Write results to CSV
  # Validated
  write.xlsx(file = file.path(dge_ctr_folder, paste0('validated_DGE_mirnas_targets_', contrast_name, ".xlsx")), multimir_validated.df, rowNames = FALSE)
  write.xlsx(file = file.path(dge_ctr_folder,paste0('validated_DGE_mirnas_targets_summary_', contrast_name, ".xlsx")), multimir_validated_summary, rowNames = FALSE)

  # Predicted
  write.xlsx(file = file.path(dge_ctr_folder, paste0('predicted_DGE_mirnas_targets_', contrast_name, ".xlsx")), multimir_predicted.df, rowNames = FALSE)
  write.xlsx(file = file.path(dge_ctr_folder, paste0('predicted_DGE_mirnas_targets_summary_', contrast_name, ".xlsx")), multimir_predicted_summary, rowNames = FALSE)

  # All
  #write.xlsx(file = file.path(dge_ctr_folder, paste0('all_DGE_mirnas_targets_', contrast_name, ".xlsx")), multimir_all.df, rowNames = FALSE)
  #write.xlsx(file = file.path(dge_ctr_folder, paste0('all_DGE_mirnas_targets_summary_', contrast_name, ".xlsx")), multimir_all_summary, rowNames = FALSE)

  results_mirna_targets[[contrast_name]] <- list(
    enrichr_input = enrichr_input,
    multimir_validated_summary = multimir_validated_summary,
    multimir_predicted_summary = multimir_predicted_summary
  )

  if (nrow(upreg_mirnas) > 0) {
    results_mirna_targets[[contrast_name]]$upreg_mir_targets <- upreg_mir_targets
  }

  if (nrow(downreg_mirnas) > 0) {
    results_mirna_targets[[contrast_name]]$downreg_mir_targets <- downreg_mir_targets
  }

})

# Retrieve miRNA-target interactions associated with a given drug or disease (if provided)
if (!is.null(drug)) {
  drug_split <- unlist(strsplit(drug, "/"))
  mir_drug <- get_multimir(disease.drug = drug_split, table = 'disease.drug')
  mir_drug.df <- as.data.frame(unique(mir_drug@data))

  if (!nrow(mir_drug.df) == 0 && !ncol(mir_drug.df) == 0) {
    write.xlsx(file = 'drugs_association.xlsx', mir_drug.df, rowNames = FALSE)
  } else {
    print('No results available in database.')
  }
}

if (!is.null(disease)) {
    disease_split <- unlist(strsplit(disease, "/"))
    mir_disease <- get_multimir(disease.drug = disease_split, table = 'disease.drug')
    mir_disease.df <- as.data.frame(unique(mir_disease@data))

    if (!nrow(mir_disease.df) == 0 && !ncol(mir_disease.df) == 0) {
      write.xlsx(file = 'disease_association.xlsx', mir_disease.df, rowNames = FALSE)
    } else {
      print('No results available in database.')
    }
}

save(results_mirna_targets, dge_res, file="miRNAs_targets.RData")
