#!/usr/bin/env Rscript

library(optparse)
library(tidyr)
library(vroom)
library(dplyr)

option_list <- list(
    make_option(c("--input"), type = "character", help = "Input CSV file in long format", metavar = "character"),
    make_option(c("--output"), type = "character", help = "Output CSV file in wide format", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Read CSV with vroom
long_data <- vroom::vroom(opt$input, delim = ",",
        col_types = c(
        Counts = "d",
        .default = "c"
    ))

# Transform to wide format
wide_data <- long_data %>%
    pivot_wider(
        names_from = Sample_ID,
        values_from = Counts,
        values_fill = 0
    )

# Export wide format
vroom_write(wide_data, opt$output, delim = "\t")
