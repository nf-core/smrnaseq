#!/usr/bin/env Rscript

library(optparse)
library(tidyr)
library(vroom)

option_list <- list(
    make_option(c("--input"), type = "character", help = "Input TSV file", metavar = "character"),
    make_option(c("--output"), type = "character", help = "Output CSV file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Read CSV with vroom
data <- vroom::vroom(opt$input, delim = "\t", col_types = c(.default = "c"))

last_col <- names(data)[ncol(data)]

# Convert from wide to long format
long_data <- data %>%
    pivot_longer(
        cols = last_col,
        names_to = "Sample_ID",
        values_to = "Counts"
    )

vroom_write(long_data, opt$output, delim = ",")

