#!/usr/bin/env Rscript
library(data.table)
# Command line arguments
args = commandArgs(trailingOnly=TRUE)

input <- as.character(args[1:length(args)])

df = read.delim(input[1], sep = "\t")
counts = as.data.table(df[!duplicated(df[["UID"]]),c(3, 13:ncol(df))])
mirna = counts[, lapply(.SD, sum), by = miRNA]
write.table(mirna, file.path(dirname(input[1]), "mirna.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
