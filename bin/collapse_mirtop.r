#!/usr/bin/env Rscript

# Written by Lorena Pantano and released under the MIT license. See LICENSE https://github.com/nf-core/smrnaseq/blob/master/LICENSE for details.

library(data.table)
# Command line arguments
args = commandArgs(trailingOnly=TRUE)

input <- as.character(args[1:length(args)])

# check.names = FALSE keeps sample names such as "2-treated" unchanged (R would rename it to "X2.treated").
df = read.delim(input[1], sep = "\t", check.names = FALSE)
dt = as.data.table(df)
# Columns 1-12 are the mirtop annotation (UID, Read, miRNA, Variant, iso_*); sample counts follow.
sample_cols = colnames(dt)[13:ncol(dt)]

# 1. Sum all rows of the same isomiR (UID) for the same miRNA. The same sequence can be reported
#    in several rows that differ only in precursor-dependent annotation columns (e.g. iso_*_nt);
#    these rows carry counts for different samples and must not be dropped (#174).
counts = dt[, lapply(.SD, sum, na.rm = TRUE), by = .(UID, miRNA), .SDcols = sample_cols]
# 2. A UID assigned to more than one miRNA (multi-mapping, e.g. identical mature sequences) is
#    counted only once, for the first miRNA it appears with, to avoid counting the same reads twice.
counts = counts[!duplicated(counts[["UID"]])]
# 3. Sum isomiR counts per miRNA.
mirna = counts[, lapply(.SD, sum), by = miRNA, .SDcols = sample_cols]
write.table(mirna, file.path(dirname(input[1]), "mirna.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
