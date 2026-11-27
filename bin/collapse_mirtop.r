#!/usr/bin/env Rscript

# Written by Lorena Pantano and released under the MIT license. See LICENSE https://github.com/nf-core/smrnaseq/blob/master/LICENSE for details.

# Collapse the joined mirtop isomiR table into a miRNA count table (mirna.tsv).
# See docs/output.md (mirtop section) for a description of the rules applied here.

library(data.table)
# Command line arguments
args = commandArgs(trailingOnly=TRUE)

input <- as.character(args[1:length(args)])

# check.names = FALSE keeps sample names such as "2-treated" unchanged (R would rename it to "X2.treated").
df = read.delim(input[1], sep = "\t", check.names = FALSE)
dt = as.data.table(df)
# Columns 1-12 are the mirtop annotation (UID, Read, miRNA, Variant, iso_*); sample counts follow.
sample_cols = colnames(dt)[13:ncol(dt)]
output = file.path(dirname(input[1]), "mirna.tsv")

# 1. Sum all rows of the same isomiR (UID) for the same miRNA. The same sequence can be reported
#    in several rows that differ only in precursor-dependent annotation columns (e.g. iso_*_nt);
#    these rows carry counts for different samples and must not be dropped (#174).
counts = dt[, lapply(.SD, sum, na.rm = TRUE), by = .(UID, miRNA), .SDcols = sample_cols]
mirna_order = unique(counts[["miRNA"]])

if (nrow(counts) == 0) {
    write.table(counts[, c("miRNA", sample_cols), with = FALSE], output, quote=FALSE, sep="\t", row.names=FALSE)
    quit(save = "no")
}

long = melt(counts, id.vars = c("UID", "miRNA"), measure.vars = sample_cols,
    variable.name = "sample", value.name = "count", variable.factor = FALSE)

# 2. A UID assigned to several miRNAs (ambiguous: the same reads are reported once per candidate
#    miRNA) is split between its candidates in proportion, per sample, to the candidates'
#    unambiguous counts (sum over the UIDs assigned only to that miRNA in that sample). If no
#    candidate has unambiguous counts in that sample, the UID is split equally. Each UID is
#    counted exactly once per sample, so sample totals are preserved. UIDs assigned to a single
#    miRNA are unaffected.
long[, n_candidates := .N, by = .(UID, sample)]
weights = long[n_candidates == 1, .(weight = sum(count)), by = .(miRNA, sample)]
long = merge(long, weights, by = c("miRNA", "sample"), all.x = TRUE, sort = FALSE)
long[is.na(weight), weight := 0]
long[, uid_count := max(count), by = .(UID, sample)]
long[, share := if (sum(weight) > 0) weight / sum(weight) else 1 / .N, by = .(UID, sample)]
long[, assigned := uid_count * share]

# 3. Sum per miRNA and round to integer counts with the largest remainder method, so that each
#    sample total stays exactly equal to the number of reads assigned to isomiRs. Ties are broken
#    by the order in which the miRNAs first appear in the input.
mirna = long[, .(value = sum(assigned)), by = .(miRNA, sample)]
mirna[, order := match(miRNA, mirna_order)]
setorder(mirna, sample, order)
mirna[, final := {
    target = round(sum(value))
    base = floor(value + 1e-9)
    remainder = round(value - base, 9)
    extra = target - sum(base)
    base + as.integer(frank(-remainder, ties.method = "first") <= extra)
}, by = sample]
mirna[, final := as.integer(final)]

wide = dcast(mirna, miRNA ~ sample, value.var = "final", fill = 0L)
wide = wide[match(mirna_order, wide[["miRNA"]])]
setcolorder(wide, c("miRNA", sample_cols))
write.table(wide, output, quote=FALSE, sep="\t", row.names=FALSE)
