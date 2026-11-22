#!/usr/bin/env Rscript

# Build a miRNA identifier table (name -> accession, precursor) from the miRBase or MirGeneDB
# GFF used for the run, plus the mature and hairpin FASTA files used for quantification.
# Released under the MIT license. See LICENSE https://github.com/nf-core/smrnaseq/blob/master/LICENSE for details.
#
# Usage: mirna_ids.r --gff <gff> --mature <fa> --hairpin <fa> --species <sps> --database <db> --output <tsv>
#
# Names are the identifiers the pipeline reports:
#   - miRBase GFF3: `Name` attribute (mirtop) == first word of the FASTA headers (edgeR)
#   - MirGeneDB GFF: `ID` attribute (mirtop) == FASTA headers (edgeR); star arms end in `*`
# Accessions come from the `Alias` attribute (MIMAT/MI). Precursors come from `Derives_from`
# (miRBase) or, when absent (MirGeneDB), from the hairpin whose coordinates contain the mature.
# FASTA names that are not in the GFF are added with empty annotation, so that every
# identifier in the edgeR and mirtop outputs is listed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2 != 0) stop("Arguments must be given as --key value pairs")
opt <- list()
for (i in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[i])]] <- args[i + 1]
for (key in c("gff", "mature", "hairpin", "species", "database", "output")) {
    if (is.null(opt[[key]])) stop(paste0("Missing argument --", key))
}

read_lines <- function(path) {
    # gzfile() also reads uncompressed files
    con <- gzfile(path, "r")
    on.exit(close(con))
    readLines(con, warn = FALSE)
}

# Same species filter as PARSE_FASTA_MIRNA (`.*<species>-.*`)
species <- opt$species
keep_species <- function(names) {
    if (!nzchar(species) || species == "unknown") return(rep(TRUE, length(names)))
    grepl(paste0(species, "-"), names, fixed = TRUE)
}

join_unique <- function(x) {
    x <- unique(x[!is.na(x) & nzchar(x)])
    paste(x[order(x, method = "radix")], collapse = ",")
}

# --- GFF ---------------------------------------------------------------------
lines <- read_lines(opt$gff)
lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
fields <- strsplit(lines, "\t", fixed = TRUE)
fields <- fields[lengths(fields) >= 9]
col <- function(i) vapply(fields, function(f) f[i], character(1))
gff <- data.frame(
    seqid      = col(1),
    feature    = col(3),
    start      = as.integer(col(4)),
    end        = as.integer(col(5)),
    strand     = col(7),
    attributes = col(9),
    stringsAsFactors = FALSE
)
hairpin_features <- c("miRNA_primary_transcript", "pre_miRNA")
gff <- gff[gff$feature %in% c("miRNA", hairpin_features), , drop = FALSE]

attr_value <- function(attributes, key) {
    prefix <- paste0(key, "=")
    vapply(strsplit(attributes, ";", fixed = TRUE), function(kv) {
        kv <- trimws(kv)
        hit <- kv[startsWith(kv, prefix)]
        if (length(hit)) substring(hit[1], nchar(prefix) + 1) else ""
    }, character(1))
}
gff$id        <- attr_value(gff$attributes, "ID")
name          <- attr_value(gff$attributes, "Name")
gff$name      <- ifelse(nzchar(name), name, gff$id)
gff$accession <- attr_value(gff$attributes, "Alias")
gff$derives   <- attr_value(gff$attributes, "Derives_from")
gff$type      <- ifelse(gff$feature == "miRNA", "mature", "hairpin")
gff <- gff[nzchar(gff$name) & keep_species(gff$name), , drop = FALSE]

hairpins <- gff[gff$type == "hairpin", , drop = FALSE]
matures  <- gff[gff$type == "mature", , drop = FALSE]

precursor_index <- function(i) {
    if (nzchar(matures$derives[i])) {
        which(hairpins$id %in% strsplit(matures$derives[i], ",", fixed = TRUE)[[1]])
    } else {
        which(hairpins$seqid == matures$seqid[i] & hairpins$strand == matures$strand[i] &
              hairpins$start <= matures$start[i] & hairpins$end >= matures$end[i])
    }
}
# Mature -> precursor links, one row per locus and precursor
links <- do.call(rbind, c(
    list(data.frame(name = character(), precursor = character(), precursor_accession = character(),
                    stringsAsFactors = FALSE)),
    lapply(seq_len(nrow(matures)), function(i) {
        idx <- precursor_index(i)
        if (!length(idx)) return(NULL)
        data.frame(name = matures$name[i], precursor = hairpins$name[idx],
                   precursor_accession = hairpins$accession[idx], stringsAsFactors = FALSE)
    })
))

# One row per name and type (a miRNA can have several genomic loci). Precursor names and
# accessions are listed in the same order, sorted by precursor name.
if (nrow(gff)) {
    keys <- unique(gff[, c("type", "name")])
    table <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
        rows <- gff[gff$type == keys$type[k] & gff$name == keys$name[k], , drop = FALSE]
        pre <- if (keys$type[k] == "mature") unique(links[links$name == keys$name[k], c("precursor", "precursor_accession")]) else links[0, c("precursor", "precursor_accession")]
        pre <- pre[order(pre$precursor, pre$precursor_accession, method = "radix"), , drop = FALSE]
        data.frame(
            name                = keys$name[k],
            accession           = join_unique(unlist(strsplit(rows$accession, ",", fixed = TRUE))),
            type                = keys$type[k],
            precursor           = paste(pre$precursor, collapse = ","),
            precursor_accession = paste(pre$precursor_accession, collapse = ","),
            stringsAsFactors    = FALSE
        )
    }))
} else {
    table <- data.frame(name = character(), accession = character(), type = character(),
                        precursor = character(), precursor_accession = character(),
                        stringsAsFactors = FALSE)
}

# --- FASTA names that are not in the GFF -------------------------------------
fasta_names <- function(path) {
    headers <- grep("^>", read_lines(path), value = TRUE)
    names <- sub("[[:space:]].*$", "", substring(headers, 2))
    unique(names[nzchar(names) & keep_species(names)])
}
for (type in c("mature", "hairpin")) {
    names <- fasta_names(opt[[type]])
    missing <- setdiff(names, table$name[table$type == type])
    if (length(missing)) {
        table <- rbind(table, data.frame(
            name = missing, accession = "", type = type, precursor = "", precursor_accession = "",
            stringsAsFactors = FALSE
        ))
    }
}

table$database <- rep(opt$database, nrow(table))
table <- table[order(table$type, table$name, method = "radix"), , drop = FALSE]
write.table(table, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
