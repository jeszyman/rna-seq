#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
salmon_str = args[1]
txdb = args[2]
out_txi = args[3]

# Load libraries
library(tximport)
library(AnnotationDbi)

txdb = loadDb(txdb)

# Make salmon file list
salmon_vect = unlist(strsplit(salmon_str, " "))
names(salmon_vect) = substr(gsub("^.*lib", "lib", salmon_vect), 1, 6)

# Make gene annotation
k = keys(txdb, keytype = "TXNAME")
tx2gene = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Make txi object
txi = tximport(salmon_vect, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T, ignoreAfterBar = T)

# Save txi object
saveRDS(txi, file = out_txi)
