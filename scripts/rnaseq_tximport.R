#!/usr/bin/env Rscript
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

library(argparse)
library(tximport)

parser <- ArgumentParser(description = "Import Salmon quant.sf files via tximport")
parser$add_argument("--quant-files", nargs = "+", required = TRUE,
                    help = "Space-separated Salmon quant.sf file paths")
parser$add_argument("--gtf", required = TRUE,
                    help = "GTF annotation file for tx2gene mapping")
parser$add_argument("--out-rds", required = TRUE,
                    help = "Output RDS file path for tximport object")
args <- parser$parse_args()

# Build tx2gene from GTF
gtf <- rtracklayer::import(args$gtf)
tx2gene <- unique(data.frame(
  TXNAME = gtf$transcript_id[!is.na(gtf$transcript_id)],
  GENEID = gtf$gene_id[!is.na(gtf$transcript_id)]
))

# Name quant files by library ID (extract from path)
quant_files <- args$quant_files
names(quant_files) <- gsub(".*/(lib[0-9]+)\\..*", "\\1", quant_files)

# Run tximport
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)

# Save
saveRDS(txi, file = args$out_rds)
