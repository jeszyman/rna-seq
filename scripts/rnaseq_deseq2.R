#!/usr/bin/env Rscript
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

library(argparse)
library(tidyverse)
library(DESeq2)
library(tximport)

parser <- ArgumentParser(description = "Differential expression with DESeq2")
parser$add_argument("--counts", required = TRUE, help = "Count matrix (TSV or RDS)")
parser$add_argument("--input-type", required = TRUE, choices = c("featurecounts", "tximport"))
parser$add_argument("--sample-tsv", required = TRUE, help = "Sample metadata TSV")
parser$add_argument("--design", required = TRUE, help = "Design formula string")
parser$add_argument("--contrast", nargs = 3, required = TRUE,
                    help = "Contrast: factor level_test level_ref")
parser$add_argument("--out-tsv", required = TRUE, help = "Results TSV output")
parser$add_argument("--out-volcano", required = TRUE, help = "Volcano plot PDF")
parser$add_argument("--out-ma", required = TRUE, help = "MA plot PDF")
args <- parser$parse_args()

samples <- read_tsv(args$sample_tsv, show_col_types = FALSE)
formula <- as.formula(args$design)

if (args$input_type == "featurecounts") {
  cts <- read.delim(args$counts)
  count_mat <- as.matrix(cts[, -1])
  rownames(count_mat) <- cts[[1]]
  dds <- DESeqDataSetFromMatrix(countData = count_mat,
                                colData = samples,
                                design = formula)
} else {
  txi <- readRDS(args$counts)
  dds <- DESeqDataSetFromTximport(txi, colData = samples, design = formula)
}

dds <- DESeq(dds)
res <- results(dds, contrast = args$contrast)
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

write_tsv(res_df, args$out_tsv)

# Volcano plot
pdf(args$out_volcano)
plot(res_df$log2FoldChange, -log10(res_df$padj),
     pch = 20, cex = 0.5,
     xlab = "log2 Fold Change", ylab = "-log10(padj)",
     main = paste("DESeq2:", args$contrast[2], "vs", args$contrast[3]))
abline(h = -log10(0.05), col = "red", lty = 2)
dev.off()

# MA plot
pdf(args$out_ma)
plotMA(res, main = paste("DESeq2:", args$contrast[2], "vs", args$contrast[3]))
dev.off()
