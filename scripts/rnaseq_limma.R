#!/usr/bin/env Rscript
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

library(argparse)
library(tidyverse)
library(edgeR)
library(limma)

parser <- ArgumentParser(description = "Differential expression with limma-voom")
parser$add_argument("--counts", required = TRUE)
parser$add_argument("--input-type", required = TRUE, choices = c("featurecounts", "tximport"))
parser$add_argument("--sample-tsv", required = TRUE)
parser$add_argument("--design", required = TRUE)
parser$add_argument("--contrast", nargs = 3, required = TRUE)
parser$add_argument("--out-tsv", required = TRUE)
parser$add_argument("--out-volcano", required = TRUE)
parser$add_argument("--out-ma", required = TRUE)
args <- parser$parse_args()

samples <- read_tsv(args$sample_tsv, show_col_types = FALSE)
formula <- as.formula(args$design)

if (args$input_type == "featurecounts") {
  cts <- read.delim(args$counts)
  count_mat <- as.matrix(cts[, -1])
  rownames(count_mat) <- cts[[1]]
  y <- DGEList(counts = count_mat)
} else {
  txi <- readRDS(args$counts)
  y <- DGEList(counts = txi$counts)
  norm_mat <- txi$length
  norm_mat <- norm_mat / exp(rowMeans(log(norm_mat)))
  norm_counts <- txi$counts / norm_mat
  eff_lib <- calcNormFactors(norm_counts) * colSums(norm_counts)
  norm_mat <- sweep(norm_mat, 2, eff_lib, "*")
  y <- scaleOffset(y, log(norm_mat))
}

design <- model.matrix(formula, data = samples)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

v <- voom(y, design)
fit <- lmFit(v, design)

# Build contrast
contrast_factor <- args$contrast[1]
contrast_test <- args$contrast[2]
contrast_ref <- args$contrast[3]
coef_test <- paste0(contrast_factor, contrast_test)
coef_ref <- paste0(contrast_factor, contrast_ref)

# If reference level is the intercept, coef_ref won't exist — use coef_test directly
if (coef_test %in% colnames(design)) {
  contrast_vec <- makeContrasts(contrasts = coef_test, levels = design)
  fit2 <- contrasts.fit(fit, contrast_vec)
} else {
  fit2 <- fit
}

fit2 <- eBayes(fit2)
res <- topTable(fit2, number = Inf, sort.by = "P")

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  rename(log2FoldChange = logFC, padj = adj.P.Val, pvalue = P.Value) %>%
  arrange(padj)

write_tsv(res_df, args$out_tsv)

# Volcano
pdf(args$out_volcano)
plot(res_df$log2FoldChange, -log10(res_df$padj),
     pch = 20, cex = 0.5,
     xlab = "log2 Fold Change", ylab = "-log10(adj.P.Val)",
     main = paste("limma:", contrast_test, "vs", contrast_ref))
abline(h = -log10(0.05), col = "red", lty = 2)
dev.off()

# MA
pdf(args$out_ma)
plotMD(fit2, main = paste("limma:", contrast_test, "vs", contrast_ref))
abline(h = 0, col = "blue")
dev.off()
