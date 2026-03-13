#!/usr/bin/env Rscript
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

library(argparse)
library(tidyverse)
library(edgeR)
library(cowplot)
library(ggrepel)

parser <- ArgumentParser(description = "RNA-seq exploratory data analysis")
parser$add_argument("--counts", required = TRUE, help = "Count matrix path (TSV or RDS)")
parser$add_argument("--input-type", required = TRUE,
                    choices = c("featurecounts", "tximport"),
                    help = "Input format type")
parser$add_argument("--sample-tsv", required = TRUE, help = "Sample metadata TSV")
parser$add_argument("--design", required = TRUE, help = "Design formula string")
parser$add_argument("--out-pca", required = TRUE, help = "PCA plot PDF output")
parser$add_argument("--out-rds", required = TRUE, help = "EDA results RDS output")
args <- parser$parse_args()

samples <- read_tsv(args$sample_tsv, show_col_types = FALSE)

if (args$input_type == "featurecounts") {
  # Read merged featureCounts matrix (Geneid + lib columns)
  fc <- read.delim(args$counts, comment.char = "#")
  counts <- as.matrix(fc[, -1])
  rownames(counts) <- fc[[1]]
  y <- DGEList(counts = counts)
} else {
  txi <- readRDS(args$counts)
  counts <- txi$counts
  # Apply tximport length offset
  norm_mat <- txi$length
  norm_mat <- norm_mat / exp(rowMeans(log(norm_mat)))
  norm_counts <- counts / norm_mat
  eff_lib <- calcNormFactors(norm_counts) * colSums(norm_counts)
  norm_mat <- sweep(norm_mat, 2, eff_lib, "*")
  norm_mat <- log(norm_mat)
  y <- DGEList(counts)
  y <- scaleOffset(y, norm_mat)
}

# Build design matrix
formula <- as.formula(args$design)
design <- model.matrix(formula, data = samples)

# Filter lowly expressed genes
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize and compute logCPM
y <- calcNormFactors(y)
logCPM <- cpm(y, prior.count = 2, log = TRUE)

# PCA
pca <- prcomp(t(logCPM))
pve_pc1 <- round(100 * summary(pca)$importance[2, 1])
pve_pc2 <- round(100 * summary(pca)$importance[2, 2])

# Determine color factor (last term in design formula)
factor_vec <- all.vars(formula)
color_var <- tail(factor_vec, 1)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("library_id") %>%
  left_join(samples, by = "library_id")

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_var]], label = library_id)) +
  geom_point(size = 4) +
  geom_text_repel() +
  xlab(paste0("PC1 (", pve_pc1, "% variance)")) +
  ylab(paste0("PC2 (", pve_pc2, "% variance)")) +
  scale_color_discrete(name = color_var) +
  theme_cowplot() +
  theme(legend.position = "bottom")

# Save
save_plot(pca_plot, filename = args$out_pca)
saveRDS(list(design = design, logCPM = logCPM, pca = pca, y = y),
        file = args$out_rds)
