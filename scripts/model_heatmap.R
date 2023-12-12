#!/usr/bin/env Rscript

# ---   Setup   --- #
# ----------------- #

## ---   Load Packages   --- ##
## ------------------------- ##

library(optparse)
library(ComplexHeatmap)
library(tidyverse)
library(ggsci)

## ---   Load Inputs   --- ##
## ----------------------- ##

option_list <- list(
  make_option(c("--dge_rds"), type = "character", default = "~/cards/analysis/rna/models/combat/hrt/edger_dge.rds"),
  make_option(c("--fct_str"), type = "character", default = "cohort sex run"),
  make_option(c("--libs_rds"), type = "character", default = "~/cards/data-model/lists/libraries_full.rds"),
  make_option(c("--nrow"), type = "character", default = "100"),
  make_option(c("--pdf_out"), type = "character", default = "/tmp/test.pdf")
)

opts <- parse_args(OptionParser(option_list = option_list))

list_of_options <- names(opts)
for (opt_name in list_of_options) {
  assign(opt_name, opts[[opt_name]], envir = .GlobalEnv)
}

libs <- readRDS(libs_rds)
dge <- readRDS(dge_rds)
factor_vec = strsplit(fct_str, " ")[[1]]

# ---   Heatmap   --- #
# ------------------- #

# Log CPM values
logcpm <- edgeR::cpm(dge,
                     normalized.lib.sizes = TRUE,
                     log = TRUE,
                     prior.count = 2)
gene_variances <- apply(logcpm, 1, var)
top_genes_var <- order(gene_variances, decreasing = TRUE)[1:nrow]
selected_genes <- logcpm[top_genes_var, ]
logcpm=selected_genes

# Z score
mat <- as.matrix(logcpm[,-1])
rowz <-t(apply(mat, 1, scale))
colnames(rowz) <- colnames(mat)

# Make column annotations
col_data <- data.frame(library = colnames(mat)) %>%
  left_join(libs, by = "library") %>%
  dplyr::select(all_of(factor_vec))
col_anno <- HeatmapAnnotation(df = col_data)

heat <- Heatmap(rowz,
                name = "Z-score",
                show_row_names = FALSE,
                top_annotation = col_anno,
                column_names_rot = 45,
                column_names_side = "bottom")

pdf(file = pdf_out, width = 10, height = 5)
draw(heat)
dev.off()
