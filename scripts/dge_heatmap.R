#!/usr/bin/env Rscript

###################################################
###   Heatmap Of Differential Gene Expression   ###
###################################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
dge_tsv= args[1]
logcpm_tsv = args[2]
out_pdf = args[3]

# Load required packages, data, and functions
#dge_tsv = "~/cards/analysis/rna/contrasts/sod2_v_fvb_at_sham/sod2_v_fvb_at_sham.tsv"
#logcpm_tsv = "~/cards/analysis/rna/models/sod_kept/sod_kept_cpm.tsv"

dge_tsv = "~/cards/analysis/rna/models/sod_kept/edger_dge.rds"
formula_str = "~0 + cohort + sex"

library(ComplexHeatmap)
library(tidyverse)
library(ggsci)

dge = read_tsv(dge_tsv)


factor_str = gsub("(~0 \\+)|\\s*\\*\\s*|\\s*\\+\\s*", " ", formula_str)
factor_str = trimws(factor_str)
factor_vec = strsplit(factor_str, " ")[[1]]


genes = dge %>% filter(rank < 2001) %>% pull(ensembl_gene_id)
logcpm = read_tsv(logcpm_tsv) %>%
  filter(ensembl %in% genes)
mat = as.matrix(logcpm[,-1])

rowz =t(apply(mat, 1, scale))
colnames(rowz) = colnames(mat)

libraries_full = readRDS(libraries_full_rds)

libs = data.frame(library = colnames(mat)) %>%
  left_join(libraries_full, by = "library") %>%
  dplyr::select(all_of(factor_vec))
ha = HeatmapAnnotation(df = libs)

pdf(out_pdf)
ht = Heatmap(rowz, top_annotation = ha,
             show_column_names = TRUE,
             column_names_rot = 45,
             column_names_side = "bottom")

draw(ht)

dev.off()
