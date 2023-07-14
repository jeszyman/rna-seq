#!/usr/bin/env Rscript

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
cohorts_str = args[1]
design_rds = args[2]
fit_rds = args[3]
res_tsv = args[4]

# Load required packages, data, and functions
library(edgeR)
library(tidyverse)

design = readRDS(design_rds)
dge = readRDS(dge_rds)
fit = readRDS(fit_rds)

cohorts_vec = strsplit(cohorts_str, " ")[[1]]
contrast_string <- paste(cohorts_vec[[1]], "-", cohorts_vec[[2]])

contrast <- makeContrasts(eval(parse(text = contrast_string)), levels=design)
qlf = glmQLFTest(fit, contrast = contrast)
res = as.data.frame(topTags(qlf, n = Inf)) %>% rownames_to_column(var = "ensembl") %>% as_tibble()

write_tsv(res, file = res_tsv)
