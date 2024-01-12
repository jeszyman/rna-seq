#!/usr/bin/env Rscript

# ---   Setup   --- #
# ----------------- #

## ---   Load Packages   --- ##
## ------------------------- ##

library(edgeR)
library(optparse)
library(tidyverse)

## ---   Load Inputs   --- ##
## ----------------------- ##

option_list <- list(
  make_option(c("--annotation_tsv"), type = "character", default = "~/cards/ref/mm10_wtrans_annotation.tsv"),
  make_option(c("--cohorts_str"), type = "character", default = "cohortir1d cohortsham"),
  make_option(c("--design_rds"), type = "character", default = "~/cards/analysis/rna/models/combat/hrt/design.rds"),
  make_option(c("--fit_rds"), type = "character", default = "~/cards/analysis/rna/models/combat/hrt/edger_fit.rds"),
  make_option(c("--res_tsv"), type = "character", default = "/tmp/test.tsv"))

opts <- parse_args(OptionParser(option_list = option_list))

list_of_options <- names(opts)
for (opt_name in list_of_options) {
  assign(opt_name, opts[[opt_name]], envir = .GlobalEnv)
}

design = readRDS(design_rds)
fit = readRDS(fit_rds)
annotation = read_tsv(annotation_tsv)


cohorts_vec = strsplit(cohorts_str, " ")[[1]]
contrast_string <- paste(cohorts_vec[[1]], "-", cohorts_vec[[2]])

contrast <- makeContrasts(eval(parse(text = contrast_string)), levels=design)

qlf = glmQLFTest(fit, contrast = contrast)

res =
  as.data.frame(topTags(qlf, n = Inf)) %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  as_tibble() %>%
  left_join(annotation, by = "ensembl_gene_id") %>%
  mutate(sign = sign(logFC)) %>%
  mutate(score = sign * -log10(PValue)) %>%
  mutate(rank = rank(-score, ties.method = "random"))

write_tsv(res, file = res_tsv)
