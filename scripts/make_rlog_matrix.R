#!/usr/bin/env Rscript

# Given an edgeR DGE object, make normalized counts by rlog for time-course comparison

# ---   Load Packages   --- #
# ------------------------- #

packages <- c("DESeq2", "optparse", "tidyverse")
sapply(packages, require, character.only = TRUE, quietly = TRUE)

# ---   Load Inputs   --- #
# ----------------------- #

option_list <- list(
  make_option(c("--dge_rds"), type = "character"),
  make_option(c("--formula_string"), type = "character"),
  make_option(c("--libraries_full_rds"), type = "character"),
  make_option(c("--out_mat_file"), type = "character")
  )

opts <- parse_args(OptionParser(option_list = option_list))

#opts$dge_rds = "~/cards/analysis/rna/models/combat/hrt/edger_dge.rds"
#opts$libraries_full_rds = "~/cards/data-model/lists/libraries_full.rds"
#opts$formula_string = "~ 0 + cohort"
#opts$out_mat_file = "/tmp/thetest.tsv"

inputs_list = function(opts){
  dge = readRDS(opts$dge_rds)
  formula_string = opts$formula_string
  libraries_full = readRDS(opts$libraries_full_rds)
  return(list(
    dge = dge,
    formula_string = formula_string,
    libraries_full = libraries_full)
  )
}

# ---   Main   --- #
# ---------------- #

main = function(opts){

  # Process inputs
  inputs = inputs_list(opts)

  # Make matrix
  mat = make_rlog_matrix(inputs$dge, inputs$libraries_full, inputs$formula_string)

  # Return a list of outputs
  return(mat)
}

# ---   Functions   --- #
# --------------------- #

make_rlog_matrix = function(dge, libraries_full, formula_string){
  mat <- matrix(as.integer(dge$counts), nrow = nrow(dge$counts), ncol = ncol(dge$counts), dimnames = dimnames(dge$counts))
  libs = data.frame(library = colnames(mat)) %>%
    left_join(libraries_full) %>% droplevels(.)
  deseq = DESeqDataSetFromMatrix(mat, libs, formula(formula_string))
  rlog_counts <- rlog(deseq, blind = TRUE)
  rlog_counts_matrix <- assay(rlog_counts)
  return(rlog_counts_matrix)
}

# ---   Run   --- #
# --------------- #

out_mat = main(opts)

saveRDS(out_mat, file = opts$out_mat_file)
