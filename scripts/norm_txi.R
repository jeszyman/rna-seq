#!/usr/bin/env Rscript

#######################
###   Human Edger   ###
#######################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
design_rds = args[1]
txi_rds = args[2]
dge_rds = args[3]
glm_rds = args[4]
logcpm_tsv = args[5]

# Load required packages, data, and functions
library(edgeR)
library(tidyverse)

design = readRDS(design_rds)
txi = readRDS(txi_rds)

# Make a DGE List
#  See https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
make_dge_list = function(txi, design){
  cts <- txi$counts
  normMat <- txi$length
  # Obtaining per-observation scaling factors for length, adjusted to avoid
  # changing the magnitude of the counts.
  normMat <- normMat/exp(rowMeans(log(normMat)))
  normCts <- cts/normMat
  # Computing effective library sizes from scaled counts, to account for
  # composition biases between samples.
  eff.lib <- calcNormFactors(normCts) * colSums(normCts)
  # Combining effective library sizes with the length factors, and calculating
  # offsets for a log-link GLM.
  normMat <- sweep(normMat, 2, eff.lib, "*")
  normMat <- log(normMat)
  # Creating a DGEList object for use in edgeR.
  y <- DGEList(cts)
  keep = filterByExpr(y, design)
  y = y[keep, ]
  return(y)
}

y = make_dge_list(txi, design)

logcpm = edgeR::cpm(y, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)

y = estimateDisp(y,design)
fit <- glmQLFit(y,design)

saveRDS(y, dge_rds)
saveRDS(fit, glm_rds)

logcpm %>% as.data.frame(.) %>% rownames_to_column(var = "ensembl") %>% as_tibble() %>% write_tsv(., file = logcpm_tsv)
