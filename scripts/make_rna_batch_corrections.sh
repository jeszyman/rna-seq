#!/usr/bin/env Rscript

# ---   Setup   --- #
# ----------------- #

## ---   Load Packages   --- ##
## ------------------------- ##

library(optparse)
library(cowplot)
library(edgeR)
library(gridExtra)
library(ggrepel)
library(patchwork)
library(sva)
library(tidyverse)

## ---   Load Inputs   --- ##
## ----------------------- ##

option_list <- list(
  make_option(c("--batch_var"), type = "character", default = "run"),
  make_option(c("--covars"), type = "character", default = "cohort"),
  make_option(c("--design_rds"), type = "character", default = "~/cards/analysis/rna/models/hrt/design.rds"),
  make_option(c("--dge_rds"), type = "character", default = "/tmp/dge.rds"),
  make_option(c("--glm_rds"), type = "character", default = "/tmp/glm.rds"),
  make_option(c("--libraries_full_rds"), type = "character", default = "~/cards/data-model/lists/libraries_full.rds"),
  make_option(c("--adjusted_rds"), type = "character", default = "/tmp/adjusted.rds"),
  make_option(c("--pdf"), type = "character", default = "/tmp/adjusted.pdf"),
  make_option(c("--txi_rds"), type = "character", default = "~/cards/analysis/rna/models/hrt/txi.rds")
)

opts <- parse_args(OptionParser(option_list = option_list))

list_of_options <- names(opts)
for (opt_name in list_of_options) {
  assign(opt_name, opts[[opt_name]], envir = .GlobalEnv)
}

design <- readRDS(design_rds)
libraries_full <- readRDS(libraries_full_rds)
txi <- readRDS(txi_rds)
counts <- txi$counts

batch <- data.frame(library = colnames(counts)) %>%
  left_join(libraries_full, by = "library") %>%
  pull(!!sym(batch_var)) %>%
  as.factor() %>%
  as.integer()

covars_vect <- unlist(strsplit(covars, " "))

# ---   Perform Combatseq   --- #
# ----------------------------- #

# Check if there is only one covariate
if (length(covars_vect) == 1) {
  covar <- data.frame(library = colnames(counts)) %>%
    left_join(libraries_full, by = "library") %>%
    pull(!!sym(covars_vect[1])) %>%
    as.factor() %>%
    as.integer()
  adjusted <- ComBat_seq(counts,
                         batch = batch,
                         group = covar)
} else { # If there are multiple covariates
  covar_list <- list() # Initialize an empty list to store covariate vectors
  for (covar_name in covars_vect) {
    covar <- data.frame(library = colnames(counts)) %>%
      left_join(libraries_full, by = "library") %>%
      pull(!!sym(covar_name)) %>%
      as.factor() %>%
      as.integer()
    covar_list[[covar_name]] <- covar
  }
  # Combine all covariate vectors into a matrix
  covar_mat <- do.call(cbind, covar_list)
  adjusted <- ComBat_seq(counts,
                         batch = batch,
                         group = NULL,
                         covar_mod = covar_mat)
}

# ---   Normalize Adjusted Counts Across Exp Design   --- #
# ------------------------------------------------------- #

make_dge_list <- function( eff_gene_len, counts_mat, design ) {
  normMat <- eff_gene_len
  # Obtaining per-observation scaling factors for length, adjusted to avoid
  # changing the magnitude of the counts.
  normMat <- normMat/exp(rowMeans(log(normMat)))
  normCts <- counts_mat/normMat
  # Computing effective library sizes from scaled counts, to account for
  # composition biases between samples.
  eff.lib <- calcNormFactors(normCts) * colSums(normCts)
  # Combining effective library sizes with the length factors, and calculating
  # offsets for a log-link GLM.
  normMat <- sweep(normMat, 2, eff.lib, "*")
  normMat <- log(normMat)
  # Creating a DGEList object for use in edgeR.
  y <- DGEList(counts_mat)
  y <- scaleOffset(y, normMat)
  keep <- filterByExpr(y, design)
  y <- y[keep, ]
  return(y)
}

y <- make_dge_list(txi$length, adjusted, design)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)

logcpm_adjusted <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)

logcpm <- edgeR::cpm(make_dge_list(txi$length, txi$counts, design), normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)

saveRDS(y, dge_rds)
saveRDS(fit, file = glm_rds)

# ---   Make PCAs   --- #
# --------------------- #

find_pca_extents <- function(logcpm1, logcpm2) {
  pca1 <- prcomp(t(as.matrix(logcpm1[,-1])))
  pca2 <- prcomp(t(as.matrix(logcpm2[,-1])))

  combined <- rbind(as.data.frame(pca1$x), as.data.frame(pca2$x))
  xlims <- range(combined$PC1, na.rm = TRUE)
  ylims <- range(combined$PC2, na.rm = TRUE)

  return(list(xlims = xlims, ylims = ylims))
}

extents <- find_pca_extents(logcpm, logcpm_adjusted)

make_pca <- function(logcpm, libraries_full, covars_vect, xlims, ylims) {
  pca <- prcomp(t(as.matrix(logcpm[,-1])))
  pve_pc1 <- round(100 * summary(pca)$importance[2,1])
  pve_pc2 <- round(100 * summary(pca)$importance[2,2])
  pca_plot <- as.data.frame(pca$x) %>%
    rownames_to_column(var = "library") %>%
    left_join(libraries_full, by = "library") %>%
    ggplot(aes(x = PC1, y = PC2, label = library, color = !!sym(covars_vect[[1]]))) +
    geom_point(size = 4) +
    geom_text_repel() +
    xlab(paste("PC1, ", pve_pc1, "% variance explained", sep ="")) +
    ylab(paste("PC2, ", pve_pc2, "% variance explained", sep ="")) +
    coord_fixed(ratio = 1) +
    xlim(xlims) +
    ylim(ylims)

  return(pca_plot)
}

pca <- make_pca(logcpm, libraries_full, covars_vect, extents$xlims, extents$ylims)
adjusted_pca <- make_pca(logcpm_adjusted, libraries_full, covars_vect, extents$xlims, extents$ylims)

combined_plot <- (pca + ggtitle("Unadjusted")) + (adjusted_pca + ggtitle("Adjusted")) +
                 plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

pdf("/tmp/test.pdf", width = 14, height =7, title = "combined plot")
combined_plot
dev.off()
