#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
cpm_tsv = args[1]
formula = args[2]
libraries_full_rds = args[3]
out_png = args[4]
out_svg = args[5]

factor_str = gsub("(~0 \\+)|\\s*\\*\\s*|\\s*\\+\\s*", " ", formula)
factor_str = trimws(factor_str)

factor_vec = strsplit(factor_str, " ")[[1]]
factor_vec = factor_vec[!factor_vec == "+"]

library(cowplot)
library(ggrepel)
library(tidyverse)

cpm = read_tsv(cpm_tsv)
libraries_full = readRDS(libraries_full_rds)

pca = prcomp(t(as.matrix(cpm[,-1])))

(pve_pc1=round(100*summary(pca)$importance[2,1]))

(pve_pc2=round(100*summary(pca)$importance[2,2]))

plot = as.data.frame(pca$x) %>%
  rownames_to_column(var = "library") %>%
  left_join(libraries_full, by = "library") %>%
  ggplot(., aes(x = PC1, y = PC2, color = get(factor_vec[[1]]), label = library)) +
  geom_point(size = 4) +
  geom_text_repel() +
  scale_color_discrete(name = factor_vec[[1]]) +
  xlab(paste("PC1, ", pve_pc1, "% variance explained", sep ="")) +
  ylab(paste("PC2, ", pve_pc2, "% variance explained", sep ="")) +
  coord_fixed(ratio = 1)

if (length(factor_vec) >= 2 && !is.null(factor_vec[[2]])) {
  plot = plot +
    aes(shape = get(factor_vec[[2]])) +
    scale_shape_discrete(name = factor_vec[[2]])
}

ggsave(filename = out_png, plot = plot, device = "png", width = 8, height = 6)
ggsave(filename = out_svg, plot = plot, device = "svg", width = 8, height = 6)
