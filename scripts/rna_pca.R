
## - Rscript

#!/usr/bin/env Rscript
dds_human_rds = "~/cards/analysis/rna/de/human_icell_wtrans_hg38_singf_lrt_deseq.rds"
dds_mouse_wtrans_rds = "~/cards/analysis/rna/de/mouse_hrt_wtrans_mm10_day_deseq.rds"
dds_nuc_rds = "~/cards/analysis/rna/de/mouse_nuc_bulk_mm10_singf_lrt_deseq.rds"
libraries_full_rds = "~/cards/data-model/libraries_full.rds"

########################
###   Rna-Seq Pcas   ###
########################

# Load required packages, data, and functions
library(DESeq2)
library(ggrepel)
library(tidyverse)

dds_human = readRDS(dds_human_rds)
dds_mouse_wtrans = readRDS(dds_mouse_wtrans_rds)
dds_nuc = readRDS(dds_nuc_rds)
libraries_full = readRDS(libraries_full_rds)

pca_list = function(dds, full_libs){
  vsd = vst(dds, blind = FALSE)
  mat = assay(vsd)
  pca = prcomp(t(mat))
  pve_pc1 = round(100*summary(pca)$importance[2,1])
  pve_pc2 = round(100*summary(pca)$importance[2,2])
  pca_plot = as.data.frame(pca$x) %>%
    rownames_to_column(var = "library") %>%
    left_join(full_libs, by = "library") %>% as_tibble()
  pca_list = list(pca_plot, pve_pc1, pve_pc2)
  names(pca_list) = c("tib", "pve1", "pve2")
  return(pca_list)
  }

#########1#########2#########3#########4#########5#########6#########7#########8

pca = pca_list(dds_nuc, libraries_full)
pca_plot =
  ggplot(pca$tib, aes(x = PC1, y = PC2, shape = gy, color = gy)) +
  geom_point(size = 4) +
  xlab(paste("PC1, ", pca$pve1, "% variance explained", sep ="")) +
  ylab(paste("PC2, ", pca$pve2, "% variance explained", sep =""))

ggsave(pca_plot, file = "~/cards/results/rna/nuc/nuc_rna_pca.pdf", width = 8, height = 6)

#########1#########2#########3#########4#########5#########6#########7#########8

pca = pca_list(dds_mouse_wtrans, libraries_full)
pca_plot =
  ggplot(pca$tib, aes(x = PC1, y = PC2, shape = run, color = post_ir_d)) +
  geom_point(size = 4) +
  xlab(paste("PC1, ", pca$pve1, "% variance explained", sep ="")) +
  ylab(paste("PC2, ", pca$pve2, "% variance explained", sep =""))

ggsave(pca_plot, file = "~/cards/results/rna/wtrans/mouse_wtrans_pca.pdf", width = 8, height = 6)

#########1#########2#########3#########4#########5#########6#########7#########8

dds = dds_mouse_wtrans
dds = dds[, dds$run=="Rentschler_s4630_MGI0042"]
dds$post_ir_d = droplevels(dds$post_ir_d)
pca = pca_list(dds, libraries_full)

pca_plot =
  ggplot(pca$tib, aes(x = PC1, y = PC2, shape = gy, color = post_ir_d)) +
  geom_point(size = 4) +
  xlab(paste("PC1, ", pca$pve1, "% variance explained", sep ="")) +
  ylab(paste("PC2, ", pca$pve2, "% variance explained", sep =""))

ggsave(pca_plot, file = "~/cards/results/rna/wtrans/mouse_wtrans_s4630_pca.pdf", width = 8, height = 6)

#########1#########2#########3#########4#########5#########6#########7#########8

dds = dds_mouse_wtrans
dds = dds[, dds$run=="Rentschler_s4730_MGI0070"]
dds$post_ir_d = droplevels(dds$post_ir_d)
pca = pca_list(dds, libraries_full)

pca_plot =
  ggplot(pca$tib, aes(x = PC1, y = PC2, shape = gy, color = post_ir_d)) +
  geom_point(size = 4) +
  xlab(paste("PC1, ", pca$pve1, "% variance explained", sep ="")) +
  ylab(paste("PC2, ", pca$pve2, "% variance explained", sep =""))

ggsave(pca_plot, file = "~/cards/results/rna/wtrans/mouse_wtrans_s4730_pca.pdf", width = 8, height = 6)

#########1#########2#########3#########4#########5#########6#########7#########8

pca = pca_list(dds_human, libraries_full)
pca_plot =
  ggplot(pca$tib, aes(x = PC1, y = PC2, shape = gy, color = post_ir_d)) +
  geom_point(size = 4) +
  xlab(paste("PC1, ", pca$pve1, "% variance explained", sep ="")) +
  ylab(paste("PC2, ", pca$pve2, "% variance explained", sep =""))

ggsave(pca_plot, file = "~/cards/results/rna/wtrans/human_wtrans_pca.pdf", width = 8, height = 6)
