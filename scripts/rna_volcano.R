#!/usr/bin/env Rscript

################################
###   Rna-Seq Volcano Plot   ###
################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
dge_tsv = args[1]
out_pdf = args[2]

library(tidyverse)
library(ggrepel)

dge = read_tsv(dge_tsv)

table =
  dge %>%
  mutate(FDR_filtered = ifelse(abs(logFC) > 2, FDR, NA)) %>%
  # Create a column with the rank of each row, when ordered by FDR_filtered
  mutate(rank = rank(FDR_filtered, na.last = "keep")) %>%
  # Create the 'tolabel' column based on the rank column
  mutate(label = ifelse(rank <= 10 & !is.na(rank), external_gene_name,NA)) %>%
  mutate(sig = ifelse(FDR < 0.05 & abs(logFC) > 1, "Sig", "Not sig"))

plot = ggplot(table, aes(x = logFC, y = -log10(FDR), label = label)) +
  geom_point(aes(color = sig)) +
  scale_color_discrete(guide = "none") +
  geom_vline(xintercept = c(1,-1), linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_label_repel(box.padding = 1, show.legend = F) +
  theme_minimal() +
  ylab(expression(paste(-log["10"]*" ",italic("p")))) +
  xlab("Log-fold Change") +
  ggtitle("")

ggsave(plot, file = out_pdf)
