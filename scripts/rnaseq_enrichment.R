#!/usr/bin/env Rscript
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

library(argparse)
library(tidyverse)
library(clusterProfiler)

parser <- ArgumentParser(description = "GO/KEGG/GSEA functional enrichment")
parser$add_argument("--de-results", required = TRUE, help = "DE results TSV (gene, log2FC, padj)")
parser$add_argument("--orgdb", required = TRUE, help = "OrgDb package name (e.g., org.Mm.eg.db)")
parser$add_argument("--out-go-tsv", required = TRUE)
parser$add_argument("--out-go-plot", required = TRUE)
parser$add_argument("--out-kegg-tsv", required = TRUE)
parser$add_argument("--out-kegg-plot", required = TRUE)
args <- parser$parse_args()

library(args$orgdb, character.only = TRUE)

de <- read_tsv(args$de_results, show_col_types = FALSE)

# Significant genes for ORA
sig_genes <- de %>% filter(padj < 0.05) %>% pull(gene)

# All genes ranked by log2FC for GSEA
gene_list <- de$log2FoldChange
names(gene_list) <- de$gene
gene_list <- sort(gene_list, decreasing = TRUE)

# GO over-representation analysis
go_res <- tryCatch(
  enrichGO(gene = sig_genes, OrgDb = args$orgdb,
           keyType = "ENSEMBL", ont = "BP",
           pAdjustMethod = "BH", pvalueCutoff = 0.05),
  error = function(e) NULL
)

if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
  write_tsv(as.data.frame(go_res), args$out_go_tsv)
  pdf(args$out_go_plot)
  print(dotplot(go_res, showCategory = 20, title = "GO Biological Process"))
  dev.off()
} else {
  write_tsv(tibble(message = "No significant GO terms found"), args$out_go_tsv)
  pdf(args$out_go_plot)
  plot.new(); text(0.5, 0.5, "No significant GO terms")
  dev.off()
}

# KEGG (requires ENTREZ IDs)
kegg_res <- tryCatch({
  entrez_map <- bitr(sig_genes, fromType = "ENSEMBL", toType = "ENTREZID",
                     OrgDb = args$orgdb)
  enrichKEGG(gene = entrez_map$ENTREZID,
             organism = ifelse(grepl("Mm", args$orgdb), "mmu", "hsa"),
             pvalueCutoff = 0.05)
}, error = function(e) NULL)

if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  write_tsv(as.data.frame(kegg_res), args$out_kegg_tsv)
  pdf(args$out_kegg_plot)
  print(barplot(kegg_res, showCategory = 20, title = "KEGG Pathways"))
  dev.off()
} else {
  write_tsv(tibble(message = "No significant KEGG pathways found"), args$out_kegg_tsv)
  pdf(args$out_kegg_plot)
  plot.new(); text(0.5, 0.5, "No significant KEGG pathways")
  dev.off()
}
