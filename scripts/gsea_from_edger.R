# Script to make a gsea table from edgeR results and msigdb pathways

# For unit testing
## results_file = "/mnt/ris/jschwarz/Active/cardiac-radiobiology/analysis/wtrans/human_ir2w-sham_edger_results.tsv"
## msig_str = "human_C2_CP:KEGG"
## enrichment_file = "/tmp/test.tsv"

# Command line arguments
args = commandArgs(trailingOnly = TRUE)

results_tsv = args[1]
msig_str = args[2]
enrichment_tsv = args[3]
enrichment_xlsx = args[4]

# Load necessary libraries
library(fgsea)
library(msigdbr)
library(tidyverse)
library(writexl)

results = read_tsv(results_tsv)
rank_tib = results %>%
  filter(F > 0) %>%
  mutate(SN = abs(logFC) / sqrt(1/F)) %>%
  arrange(desc(SN))
rank = rank_tib %>% pull(rank)
names(rank) = rank_tib$ensembl_gene_id
rank = rank[!duplicated(names(rank))]

msig_vect = as.character(unlist(strsplit(msig_str, "_")))

make_pathway_set = function(msig){
  # Pull in pathway set from MSigDb
  if (length(msig) == 2) {
    tib = msigdbr(msig[1], msig[2])
  } else {
    tib = msigdbr(msig[1], msig[2], msig[3])
  }
  pathways = split(as.character(tib$ensembl_gene), tib$gs_name)
  return(pathways)
}

pathway_set = make_pathway_set(msig_vect)

# str(head(pathway_set))

run_fgsea = function(pathways, stats){
  gsea = fgseaMultilevel(pathways = pathways,
                         stats = stats,
                         scoreType = "pos")
  gsea = as_tibble(gsea) %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ",")) %>% arrange(padj, pval, -ES)
  return(gsea)
}

enrichment = run_fgsea(pathway_set, rank)
#enrichment

write_tsv(enrichment, enrichment_tsv)
write_xlsx(enrichment, enrichment_xlsx)
