# For unit testing
#gsea_file_path = "/mnt/ris/jschwarz/Active/cardiac-radiobiology/analysis/wtrans"
#gsea_file_pattern = "_gsea"
#gsea_xlsx = "/tmp/gsea.xlsx"

# Command line arguments
args = commandArgs(trailingOnly = TRUE)
gsea_file_path = args[1]
gsea_file_pattern = args[2]
gsea_xlsx = args[3]

library(tidyverse)
library(openxlsx)

# Aggregate gsea results to single table
gsea_files = list.files(gsea_file_path, gsea_file_pattern, full.names = TRUE)

names(gsea_files) = list.files(gsea_file_path, gsea_file_pattern, full.names = FALSE)

gsea_dfs = lapply(gsea_files, read_tsv)

gsea = bind_rows(gsea_dfs, .id = "gsea") %>%
  mutate(species = gsub("_.*$", "", gsea)) %>%
  mutate(cohort = gsub(".*_(.*)_gsea.*", "\\1", gsea)) %>%
  mutate(pathway_set = gsub("_.*$","", pathway)) %>%
  select(species, cohort, pathway_set, everything()) %>%
  select(!gsea)
gsea

write.xlsx(gsea, gsea_xlsx)
