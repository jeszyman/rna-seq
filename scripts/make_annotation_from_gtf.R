#!/usr/bin/env Rscript

#########################################
###   Make Annotate From A Gtf File   ###
#########################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
gtf_file = args[1]
bmart_dataset = args[2]
tsv = args[3]

#  "~/cards/ref/mm10.ensGene.gtf.gz"
#bmart_dataset = "mmusculus_gene_ensembl"
#tsv =

# Load required packages, data, and functions

library(biomaRt)
library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)

# Load the GTF file
gtf <- rtracklayer::import(gtf_file)

annotation = data.frame(ensembl_gene_id = gtf$gene_id) %>% distinct(ensembl_gene_id, .keep_all = TRUE)

mart = useMart("ensembl")
mart = useDataset(bmart_dataset, mart)

names = getBM(
  filters = "ensembl_gene_id",
  attributes=c("ensembl_gene_id",
               "entrezgene_id",
               "description",
               "external_gene_name",
               "gene_biotype"),
  values = annotation$ensembl_gene_id,
  mart = mart,
  uniqueRows = T)

names =
  names %>% group_by(ensembl_gene_id) %>% slice_head(n = 1)

write_tsv(names, file = tsv)
