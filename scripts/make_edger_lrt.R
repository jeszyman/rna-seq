#!/usr/bin/env Rscript
#!/usr/bin/env Rscript

########################################################
###   Make Edger Contrast By Likelihood Ratio Test   ###
########################################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
design_rds = args[1]
fit_rds = args[2]
res_tsv = args[3]

# Load required packages, data, and functions
design = readRDS(design_rds)
fit = readRDS(fit_rds)

lrt = glmLRT(fit, coef = 2)


res =
  as.data.frame(topTags(lrt, n = Inf)) %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  as_tibble() %>%
  left_join(annotation, by = "ensembl_gene_id") %>%
  mutate(sign = sign(logFC)) %>%
  mutate(score = sign * -log10(PValue)) %>%
  mutate(rank = rank(-score, ties.method = "random"))

write_tsv(res, file = res_tsv)
