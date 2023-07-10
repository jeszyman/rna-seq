args = commandArgs(trailingOnly = TRUE)
gtf = args[1]
build = args[2]

library(GenomicFeatures)

txdb = makeTxDbFromGFF(gtf)

saveDb(txdb, file = paste0("/tmp/",build,"_protein.txdb"))
