#!/usr/bin/env bash
in_gtf="${1}"
out_gtf="${2}"

zcat $in_gtf | awk '$0 ~ /gene_biotype "protein_coding"|gene_biotype "lncRNA"/' | gzip > $out_gtf
