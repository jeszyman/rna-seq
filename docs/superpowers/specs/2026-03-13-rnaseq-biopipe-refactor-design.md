# RNA-seq Biopipe Module Refactor — Design Spec

## Overview

Full rewrite of the rna-seq repository to match established biopipe module patterns (emseq, frag). Covers FASTQ through differential expression and functional enrichment, with dual alignment paths (STAR + Salmon), three DE tools (DESeq2, edgeR, limma), and test data from two organisms (mouse + human).

## Org File Structure

```
* RNA-seq Bioinformatics Processes                              :biopipe:
** README
   - Changelog (versioned per wf/rnaseq/vX.Y.Z convention)
** Repository setup and administration
   *** Org update                    (→ tools/shell/org_update.sh)
** Conda environments
   *** Alignment                     (→ config/rnaseq-align-conda-env.yaml)
   *** BAM QC                        (→ config/rnaseq-qc-conda-env.yaml)
   *** R                             (→ config/rnaseq-r-env.yaml)
** RNA-seq processing
   *** Preamble
   *** FASTQ QC
       **** fastp                    (rule rnaseq_fastp)
       **** FastQC                   (rule rnaseq_fastqc)
   *** Alignment: STAR
       **** STAR index               (rule rnaseq_star_index)
       **** STAR align               (rule rnaseq_star_align)
       **** featureCounts            (rule rnaseq_featurecounts)
   *** Alignment: Salmon
       **** Salmon index             (rule rnaseq_salmon_index)
       **** Salmon quant             (rule rnaseq_salmon_quant)
       **** tximport                 (rule rnaseq_tximport)
   *** BAM QC (STAR path only)
       **** RSeQC                    (rule rnaseq_rseqc)
       **** Qualimap                 (rule rnaseq_qualimap)
   *** QC aggregation
       **** MultiQC                  (rule rnaseq_multiqc)
   *** Exploratory data analysis
       **** EDA                      (rule rnaseq_eda)
   *** Differential expression
       **** DESeq2                   (rule rnaseq_deseq2)
       **** edgeR                    (rule rnaseq_edger)
       **** limma-voom               (rule rnaseq_limma)
   *** Functional enrichment
       **** Enrichment               (rule rnaseq_enrichment)
** Testing
   *** Test data provenance          (→ tools/get_test_data.sh)
   *** Test configuration            (→ config/test.yaml)
   *** Test wrapper                  (→ workflows/test.smk)
   *** GitHub CI/CD
       **** Generate data workflow   (→ .github/workflows/generate-data.yaml)
       **** Pipeline test workflow   (→ .github/workflows/test.yaml)
** [Existing background/reference text — retained, reorganized]
** General reference
** Development
```

## Directory Layout

```
rna-seq/
├── .github/workflows/
│   ├── generate-data.yaml
│   └── test.yaml
├── config/
│   ├── rnaseq-align-conda-env.yaml
│   ├── rnaseq-qc-conda-env.yaml
│   ├── rnaseq-r-env.yaml
│   └── test.yaml
├── data/
│   └── test-samples.tsv
├── scripts/
│   ├── rnaseq_tximport.R
│   ├── rnaseq_eda.R
│   ├── rnaseq_deseq2.R
│   ├── rnaseq_edger.R
│   ├── rnaseq_limma.R
│   └── rnaseq_enrichment.R
├── tools/
│   ├── get_test_data.sh
│   └── shell/org_update.sh
├── tests/
│   └── full/
│       └── inputs/          # Committed: FASTQs, FASTA subsets, GTF subsets
├── workflows/
│   ├── rnaseq.smk           # Modular — no rule all, no preamble
│   └── test.smk             # Preamble, config, rule all, includes rnaseq.smk
├── resources/
├── rna-seq.org
├── CLAUDE.md
└── README.md
```

## Conda Environments

### rnaseq-align-conda-env.yaml
Alignment + FASTQ processing section:
- STAR, Salmon, fastp, FastQC, samtools, subread (featureCounts), MultiQC

### rnaseq-qc-conda-env.yaml
BAM QC section (STAR path only):
- RSeQC, Qualimap

### rnaseq-r-env.yaml
All R/Bioconductor work across pipeline:
- DESeq2, edgeR, limma, tximport, clusterProfiler
- tidyverse, ggrepel, cowplot, pheatmap
- org.Mm.eg.db, org.Hs.eg.db (annotation DBs for enrichment)

## Snakemake Architecture

### Modular Design

`workflows/rnaseq.smk` contains all `rnaseq_*` rules. No `rule all`, no preamble, no config parsing. Designed to be included by a wrapper.

`workflows/test.smk` is the test wrapper:
- Loads `config/test.yaml`
- Parses config into variables (directory paths, conda env paths, etc.)
- Parses `de_map` from config
- Defines `rule all` with all target outputs
- Defines input symlink rules
- `include: "rnaseq.smk"`

### Configuration (test.yaml)

```yaml
main-data-dir: tests/full
sample-tsv: data/test-samples.tsv

refs:
  mm10_chr19:
    fasta: tests/full/inputs/mm10_chr19.fa.gz
    gtf: tests/full/inputs/mm10_chr19.gtf.gz
    star_sa_bases: 11    # --genomeSAindexNbases for small genome
  hg38_chr22:
    fasta: tests/full/inputs/hg38_chr22.fa.gz
    gtf: tests/full/inputs/hg38_chr22.gtf.gz
    star_sa_bases: 11

envs:
  align: config/rnaseq-align-conda-env.yaml
  qc: config/rnaseq-qc-conda-env.yaml
  r: config/rnaseq-r-env.yaml

repos:
  rnaseq: .

de-map:
  mouse_test:
    libs: ["lib001", "lib002", "lib003", "lib004"]
    ref_name: mm10_chr19
    align_method: ["star", "salmon"]
    design: "~ group"
    contrast: ["group", "exp", "ctrl"]
    tools: ["deseq2", "edgeR", "limma"]
    orgdb: "org.Mm.eg.db"
  human_test:
    libs: ["lib005", "lib006", "lib007", "lib008"]
    ref_name: hg38_chr22
    align_method: ["star", "salmon"]
    design: "~ group"
    contrast: ["group", "tumor", "normal"]
    tools: ["deseq2", "edgeR", "limma"]
    orgdb: "org.Hs.eg.db"
```

### Wildcards

- `{library_id}` — per-sample (fastp, fastqc, alignment, quantification)
- `{ref_name}` — reference genome (mm10_chr19, hg38_chr22)
- `{align_method}` — "star" or "salmon" (file name designator, like emseq)
- `{experiment}` — DE experiment key from de-map
- `{processing}` — "raw" or "trimmed" (FastQC)

### File Naming Convention

Method designator in filenames (like emseq's `{align_method}`):
```
{D_RNASEQ}/fastqs/{library_id}.{processing}.R{read}.fastq.gz
{D_RNASEQ}/bams/{library_id}.{ref_name}.star.bam
{D_RNASEQ}/quant/{library_id}.{ref_name}.salmon/quant.sf
{D_RNASEQ}/counts/{experiment}.{align_method}.txi.rds
{D_RNASEQ}/counts/{experiment}.{align_method}.featurecounts.tsv
{D_RNASEQ}/de/{experiment}.{align_method}.deseq2_results.tsv
{D_RNASEQ}/de/{experiment}.{align_method}.edger_results.tsv
{D_RNASEQ}/de/{experiment}.{align_method}.limma_results.tsv
{D_RNASEQ}/enrichment/{experiment}.{align_method}.{de_tool}.go.tsv
{D_RNASEQ}/qc/multiqc.html
{D_RNASEQ}/eda/{experiment}.{align_method}.pca.pdf
{D_RNASEQ}/eda/{experiment}.{align_method}.eda.rds
```

### de_map Dispatch Pattern

Same as emseq's meth_map:
```python
de_map = config["de-map"]

# In rules — resolve per-experiment params:
rule rnaseq_deseq2:
    input:
        counts = lambda wc: get_count_matrix(wc, de_map),
        sample_tsv = SAMPLE_TSV,
    params:
        design = lambda wc: de_map[wc.experiment]["design"],
        contrast = lambda wc: " ".join(de_map[wc.experiment]["contrast"]),
        orgdb = lambda wc: de_map[wc.experiment]["orgdb"],
    ...

# In rule all — expand across experiments and methods:
expand("{D_RNASEQ}/de/{experiment}.{align_method}.deseq2_results.tsv",
       experiment=[e for e in de_map if "deseq2" in de_map[e]["tools"]],
       align_method=lambda: de_map[e]["align_method"]),
```

### Rules Summary (16 rules in rnaseq.smk)

| # | Rule | Conda | Input | Output |
|---|------|-------|-------|--------|
| 1 | `rnaseq_fastp` | align | raw FASTQ R1+R2 | trimmed FASTQ R1+R2 + JSON |
| 2 | `rnaseq_fastqc` | align | raw/trimmed FASTQ | HTML + zip |
| 3 | `rnaseq_star_index` | align | FASTA + GTF | STAR genome dir |
| 4 | `rnaseq_star_align` | align | trimmed FASTQ + index | sorted BAM |
| 5 | `rnaseq_featurecounts` | align | BAM + GTF | count matrix TSV |
| 6 | `rnaseq_salmon_index` | align | FASTA + GTF (gentrome) | Salmon index dir |
| 7 | `rnaseq_salmon_quant` | align | trimmed FASTQ + index | quant.sf |
| 8 | `rnaseq_tximport` | r | quant.sf files | txi.rds |
| 9 | `rnaseq_rseqc` | qc | BAM | gene body coverage, read distribution |
| 10 | `rnaseq_qualimap` | qc | BAM + GTF | Qualimap report dir |
| 11 | `rnaseq_multiqc` | align | all QC outputs | multiqc.html |
| 12 | `rnaseq_eda` | r | count matrix + sample TSV | PCA PDF + RDS |
| 13 | `rnaseq_deseq2` | r | count matrix + de-map params | results TSV + volcano/MA PDFs |
| 14 | `rnaseq_edger` | r | count matrix + de-map params | results TSV + volcano/MA PDFs |
| 15 | `rnaseq_limma` | r | count matrix + de-map params | results TSV + volcano/MA PDFs |
| 16 | `rnaseq_enrichment` | r | DE results TSV | GO/KEGG tables + dot/bar PDFs |

### Conditional Logic

- **STAR path**: fastp → FastQC → STAR index → STAR align → featureCounts → RSeQC → Qualimap → EDA → DE → enrichment
- **Salmon path**: fastp → FastQC → Salmon index → Salmon quant → tximport → EDA → DE → enrichment
- **RSeQC + Qualimap** only fire for STAR path (no Salmon BAMs)
- **MultiQC** aggregates whatever QC outputs exist

## Test Data

### Mouse (mm10, chr19)
- Reference: chr19 FASTA subset + chr19 GTF subset from UCSC/Ensembl
- 4 paired-end libraries from a public GEO dataset, subsampled to ~50-100k read pairs
- 2 groups: ctrl (lib001, lib002) / exp (lib003, lib004)

### Human (hg38, chr22)
- Reference: chr22 FASTA subset + chr22 GTF subset (consistent with emseq/frag)
- 4 paired-end libraries from a public GEO dataset, subsampled to ~50-100k read pairs
- 2 groups: normal (lib005, lib006) / tumor (lib007, lib008)

### tools/get_test_data.sh
1. Downloads reference FASTA + GTF (full genome, then subsets to target chromosome)
2. Downloads FASTQs from SRA/ENA
3. Subsamples reads (seqtk sample)
4. Validates outputs (gzip, FASTQ format)
5. Documents provenance (GEO accessions, commands used)

### Committed to repo
- `tests/full/inputs/` — subsampled FASTQs, chromosome-subset FASTA + GTF
- `data/test-samples.tsv` — sample manifest
- Estimated ~30-40MB total

### NOT committed
- Pipeline outputs (BAMs, counts, DE results, plots, indices)
- Governed by `.gitignore`

## R Scripts

All tangled from org to `scripts/`, all use argparse for CLI args:

### rnaseq_tximport.R
- Input: quant.sf file paths, GTF path
- Output: txi.rds (tximport list object)
- Builds tx2gene from GTF, runs tximport

### rnaseq_eda.R
- Input: count matrix (featureCounts TSV or txi.rds), sample metadata TSV, design formula
- Output: PCA plot PDF, EDA RDS (design, logCPM, PCA, DGEList)
- edgeR TMM normalization, filterByExpr, PCA on logCPM

### rnaseq_deseq2.R
- Input: count matrix, sample TSV, design formula, contrast vector
- Output: results TSV (gene, log2FC, padj, baseMean), volcano PDF, MA PDF
- DESeqDataSet → DESeq() → results() with contrast

### rnaseq_edger.R
- Input: count matrix, sample TSV, design formula, contrast vector
- Output: results TSV, volcano PDF, MA PDF
- DGEList → calcNormFactors → estimateDisp → glmQLFit → glmQLFTest

### rnaseq_limma.R
- Input: count matrix, sample TSV, design formula, contrast vector
- Output: results TSV, volcano PDF, MA PDF
- voom → lmFit → contrasts.fit → eBayes → topTable

### rnaseq_enrichment.R
- Input: DE results TSV (log2FC + padj), orgdb name
- Output: GO enrichment table + dotplot PDF, KEGG enrichment table + barplot PDF
- clusterProfiler enrichGO/enrichKEGG + gseGO/gseKEGG (GSEA)

## CI/CD

### .github/workflows/generate-data.yaml
- Runs `tools/get_test_data.sh`
- Validates expected output files exist in `tests/full/inputs/`

### .github/workflows/test.yaml
- Creates conda environments
- Runs `snakemake -s workflows/test.smk --configfile config/test.yaml`
- Both STAR + Salmon paths exercised via de-map
- Verifies key output files (DE results, plots, MultiQC)

## Migration Notes

- Existing background/reference text (~1600 lines) retained and reorganized under reference headings
- Existing pipeline code (2 functional rules, R scripts) replaced entirely
- `workflow/` renamed to `workflows/`
- Singularity containers replaced with conda environments
- Old config files (int_test.yaml, env.yaml, rna_env.yaml) removed
- Old test data (quant.sf files from local paths) replaced with GEO-sourced data
