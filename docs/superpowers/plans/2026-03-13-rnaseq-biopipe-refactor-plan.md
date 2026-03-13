# RNA-seq Biopipe Module Refactor Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite rna-seq repository as a production biopipe module following emseq/frag patterns, covering FASTQ through differential expression and functional enrichment.

**Architecture:** Single modular `rnaseq.smk` (no rule all, no preamble) included by `test.smk` wrapper. Config-driven `de-map` dispatch pattern for multi-experiment DE across DESeq2/edgeR/limma. Dual alignment paths (STAR + Salmon) with method designator in filenames.

**Tech Stack:** Snakemake, conda, STAR, Salmon, fastp, FastQC, featureCounts, RSeQC, Qualimap, MultiQC, R (DESeq2, edgeR, limma, tximport, clusterProfiler, tidyverse)

**Spec:** `docs/superpowers/specs/2026-03-13-rnaseq-biopipe-refactor-design.md`

**Reference implementations:**
- `~/repos/emseq/workflows/emseq.smk` + `test.smk` + `config/test.yaml` — meth_map dispatch pattern
- `~/repos/frag/workflows/frag.smk` + `test.smk` + `config/test.yaml` — SampleTable class, directory structure
- `~/repos/frag/tools/get_test_data.sh` — test data acquisition pattern

---

## File Map

### Create (new files)
- `workflows/rnaseq.smk` — modular pipeline (16 rules, no rule all, no preamble)
- `workflows/test.smk` — test wrapper (preamble, config, de_map, rule all, symlinks, include)
- `config/test.yaml` — test configuration with de-map
- `config/rnaseq-align-conda-env.yaml` — STAR, Salmon, fastp, FastQC, samtools, subread, MultiQC
- `config/rnaseq-qc-conda-env.yaml` — RSeQC, Qualimap
- `config/rnaseq-r-env.yaml` — DESeq2, edgeR, limma, tximport, clusterProfiler, tidyverse
- `scripts/rnaseq_tximport.R` — Salmon quant.sf → txi.rds
- `scripts/rnaseq_eda.R` — count matrix → PCA + EDA
- `scripts/rnaseq_deseq2.R` — DE via DESeq2
- `scripts/rnaseq_edger.R` — DE via edgeR
- `scripts/rnaseq_limma.R` — DE via limma-voom
- `scripts/rnaseq_enrichment.R` — GO/KEGG/GSEA via clusterProfiler
- `tools/get_test_data.sh` — download + subset mouse/human test data from GEO/SRA
- `tools/shell/org_update.sh` — org tangle + README export
- `data/test-samples.tsv` — 8-library sample manifest (4 mouse, 4 human)
- `.github/workflows/generate-data.yaml` — CI: test data validation
- `.github/workflows/test.yaml` — CI: full pipeline test
- `CLAUDE.md` — project-specific Claude instructions
- `tests/full/inputs/` — committed test data (FASTQs, FASTA subsets, GTF subsets)

### Remove (old files to clean up)
- `workflow/` directory (old name, replaced by `workflows/`)
- `config/int_test.yaml`, `config/env.yaml`, `config/rna_env.yaml`
- `scripts/make_salmon_txi.R`, `scripts/all_rna_eda.R`, `scripts/dumbtest.R`, `scripts/symlink_salmon.sh`
- `test/` directory (old test data)
- `test.smk` (root-level)

### Modify
- `rna-seq.org` — restructure headings, add pipeline src blocks, retain background text
- `.gitignore` — add tests/full/* with !tests/full/inputs/ pattern
- `README.md` — will be regenerated from org export

---

## Chunk 1: Repository Scaffolding and Conda Environments

### Task 1: Clean up old files and create directory structure

**Files:**
- Remove: `workflow/`, `config/int_test.yaml`, `config/env.yaml`, `config/rna_env.yaml`, `scripts/make_salmon_txi.R`, `scripts/all_rna_eda.R`, `scripts/dumbtest.R`, `scripts/symlink_salmon.sh`, `test/`, `test.smk`
- Create: `workflows/`, `tools/shell/`, `data/`, `tests/full/inputs/`, `.github/workflows/`, `resources/`

- [ ] **Step 1: Remove old pipeline files**

```bash
rm -rf workflow/ test/ test.smk
rm -f config/int_test.yaml config/env.yaml config/rna_env.yaml
rm -f scripts/make_salmon_txi.R scripts/all_rna_eda.R scripts/dumbtest.R scripts/symlink_salmon.sh
```

- [ ] **Step 2: Create new directory structure**

```bash
mkdir -p workflows tools/shell data tests/full/inputs .github/workflows resources
```

- [ ] **Step 3: Update .gitignore**

Add the frag-style test data pattern to `.gitignore`:
```
# --- test data (auto-managed by get_test_data.sh) ---
tests/full/*
!tests/full/inputs/
!tests/full/inputs/**
```

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "Remove old pipeline files and create biopipe directory structure"
```

### Task 2: Write conda environment YAMLs

**Files:**
- Create: `config/rnaseq-align-conda-env.yaml`
- Create: `config/rnaseq-qc-conda-env.yaml`
- Create: `config/rnaseq-r-env.yaml`

- [ ] **Step 1: Write alignment + FASTQ processing env**

```yaml
# config/rnaseq-align-conda-env.yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - fastp
  - fastqc
  - multiqc
  - salmon
  - samtools
  - star
  - subread
```

- [ ] **Step 2: Write BAM QC env**

```yaml
# config/rnaseq-qc-conda-env.yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - qualimap
  - rseqc
  - ucsc-gtftogenepred
```

- [ ] **Step 3: Write R env**

```yaml
# config/rnaseq-r-env.yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - bioconductor-clusterprofiler
  - bioconductor-deseq2
  - bioconductor-edger
  - bioconductor-limma
  - bioconductor-org.hs.eg.db
  - bioconductor-org.mm.eg.db
  - bioconductor-tximport
  - r-argparse
  - r-cowplot
  - r-ggrepel
  - r-pheatmap
  - r-tidyverse
```

- [ ] **Step 4: Commit**

```bash
git add config/rnaseq-align-conda-env.yaml config/rnaseq-qc-conda-env.yaml config/rnaseq-r-env.yaml
git commit -m "Add conda environment YAMLs for alignment, QC, and R"
```

### Task 3: Write test configuration YAML

**Files:**
- Create: `config/test.yaml`

- [ ] **Step 1: Write config/test.yaml**

```yaml
# config/test.yaml
# Test pipeline configuration for RNA-seq biopipe module.
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

main-data-dir: tests/full

sample-tsv-path: ${HOME}/repos/rna-seq/data/test-samples.tsv

refs:
  mm10_chr19:
    fasta: tests/full/inputs/mm10_chr19.fa.gz
    gtf: tests/full/inputs/mm10_chr19.gtf.gz
    star_sa_bases: 11
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

- [ ] **Step 2: Write data/test-samples.tsv**

```tsv
library_id	r1_basename	r2_basename	group	organism
lib001	lib001_R1.fastq.gz	lib001_R2.fastq.gz	ctrl	mouse
lib002	lib002_R1.fastq.gz	lib002_R2.fastq.gz	ctrl	mouse
lib003	lib003_R1.fastq.gz	lib003_R2.fastq.gz	exp	mouse
lib004	lib004_R1.fastq.gz	lib004_R2.fastq.gz	exp	mouse
lib005	lib005_R1.fastq.gz	lib005_R2.fastq.gz	normal	human
lib006	lib006_R1.fastq.gz	lib006_R2.fastq.gz	normal	human
lib007	lib007_R1.fastq.gz	lib007_R2.fastq.gz	tumor	human
lib008	lib008_R1.fastq.gz	lib008_R2.fastq.gz	tumor	human
```

- [ ] **Step 3: Commit**

```bash
git add config/test.yaml data/test-samples.tsv
git commit -m "Add test configuration YAML and sample manifest"
```

### Task 4: Write tools/get_test_data.sh

**Files:**
- Create: `tools/get_test_data.sh`

Follow the frag pattern (`~/repos/frag/tools/get_test_data.sh`): parse_args, clean_inputs_dir, ensure_gitignore, fetch functions, main.

- [ ] **Step 1: Write get_test_data.sh**

The script must:
1. Download mouse chr19 FASTA subset from UCSC + chr19 GTF from Ensembl/GENCODE
2. Download human chr22 FASTA subset from UCSC + chr22 GTF from GENCODE
3. Download 8 paired-end RNA-seq FASTQ files from SRA (4 mouse, 4 human)
   - Mouse: pick a simple treatment/control GEO dataset (e.g., GSE series with ≥4 samples)
   - Human: pick a tumor/normal GEO dataset
4. Subsample each FASTQ to ~50k reads using `fastq-dump -X 50000 --split-files --gzip`
5. Subset GTF to target chromosome using `awk '$1 == "chr19"'` / `awk '$1 == "chr22"'`
6. Validate all outputs
7. Rename files to lib001-lib008 convention

Structure (following frag's get_test_data.sh):
```bash
#!/usr/bin/env bash
set -euo pipefail

parse_args() {
    declare -g script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    declare -g repo_dir="$(cd "${script_dir}/.." && pwd)"
    declare -g out_dir="${repo_dir}/tests/full/inputs"

    declare -g nreads=50000
    declare -g fa_head_lines=4000000

    # Mouse refs
    declare -g mm_fa_url="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr19.fa.gz"
    declare -g mm_gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
    declare -g out_mm_fa="${out_dir}/mm10_chr19.fa.gz"
    declare -g out_mm_gtf="${out_dir}/mm10_chr19.gtf.gz"

    # Human refs
    declare -g hs_fa_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
    declare -g hs_gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"
    declare -g out_hs_fa="${out_dir}/hg38_chr22.fa.gz"
    declare -g out_hs_gtf="${out_dir}/hg38_chr22.gtf.gz"

    # Mouse SRA accessions (pick a suitable GEO dataset at implementation time)
    # These are placeholders — replace with actual accessions
    declare -g MOUSE_RUNS=( SRRxxxxxxx SRRxxxxxxx SRRxxxxxxx SRRxxxxxxx )
    declare -g MOUSE_LIBS=( lib001 lib002 lib003 lib004 )

    # Human SRA accessions (pick a suitable GEO dataset at implementation time)
    declare -g HUMAN_RUNS=( SRRxxxxxxx SRRxxxxxxx SRRxxxxxxx SRRxxxxxxx )
    declare -g HUMAN_LIBS=( lib005 lib006 lib007 lib008 )
}

main() {
    parse_args "$@"
    clean_inputs_dir
    ensure_gitignore
    fetch_mouse_ref
    fetch_human_ref
    fetch_mouse_reads
    fetch_human_reads
    validate_outputs
}
```

**IMPORTANT:** The actual SRA accession numbers must be chosen at implementation time. Use NCBI GEO search or SRA Run Selector to find:
- Mouse: paired-end RNA-seq, simple 2-group comparison (e.g., treatment vs control), ≥4 samples, Illumina platform, mouse genome (mm10/GRCm38)
- Human: paired-end RNA-seq, tumor vs normal, ≥4 samples, Illumina platform, human genome (hg38/GRCh38)
- Both should use standard library prep (polyA enrichment or rRNA depletion)
- Verify accessions are valid: `fastq-dump --split-files --gzip -X 10 SRRxxxxxxx` should succeed

**GEO search strategy:**
1. Go to https://www.ncbi.nlm.nih.gov/gds
2. Search: `"RNA-seq"[DataSet Type] AND "Mus musculus"[Organism] AND "paired"[Filter]`
3. Look for simple 2-condition experiments with 4-6 samples
4. Click through to SRA Run Selector to get SRR accessions
5. Pick 4 runs (2 per condition)

If network access is unavailable, create minimal synthetic FASTQ files for testing:
```bash
# Generate minimal synthetic paired-end FASTQs (for testing pipeline logic only)
for lib in lib001 lib002 lib003 lib004 lib005 lib006 lib007 lib008; do
  for read in R1 R2; do
    python3 -c "
import random; random.seed(42)
bases='ACGT'
for i in range(1000):
    seq=''.join(random.choice(bases) for _ in range(100))
    qual='I'*100
    print(f'@{lib}_{read}_{i}')
    print(seq)
    print('+')
    print(qual)
" | gzip > "tests/full/inputs/${lib}_${read}.fastq.gz"
  done
done
```
This produces syntactically valid FASTQs that will pass through fastp/FastQC but won't align meaningfully. Good enough for DAG validation.

Each fetch function subsamples to `$nreads` reads, renames to lib00N convention, validates gzip headers.

- [ ] **Step 2: Make executable and commit**

```bash
chmod +x tools/get_test_data.sh
git add tools/get_test_data.sh
git commit -m "Add test data acquisition script (GEO/SRA provenance)"
```

### Task 5: Write tools/shell/org_update.sh

**Files:**
- Create: `tools/shell/org_update.sh`

Follow frag's pattern (`~/repos/frag/tools/shell/org_update.sh`).

- [ ] **Step 1: Write org_update.sh**

```bash
#!/usr/bin/env bash
set -euo pipefail

# Tangle rna-seq.org and export README
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ORG_FILE="${REPO_DIR}/rna-seq.org"

emacsclient --socket-name ~/.emacs.d/server/server \
  --eval "(progn
    (find-file \"${ORG_FILE}\")
    (revert-buffer t t)
    (org-babel-tangle)
    (org-md-export-to-markdown))"
```

- [ ] **Step 2: Make executable and commit**

```bash
chmod +x tools/shell/org_update.sh
git add tools/shell/org_update.sh
git commit -m "Add org update script for tangle + README export"
```

## Chunk 2: Snakemake Workflows — FASTQ QC + Alignment

### Task 6: Write workflows/rnaseq.smk preamble and FASTQ QC rules

**Files:**
- Create: `workflows/rnaseq.smk`

- [ ] **Step 1: Write modular header + rnaseq_fastp + rnaseq_fastqc**

```python
# workflows/rnaseq.smk
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
#
# Source:  rna-seq.org
# ============================================================
#
# This is a modular snakefile, intended to be incorporated into a larger
# workflow using the "include:" directive.
# (See https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html)
#

rule rnaseq_fastp:
    message: "RNA-seq fastp FASTQ processing"
    conda: ENV_ALIGN
    input:
        r1 = f"{D_RNASEQ}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_RNASEQ}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_rnaseq_fastp.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_rnaseq_fastp.tsv"
    threads: 8
    output:
        r1   = f"{D_RNASEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2   = f"{D_RNASEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        json = f"{D_RNASEQ}/qc/{{library_id}}_fastp.json",
        html = f"{D_RNASEQ}/qc/{{library_id}}_fastp.html",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[fastp] $(date) lib={wildcards.library_id} threads={threads}"

        fastp \
          --detect_adapter_for_pe \
          --in1 "{input.r1}" --in2 "{input.r2}" \
          --out1 "{output.r1}" --out2 "{output.r2}" \
          --json "{output.json}" --html "{output.html}" \
          --thread {threads}
        """

rule rnaseq_fastqc:
    message: "RNA-seq FastQC quality check"
    conda: ENV_ALIGN
    input:
        fq = f"{D_RNASEQ}/fastqs/{{library_id}}.{{processing}}_{{read}}.fastq.gz",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{processing}}_{{read}}_rnaseq_fastqc.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{processing}}_{{read}}_rnaseq_fastqc.tsv"
    params:
        out_dir = f"{D_RNASEQ}/qc",
    threads: 2
    output:
        html = f"{D_RNASEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.html",
        zip  = f"{D_RNASEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.zip",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[fastqc] $(date) lib={wildcards.library_id} {wildcards.processing} {wildcards.read}"

        fastqc "{input.fq}" \
          --outdir "{params.out_dir}" \
          --threads {threads}
        """
```

- [ ] **Step 2: Commit**

```bash
git add workflows/rnaseq.smk
git commit -m "Add rnaseq.smk with fastp and fastqc rules"
```

### Task 7: Add STAR alignment rules to rnaseq.smk

**Files:**
- Modify: `workflows/rnaseq.smk`

- [ ] **Step 1: Append STAR index, align, and featureCounts rules**

```python
rule rnaseq_star_index:
    message: "Build STAR genome index"
    conda: ENV_ALIGN
    input:
        fasta = lambda wc: config["refs"][wc.ref_name]["fasta"],
        gtf   = lambda wc: config["refs"][wc.ref_name]["gtf"],
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_rnaseq_star_index.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_rnaseq_star_index.tsv"
    params:
        genome_dir    = f"{D_RNASEQ}/ref/star/{{ref_name}}",
        sa_index_bases = lambda wc: config["refs"][wc.ref_name].get("star_sa_bases", 14),
    threads: 8
    output:
        sa = f"{D_RNASEQ}/ref/star/{{ref_name}}/SA",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[STAR index] $(date) ref={wildcards.ref_name} threads={threads}"

        mkdir -p "{params.genome_dir}"

        STAR --runMode genomeGenerate \
          --genomeDir "{params.genome_dir}" \
          --genomeFastaFiles "{input.fasta}" \
          --sjdbGTFfile "{input.gtf}" \
          --genomeSAindexNbases {params.sa_index_bases} \
          --runThreadN {threads}
        """

rule rnaseq_star_align:
    message: "STAR alignment of RNA-seq reads"
    conda: ENV_ALIGN
    input:
        r1        = f"{D_RNASEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2        = f"{D_RNASEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        star_idx  = f"{D_RNASEQ}/ref/star/{{ref_name}}/SA",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_rnaseq_star_align.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_rnaseq_star_align.tsv"
    params:
        genome_dir = f"{D_RNASEQ}/ref/star/{{ref_name}}",
        prefix     = f"{D_RNASEQ}/bams/{{library_id}}.{{ref_name}}.star.",
    threads: 8
    output:
        bam = f"{D_RNASEQ}/bams/{{library_id}}.{{ref_name}}.star.Aligned.sortedByCoord.out.bam",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[STAR align] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"

        STAR --runMode alignReads \
          --genomeDir "{params.genome_dir}" \
          --readFilesIn "{input.r1}" "{input.r2}" \
          --readFilesCommand zcat \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix "{params.prefix}" \
          --runThreadN {threads}

        samtools index "{output.bam}"
        """

rule rnaseq_featurecounts:
    message: "featureCounts gene-level quantification"
    conda: ENV_ALIGN
    input:
        bam = f"{D_RNASEQ}/bams/{{library_id}}.{{ref_name}}.star.Aligned.sortedByCoord.out.bam",
        gtf = lambda wc: config["refs"][wc.ref_name]["gtf"],
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_rnaseq_featurecounts.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_rnaseq_featurecounts.tsv"
    threads: 4
    output:
        counts  = f"{D_RNASEQ}/counts/{{library_id}}.{{ref_name}}.star.featurecounts.tsv",
        summary = f"{D_RNASEQ}/counts/{{library_id}}.{{ref_name}}.star.featurecounts.tsv.summary",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[featureCounts] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"

        featureCounts \
          -a "{input.gtf}" \
          -o "{output.counts}" \
          -p --countReadPairs \
          -T {threads} \
          "{input.bam}"
        """
```

- [ ] **Step 2: Commit**

```bash
git add workflows/rnaseq.smk
git commit -m "Add STAR index, alignment, and featureCounts rules"
```

### Task 8: Add Salmon alignment rules to rnaseq.smk

**Files:**
- Modify: `workflows/rnaseq.smk`

- [ ] **Step 1: Append Salmon index and quant rules**

```python
rule rnaseq_salmon_index:
    message: "Build Salmon index with decoy-aware gentrome"
    conda: ENV_ALIGN
    input:
        fasta = lambda wc: config["refs"][wc.ref_name]["fasta"],
        gtf   = lambda wc: config["refs"][wc.ref_name]["gtf"],
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_rnaseq_salmon_index.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_rnaseq_salmon_index.tsv"
    params:
        index_dir = f"{D_RNASEQ}/ref/salmon/{{ref_name}}",
    threads: 8
    output:
        info = f"{D_RNASEQ}/ref/salmon/{{ref_name}}/info.json",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[Salmon index] $(date) ref={wildcards.ref_name} threads={threads}"

        tmpdir=$(mktemp -d)

        # Extract transcriptome from GTF + genome FASTA using gffread
        # For gentrome: concatenate transcriptome + genome, use genome as decoy
        grep "^>" <(zcat "{input.fasta}" 2>/dev/null || cat "{input.fasta}") \
          | cut -d " " -f 1 | sed 's/>//' > "$tmpdir/decoys.txt"

        # Build gentrome (transcripts + genome)
        gffread -w "$tmpdir/transcripts.fa" -g "{input.fasta}" "{input.gtf}"
        cat "$tmpdir/transcripts.fa" <(zcat "{input.fasta}" 2>/dev/null || cat "{input.fasta}") \
          > "$tmpdir/gentrome.fa"

        salmon index \
          -t "$tmpdir/gentrome.fa" \
          -d "$tmpdir/decoys.txt" \
          -i "{params.index_dir}" \
          -p {threads}

        rm -rf "$tmpdir"
        """

rule rnaseq_salmon_quant:
    message: "Salmon quasi-mapping quantification"
    conda: ENV_ALIGN
    input:
        r1    = f"{D_RNASEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2    = f"{D_RNASEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        index = f"{D_RNASEQ}/ref/salmon/{{ref_name}}/info.json",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_rnaseq_salmon_quant.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_rnaseq_salmon_quant.tsv"
    params:
        index_dir = f"{D_RNASEQ}/ref/salmon/{{ref_name}}",
        out_dir   = f"{D_RNASEQ}/quant/{{library_id}}.{{ref_name}}.salmon",
    threads: 8
    output:
        quant = f"{D_RNASEQ}/quant/{{library_id}}.{{ref_name}}.salmon/quant.sf",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[Salmon quant] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"

        salmon quant \
          -i "{params.index_dir}" \
          -l A \
          -1 "{input.r1}" -2 "{input.r2}" \
          -o "{params.out_dir}" \
          -p {threads} \
          --validateMappings
        """
```

- [ ] **Step 2: Commit**

```bash
git add workflows/rnaseq.smk
git commit -m "Add Salmon index and quantification rules"
```

## Chunk 3: R Scripts — tximport + EDA

### Task 9: Write scripts/rnaseq_tximport.R

**Files:**
- Create: `scripts/rnaseq_tximport.R`

- [ ] **Step 1: Write tximport R script**

```r
#!/usr/bin/env Rscript
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

library(argparse)
library(tximport)

parser <- ArgumentParser(description = "Import Salmon quant.sf files via tximport")
parser$add_argument("--quant-files", nargs = "+", required = TRUE,
                    help = "Space-separated Salmon quant.sf file paths")
parser$add_argument("--gtf", required = TRUE,
                    help = "GTF annotation file for tx2gene mapping")
parser$add_argument("--out-rds", required = TRUE,
                    help = "Output RDS file path for tximport object")
args <- parser$parse_args()

# Build tx2gene from GTF
gtf <- rtracklayer::import(args$gtf)
tx2gene <- unique(data.frame(
  TXNAME = gtf$transcript_id[!is.na(gtf$transcript_id)],
  GENEID = gtf$gene_id[!is.na(gtf$transcript_id)]
))

# Name quant files by library ID (extract from path)
quant_files <- args$quant_files
names(quant_files) <- gsub(".*/(lib[0-9]+)\\..*", "\\1", quant_files)

# Run tximport
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)

# Save
saveRDS(txi, file = args$out_rds)
```

- [ ] **Step 2: Add tximport rule to rnaseq.smk**

```python
rule rnaseq_tximport:
    message: "Import Salmon quantification via tximport"
    conda: ENV_R
    input:
        quant_files = lambda wc: expand(
            f"{D_RNASEQ}/quant/{{library_id}}.{{ref_name}}.salmon/quant.sf",
            library_id=de_map[wc.experiment]["libs"],
            ref_name=de_map[wc.experiment]["ref_name"],
        ),
        gtf = lambda wc: config["refs"][de_map[wc.experiment]["ref_name"]]["gtf"],
    log:
        cmd = f"{D_LOGS}/{{experiment}}_rnaseq_tximport.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_rnaseq_tximport.tsv"
    params:
        script = f"{R_RNASEQ}/scripts/rnaseq_tximport.R",
    threads: 1
    output:
        rds = f"{D_RNASEQ}/counts/{{experiment}}.salmon.txi.rds",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[tximport] $(date) experiment={wildcards.experiment}"

        Rscript {params.script} \
          --quant-files {input.quant_files} \
          --gtf "{input.gtf}" \
          --out-rds "{output.rds}"
        """
```

- [ ] **Step 3: Commit**

```bash
git add scripts/rnaseq_tximport.R workflows/rnaseq.smk
git commit -m "Add tximport R script and Snakemake rule"
```

### Task 10: Write scripts/rnaseq_eda.R and rule

**Files:**
- Create: `scripts/rnaseq_eda.R`
- Modify: `workflows/rnaseq.smk`

- [ ] **Step 1: Write EDA R script**

```r
#!/usr/bin/env Rscript
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

library(argparse)
library(tidyverse)
library(edgeR)
library(cowplot)
library(ggrepel)

parser <- ArgumentParser(description = "RNA-seq exploratory data analysis")
parser$add_argument("--counts", required = TRUE, help = "Count matrix path (TSV or RDS)")
parser$add_argument("--input-type", required = TRUE,
                    choices = c("featurecounts", "tximport"),
                    help = "Input format type")
parser$add_argument("--sample-tsv", required = TRUE, help = "Sample metadata TSV")
parser$add_argument("--design", required = TRUE, help = "Design formula string")
parser$add_argument("--out-pca", required = TRUE, help = "PCA plot PDF output")
parser$add_argument("--out-rds", required = TRUE, help = "EDA results RDS output")
args <- parser$parse_args()

samples <- read_tsv(args$sample_tsv, show_col_types = FALSE)

if (args$input_type == "featurecounts") {
  # Read featureCounts output (skip first line comment, gene ID in col 1, counts in col 7+)
  fc <- read.delim(args$counts, comment.char = "#")
  counts <- as.matrix(fc[, 7:ncol(fc)])
  rownames(counts) <- fc$Geneid
  y <- DGEList(counts = counts)
} else {
  txi <- readRDS(args$counts)
  counts <- txi$counts
  # Apply tximport length offset
  norm_mat <- txi$length
  norm_mat <- norm_mat / exp(rowMeans(log(norm_mat)))
  norm_counts <- counts / norm_mat
  eff_lib <- calcNormFactors(norm_counts) * colSums(norm_counts)
  norm_mat <- sweep(norm_mat, 2, eff_lib, "*")
  norm_mat <- log(norm_mat)
  y <- DGEList(counts)
  y <- scaleOffset(y, norm_mat)
}

# Build design matrix
formula <- as.formula(args$design)
design <- model.matrix(formula, data = samples)

# Filter lowly expressed genes
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize and compute logCPM
y <- calcNormFactors(y)
logCPM <- cpm(y, prior.count = 2, log = TRUE)

# PCA
pca <- prcomp(t(logCPM))
pve_pc1 <- round(100 * summary(pca)$importance[2, 1])
pve_pc2 <- round(100 * summary(pca)$importance[2, 2])

# Determine color factor (last term in design formula)
factor_vec <- all.vars(formula)
color_var <- tail(factor_vec, 1)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("library_id") %>%
  left_join(samples, by = "library_id")

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_var]], label = library_id)) +
  geom_point(size = 4) +
  geom_text_repel() +
  xlab(paste0("PC1 (", pve_pc1, "% variance)")) +
  ylab(paste0("PC2 (", pve_pc2, "% variance)")) +
  scale_color_discrete(name = color_var) +
  theme_cowplot() +
  theme(legend.position = "bottom")

# Save
save_plot(pca_plot, filename = args$out_pca)
saveRDS(list(design = design, logCPM = logCPM, pca = pca, y = y),
        file = args$out_rds)
```

- [ ] **Step 2: Add EDA rule to rnaseq.smk**

```python
def get_count_matrix(wc, de_map):
    """Resolve count matrix path based on alignment method."""
    if wc.align_method == "star":
        libs = de_map[wc.experiment]["libs"]
        ref = de_map[wc.experiment]["ref_name"]
        return [f"{D_RNASEQ}/counts/{lib}.{ref}.star.featurecounts.tsv" for lib in libs]
    elif wc.align_method == "salmon":
        return f"{D_RNASEQ}/counts/{wc.experiment}.salmon.txi.rds"

rule rnaseq_eda:
    message: "RNA-seq exploratory data analysis"
    conda: ENV_R
    input:
        counts     = lambda wc: get_count_matrix(wc, de_map),
        sample_tsv = SAMPLE_TSV,
    log:
        cmd = f"{D_LOGS}/{{experiment}}_{{align_method}}_rnaseq_eda.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_{{align_method}}_rnaseq_eda.tsv"
    params:
        script     = f"{R_RNASEQ}/scripts/rnaseq_eda.R",
        input_type = lambda wc: "featurecounts" if wc.align_method == "star" else "tximport",
        design     = lambda wc: de_map[wc.experiment]["design"],
    threads: 1
    output:
        pca = f"{D_RNASEQ}/eda/{{experiment}}.{{align_method}}.pca.pdf",
        rds = f"{D_RNASEQ}/eda/{{experiment}}.{{align_method}}.eda.rds",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[EDA] $(date) experiment={wildcards.experiment} method={wildcards.align_method}"

        Rscript {params.script} \
          --counts "{input.counts}" \
          --input-type "{params.input_type}" \
          --sample-tsv "{input.sample_tsv}" \
          --design "{params.design}" \
          --out-pca "{output.pca}" \
          --out-rds "{output.rds}"
        """
```

Note: The `get_count_matrix` helper for STAR returns a list of per-library featureCounts files. The R script needs to handle merging them. Actually, it's cleaner to have a separate rule that merges per-library featureCounts into an experiment-level count matrix first. Add a `rnaseq_merge_featurecounts` rule:

```python
rule rnaseq_merge_featurecounts:
    message: "Merge per-library featureCounts into experiment count matrix"
    conda: ENV_R
    input:
        counts = lambda wc: expand(
            f"{D_RNASEQ}/counts/{{library_id}}.{{ref_name}}.star.featurecounts.tsv",
            library_id=de_map[wc.experiment]["libs"],
            ref_name=de_map[wc.experiment]["ref_name"],
        ),
    log:
        cmd = f"{D_LOGS}/{{experiment}}_rnaseq_merge_featurecounts.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_rnaseq_merge_featurecounts.tsv"
    threads: 1
    output:
        tsv = f"{D_RNASEQ}/counts/{{experiment}}.star.featurecounts.tsv",
    run:
        import pandas as pd
        frames = []
        for f in input.counts:
            df = pd.read_csv(f, sep="\t", comment="#")
            lib_col = df.columns[-1]  # last column is the count
            lib_id = lib_col.split("/")[-1].split(".")[0]  # extract library_id
            frames.append(df[["Geneid", lib_col]].rename(columns={lib_col: lib_id}))
        merged = frames[0]
        for df in frames[1:]:
            merged = merged.merge(df, on="Geneid")
        merged.to_csv(output.tsv, sep="\t", index=False)
```

Then update `get_count_matrix`:
```python
def get_count_matrix(wc, de_map):
    if wc.align_method == "star":
        return f"{D_RNASEQ}/counts/{wc.experiment}.star.featurecounts.tsv"
    elif wc.align_method == "salmon":
        return f"{D_RNASEQ}/counts/{wc.experiment}.salmon.txi.rds"
```

- [ ] **Step 3: Commit**

```bash
git add scripts/rnaseq_eda.R workflows/rnaseq.smk
git commit -m "Add EDA R script, merge featurecounts rule, and EDA snakemake rule"
```

## Chunk 4: BAM QC + MultiQC

### Task 11: Add RSeQC, Qualimap, and MultiQC rules

**Files:**
- Modify: `workflows/rnaseq.smk`

- [ ] **Step 1: Append BAM QC and MultiQC rules**

```python
rule rnaseq_rseqc:
    message: "RSeQC gene body coverage and read distribution"
    conda: ENV_QC
    input:
        bam = f"{D_RNASEQ}/bams/{{library_id}}.{{ref_name}}.star.Aligned.sortedByCoord.out.bam",
        gtf = lambda wc: config["refs"][wc.ref_name]["gtf"],
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_rnaseq_rseqc.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_rnaseq_rseqc.tsv"
    params:
        prefix = f"{D_RNASEQ}/qc/rseqc/{{library_id}}.{{ref_name}}",
        bed    = f"{D_RNASEQ}/ref/{{ref_name}}.rseqc.bed",
    threads: 1
    output:
        read_dist = f"{D_RNASEQ}/qc/rseqc/{{library_id}}.{{ref_name}}.read_distribution.txt",
        gene_body = f"{D_RNASEQ}/qc/rseqc/{{library_id}}.{{ref_name}}.geneBodyCoverage.txt",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[RSeQC] $(date) lib={wildcards.library_id} ref={wildcards.ref_name}"

        mkdir -p "$(dirname "{params.prefix}")"

        # Convert GTF to BED12 for RSeQC
        gtfToGenePred "{input.gtf}" /dev/stdout | genePredToBed /dev/stdin "{params.bed}" || \
          awk '$3=="transcript"' "{input.gtf}" | gtf2bed > "{params.bed}"

        read_distribution.py -i "{input.bam}" -r "{params.bed}" \
          > "{output.read_dist}"

        geneBody_coverage.py -i "{input.bam}" -r "{params.bed}" \
          -o "{params.prefix}"
        """

rule rnaseq_qualimap:
    message: "Qualimap RNA-seq QC"
    conda: ENV_QC
    input:
        bam = f"{D_RNASEQ}/bams/{{library_id}}.{{ref_name}}.star.Aligned.sortedByCoord.out.bam",
        gtf = lambda wc: config["refs"][wc.ref_name]["gtf"],
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_rnaseq_qualimap.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_rnaseq_qualimap.tsv"
    params:
        out_dir = f"{D_RNASEQ}/qc/qualimap/{{library_id}}.{{ref_name}}",
    threads: 4
    output:
        report = f"{D_RNASEQ}/qc/qualimap/{{library_id}}.{{ref_name}}/qualimapReport.html",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[Qualimap] $(date) lib={wildcards.library_id} ref={wildcards.ref_name}"

        qualimap rnaseq \
          -bam "{input.bam}" \
          -gtf "{input.gtf}" \
          --outdir "{params.out_dir}" \
          --java-mem-size=4G
        """

rule rnaseq_multiqc:
    message: "MultiQC aggregation of QC reports"
    conda: ENV_ALIGN
    input:
        qc_dir = f"{D_RNASEQ}/qc",
    log:
        cmd = f"{D_LOGS}/rnaseq_multiqc.log",
    benchmark:
        f"{D_BENCHMARK}/rnaseq_multiqc.tsv"
    threads: 1
    output:
        html = f"{D_RNASEQ}/qc/multiqc.html",
    params:
        out_dir = f"{D_RNASEQ}/qc",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[MultiQC] $(date)"

        multiqc "{input.qc_dir}" \
          --outdir "{params.out_dir}" \
          --filename multiqc.html \
          --force
        """
```

- [ ] **Step 2: Commit**

```bash
git add workflows/rnaseq.smk
git commit -m "Add RSeQC, Qualimap, and MultiQC rules"
```

## Chunk 5: Differential Expression R Scripts + Rules

### Task 12: Write scripts/rnaseq_deseq2.R and rule

**Files:**
- Create: `scripts/rnaseq_deseq2.R`
- Modify: `workflows/rnaseq.smk`

- [ ] **Step 1: Write DESeq2 R script**

```r
#!/usr/bin/env Rscript
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

library(argparse)
library(tidyverse)
library(DESeq2)
library(tximport)

parser <- ArgumentParser(description = "Differential expression with DESeq2")
parser$add_argument("--counts", required = TRUE, help = "Count matrix (TSV or RDS)")
parser$add_argument("--input-type", required = TRUE, choices = c("featurecounts", "tximport"))
parser$add_argument("--sample-tsv", required = TRUE, help = "Sample metadata TSV")
parser$add_argument("--design", required = TRUE, help = "Design formula string")
parser$add_argument("--contrast", nargs = 3, required = TRUE,
                    help = "Contrast: factor level_test level_ref")
parser$add_argument("--out-tsv", required = TRUE, help = "Results TSV output")
parser$add_argument("--out-volcano", required = TRUE, help = "Volcano plot PDF")
parser$add_argument("--out-ma", required = TRUE, help = "MA plot PDF")
args <- parser$parse_args()

samples <- read_tsv(args$sample_tsv, show_col_types = FALSE)
formula <- as.formula(args$design)

if (args$input_type == "featurecounts") {
  cts <- read.delim(args$counts)
  count_mat <- as.matrix(cts[, -1])
  rownames(count_mat) <- cts[[1]]
  dds <- DESeqDataSetFromMatrix(countData = count_mat,
                                colData = samples,
                                design = formula)
} else {
  txi <- readRDS(args$counts)
  dds <- DESeqDataSetFromTximport(txi, colData = samples, design = formula)
}

dds <- DESeq(dds)
res <- results(dds, contrast = args$contrast)
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

write_tsv(res_df, args$out_tsv)

# Volcano plot
pdf(args$out_volcano)
plot(res_df$log2FoldChange, -log10(res_df$padj),
     pch = 20, cex = 0.5,
     xlab = "log2 Fold Change", ylab = "-log10(padj)",
     main = paste("DESeq2:", args$contrast[2], "vs", args$contrast[3]))
abline(h = -log10(0.05), col = "red", lty = 2)
dev.off()

# MA plot
pdf(args$out_ma)
plotMA(res, main = paste("DESeq2:", args$contrast[2], "vs", args$contrast[3]))
dev.off()
```

- [ ] **Step 2: Add DESeq2 rule to rnaseq.smk**

```python
rule rnaseq_deseq2:
    message: "Differential expression with DESeq2"
    conda: ENV_R
    input:
        counts     = lambda wc: get_count_matrix(wc, de_map),
        sample_tsv = SAMPLE_TSV,
    log:
        cmd = f"{D_LOGS}/{{experiment}}_{{align_method}}_rnaseq_deseq2.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_{{align_method}}_rnaseq_deseq2.tsv"
    params:
        script     = f"{R_RNASEQ}/scripts/rnaseq_deseq2.R",
        input_type = lambda wc: "featurecounts" if wc.align_method == "star" else "tximport",
        design     = lambda wc: de_map[wc.experiment]["design"],
        contrast   = lambda wc: " ".join(de_map[wc.experiment]["contrast"]),
    threads: 1
    output:
        tsv     = f"{D_RNASEQ}/de/{{experiment}}.{{align_method}}.deseq2_results.tsv",
        volcano = f"{D_RNASEQ}/de/{{experiment}}.{{align_method}}.deseq2_volcano.pdf",
        ma      = f"{D_RNASEQ}/de/{{experiment}}.{{align_method}}.deseq2_ma.pdf",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[DESeq2] $(date) experiment={wildcards.experiment} method={wildcards.align_method}"

        Rscript {params.script} \
          --counts "{input.counts}" \
          --input-type "{params.input_type}" \
          --sample-tsv "{input.sample_tsv}" \
          --design "{params.design}" \
          --contrast {params.contrast} \
          --out-tsv "{output.tsv}" \
          --out-volcano "{output.volcano}" \
          --out-ma "{output.ma}"
        """
```

- [ ] **Step 3: Commit**

```bash
git add scripts/rnaseq_deseq2.R workflows/rnaseq.smk
git commit -m "Add DESeq2 R script and Snakemake rule"
```

### Task 13: Write scripts/rnaseq_edger.R and rule

**Files:**
- Create: `scripts/rnaseq_edger.R`
- Modify: `workflows/rnaseq.smk`

- [ ] **Step 1: Write edgeR R script**

```r
#!/usr/bin/env Rscript
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

library(argparse)
library(tidyverse)
library(edgeR)

parser <- ArgumentParser(description = "Differential expression with edgeR")
parser$add_argument("--counts", required = TRUE)
parser$add_argument("--input-type", required = TRUE, choices = c("featurecounts", "tximport"))
parser$add_argument("--sample-tsv", required = TRUE)
parser$add_argument("--design", required = TRUE)
parser$add_argument("--contrast", nargs = 3, required = TRUE)
parser$add_argument("--out-tsv", required = TRUE)
parser$add_argument("--out-volcano", required = TRUE)
parser$add_argument("--out-ma", required = TRUE)
args <- parser$parse_args()

samples <- read_tsv(args$sample_tsv, show_col_types = FALSE)
formula <- as.formula(args$design)

if (args$input_type == "featurecounts") {
  cts <- read.delim(args$counts)
  count_mat <- as.matrix(cts[, -1])
  rownames(count_mat) <- cts[[1]]
  y <- DGEList(counts = count_mat)
} else {
  txi <- readRDS(args$counts)
  y <- DGEList(counts = txi$counts)
  # Apply tximport length offset
  norm_mat <- txi$length
  norm_mat <- norm_mat / exp(rowMeans(log(norm_mat)))
  norm_counts <- txi$counts / norm_mat
  eff_lib <- calcNormFactors(norm_counts) * colSums(norm_counts)
  norm_mat <- sweep(norm_mat, 2, eff_lib, "*")
  y <- scaleOffset(y, log(norm_mat))
}

design <- model.matrix(formula, data = samples)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

# Build contrast vector from named levels
contrast_factor <- args$contrast[1]
contrast_test <- args$contrast[2]
contrast_ref <- args$contrast[3]
coef_name <- paste0(contrast_factor, contrast_test)
qlf <- glmQLFTest(fit, coef = coef_name)

res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  rename(log2FoldChange = logFC, padj = FDR, pvalue = PValue) %>%
  arrange(padj)

write_tsv(res_df, args$out_tsv)

# Volcano
pdf(args$out_volcano)
plot(res_df$log2FoldChange, -log10(res_df$padj),
     pch = 20, cex = 0.5,
     xlab = "log2 Fold Change", ylab = "-log10(FDR)",
     main = paste("edgeR:", contrast_test, "vs", contrast_ref))
abline(h = -log10(0.05), col = "red", lty = 2)
dev.off()

# MA
pdf(args$out_ma)
plotMD(qlf, main = paste("edgeR:", contrast_test, "vs", contrast_ref))
abline(h = 0, col = "blue")
dev.off()
```

- [ ] **Step 2: Add edgeR rule (same pattern as DESeq2)**

```python
rule rnaseq_edger:
    message: "Differential expression with edgeR"
    conda: ENV_R
    input:
        counts     = lambda wc: get_count_matrix(wc, de_map),
        sample_tsv = SAMPLE_TSV,
    log:
        cmd = f"{D_LOGS}/{{experiment}}_{{align_method}}_rnaseq_edger.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_{{align_method}}_rnaseq_edger.tsv"
    params:
        script     = f"{R_RNASEQ}/scripts/rnaseq_edger.R",
        input_type = lambda wc: "featurecounts" if wc.align_method == "star" else "tximport",
        design     = lambda wc: de_map[wc.experiment]["design"],
        contrast   = lambda wc: " ".join(de_map[wc.experiment]["contrast"]),
    threads: 1
    output:
        tsv     = f"{D_RNASEQ}/de/{{experiment}}.{{align_method}}.edger_results.tsv",
        volcano = f"{D_RNASEQ}/de/{{experiment}}.{{align_method}}.edger_volcano.pdf",
        ma      = f"{D_RNASEQ}/de/{{experiment}}.{{align_method}}.edger_ma.pdf",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[edgeR] $(date) experiment={wildcards.experiment} method={wildcards.align_method}"

        Rscript {params.script} \
          --counts "{input.counts}" \
          --input-type "{params.input_type}" \
          --sample-tsv "{input.sample_tsv}" \
          --design "{params.design}" \
          --contrast {params.contrast} \
          --out-tsv "{output.tsv}" \
          --out-volcano "{output.volcano}" \
          --out-ma "{output.ma}"
        """
```

- [ ] **Step 3: Commit**

```bash
git add scripts/rnaseq_edger.R workflows/rnaseq.smk
git commit -m "Add edgeR R script and Snakemake rule"
```

### Task 14: Write scripts/rnaseq_limma.R and rule

**Files:**
- Create: `scripts/rnaseq_limma.R`
- Modify: `workflows/rnaseq.smk`

- [ ] **Step 1: Write limma-voom R script**

```r
#!/usr/bin/env Rscript
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)

library(argparse)
library(tidyverse)
library(edgeR)
library(limma)

parser <- ArgumentParser(description = "Differential expression with limma-voom")
parser$add_argument("--counts", required = TRUE)
parser$add_argument("--input-type", required = TRUE, choices = c("featurecounts", "tximport"))
parser$add_argument("--sample-tsv", required = TRUE)
parser$add_argument("--design", required = TRUE)
parser$add_argument("--contrast", nargs = 3, required = TRUE)
parser$add_argument("--out-tsv", required = TRUE)
parser$add_argument("--out-volcano", required = TRUE)
parser$add_argument("--out-ma", required = TRUE)
args <- parser$parse_args()

samples <- read_tsv(args$sample_tsv, show_col_types = FALSE)
formula <- as.formula(args$design)

if (args$input_type == "featurecounts") {
  cts <- read.delim(args$counts)
  count_mat <- as.matrix(cts[, -1])
  rownames(count_mat) <- cts[[1]]
  y <- DGEList(counts = count_mat)
} else {
  txi <- readRDS(args$counts)
  y <- DGEList(counts = txi$counts)
  norm_mat <- txi$length
  norm_mat <- norm_mat / exp(rowMeans(log(norm_mat)))
  norm_counts <- txi$counts / norm_mat
  eff_lib <- calcNormFactors(norm_counts) * colSums(norm_counts)
  norm_mat <- sweep(norm_mat, 2, eff_lib, "*")
  y <- scaleOffset(y, log(norm_mat))
}

design <- model.matrix(formula, data = samples)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

v <- voom(y, design)
fit <- lmFit(v, design)

# Build contrast
contrast_factor <- args$contrast[1]
contrast_test <- args$contrast[2]
contrast_ref <- args$contrast[3]
coef_test <- paste0(contrast_factor, contrast_test)
coef_ref <- paste0(contrast_factor, contrast_ref)

# If reference level is the intercept, coef_ref won't exist — use coef_test directly
if (coef_test %in% colnames(design)) {
  contrast_vec <- makeContrasts(contrasts = coef_test, levels = design)
  fit2 <- contrasts.fit(fit, contrast_vec)
} else {
  fit2 <- fit
}

fit2 <- eBayes(fit2)
res <- topTable(fit2, number = Inf, sort.by = "P")

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  rename(log2FoldChange = logFC, padj = adj.P.Val, pvalue = P.Value) %>%
  arrange(padj)

write_tsv(res_df, args$out_tsv)

# Volcano
pdf(args$out_volcano)
plot(res_df$log2FoldChange, -log10(res_df$padj),
     pch = 20, cex = 0.5,
     xlab = "log2 Fold Change", ylab = "-log10(adj.P.Val)",
     main = paste("limma:", contrast_test, "vs", contrast_ref))
abline(h = -log10(0.05), col = "red", lty = 2)
dev.off()

# MA
pdf(args$out_ma)
plotMD(fit2, main = paste("limma:", contrast_test, "vs", contrast_ref))
abline(h = 0, col = "blue")
dev.off()
```

- [ ] **Step 2: Add limma rule (same pattern)**

```python
rule rnaseq_limma:
    message: "Differential expression with limma-voom"
    conda: ENV_R
    input:
        counts     = lambda wc: get_count_matrix(wc, de_map),
        sample_tsv = SAMPLE_TSV,
    log:
        cmd = f"{D_LOGS}/{{experiment}}_{{align_method}}_rnaseq_limma.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_{{align_method}}_rnaseq_limma.tsv"
    params:
        script     = f"{R_RNASEQ}/scripts/rnaseq_limma.R",
        input_type = lambda wc: "featurecounts" if wc.align_method == "star" else "tximport",
        design     = lambda wc: de_map[wc.experiment]["design"],
        contrast   = lambda wc: " ".join(de_map[wc.experiment]["contrast"]),
    threads: 1
    output:
        tsv     = f"{D_RNASEQ}/de/{{experiment}}.{{align_method}}.limma_results.tsv",
        volcano = f"{D_RNASEQ}/de/{{experiment}}.{{align_method}}.limma_volcano.pdf",
        ma      = f"{D_RNASEQ}/de/{{experiment}}.{{align_method}}.limma_ma.pdf",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[limma] $(date) experiment={wildcards.experiment} method={wildcards.align_method}"

        Rscript {params.script} \
          --counts "{input.counts}" \
          --input-type "{params.input_type}" \
          --sample-tsv "{input.sample_tsv}" \
          --design "{params.design}" \
          --contrast {params.contrast} \
          --out-tsv "{output.tsv}" \
          --out-volcano "{output.volcano}" \
          --out-ma "{output.ma}"
        """
```

- [ ] **Step 3: Commit**

```bash
git add scripts/rnaseq_limma.R workflows/rnaseq.smk
git commit -m "Add limma-voom R script and Snakemake rule"
```

## Chunk 6: Functional Enrichment + Test Wrapper

### Task 15: Write scripts/rnaseq_enrichment.R and rule

**Files:**
- Create: `scripts/rnaseq_enrichment.R`
- Modify: `workflows/rnaseq.smk`

- [ ] **Step 1: Write enrichment R script**

```r
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
```

- [ ] **Step 2: Add enrichment rule to rnaseq.smk**

```python
rule rnaseq_enrichment:
    message: "Functional enrichment analysis (GO/KEGG)"
    conda: ENV_R
    input:
        de_results = f"{D_RNASEQ}/de/{{experiment}}.{{align_method}}.{{de_tool}}_results.tsv",
    log:
        cmd = f"{D_LOGS}/{{experiment}}_{{align_method}}_{{de_tool}}_rnaseq_enrichment.log",
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_{{align_method}}_{{de_tool}}_rnaseq_enrichment.tsv"
    params:
        script = f"{R_RNASEQ}/scripts/rnaseq_enrichment.R",
        orgdb  = lambda wc: de_map[wc.experiment]["orgdb"],
    threads: 1
    output:
        go_tsv    = f"{D_RNASEQ}/enrichment/{{experiment}}.{{align_method}}.{{de_tool}}.go.tsv",
        go_plot   = f"{D_RNASEQ}/enrichment/{{experiment}}.{{align_method}}.{{de_tool}}.go.pdf",
        kegg_tsv  = f"{D_RNASEQ}/enrichment/{{experiment}}.{{align_method}}.{{de_tool}}.kegg.tsv",
        kegg_plot = f"{D_RNASEQ}/enrichment/{{experiment}}.{{align_method}}.{{de_tool}}.kegg.pdf",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[enrichment] $(date) experiment={wildcards.experiment} method={wildcards.align_method} tool={wildcards.de_tool}"

        mkdir -p "$(dirname "{output.go_tsv}")"

        Rscript {params.script} \
          --de-results "{input.de_results}" \
          --orgdb "{params.orgdb}" \
          --out-go-tsv "{output.go_tsv}" \
          --out-go-plot "{output.go_plot}" \
          --out-kegg-tsv "{output.kegg_tsv}" \
          --out-kegg-plot "{output.kegg_plot}"
        """
```

- [ ] **Step 3: Commit**

```bash
git add scripts/rnaseq_enrichment.R workflows/rnaseq.smk
git commit -m "Add functional enrichment R script and Snakemake rule"
```

### Task 16: Write workflows/test.smk

**Files:**
- Create: `workflows/test.smk`

Follow the frag/emseq test.smk pattern: imports, config loading, path resolution, constants, SampleTable class, de_map, rule all, symlinks, include.

- [ ] **Step 1: Write test.smk**

```python
# workflows/test.smk
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)
# ============================================================
# RNA-SEQ FULL PIPELINE TEST WRAPPER SNAKEFILE

import os
import pandas as pd

# --- Config loading & path expansion ---
configfile: "config/test.yaml"

def resolve_config_paths(config_dict):
    for k, v in config_dict.items():
        if isinstance(v, str):
            config_dict[k] = os.path.expandvars(os.path.expanduser(v))
        elif isinstance(v, dict):
            resolve_config_paths(v)
        elif isinstance(v, list):
            config_dict[k] = [os.path.expandvars(os.path.expanduser(i))
                              if isinstance(i, str) else i for i in v]

resolve_config_paths(config)

# --- Constants ---
D_DATA      = config["main-data-dir"]
D_RNASEQ    = f"{D_DATA}/rnaseq"
D_LOGS      = f"{D_DATA}/logs"
D_BENCHMARK = f"{D_DATA}/benchmark"
D_INPUTS    = f"{D_DATA}/inputs"

ENV_ALIGN = config["envs"]["align"]
ENV_QC    = config["envs"]["qc"]
ENV_R     = config["envs"]["r"]

R_RNASEQ  = config["repos"]["rnaseq"]
SAMPLE_TSV = config["sample-tsv-path"]

# --- Sample table ---
class SampleTable:
    def __init__(self, tsv_path, selected_ids):
        df = pd.read_csv(tsv_path, sep="\t")
        missing = sorted(set(selected_ids) - set(df["library_id"]))
        if missing:
            raise ValueError(f"library_id not found in TSV: {missing}")
        self.df = df[df["library_id"].isin(selected_ids)].copy()

    @property
    def library_ids(self):
        return sorted(self.df["library_id"].unique())

    @property
    def r1_map(self):
        return dict(zip(self.df["library_id"], self.df["r1_basename"]))

    @property
    def r2_map(self):
        return dict(zip(self.df["library_id"], self.df["r2_basename"]))

# Collect all library IDs from all de-map experiments
de_map = config["de-map"]
all_libs = sorted(set(lib for exp in de_map.values() for lib in exp["libs"]))
samples = SampleTable(tsv_path=SAMPLE_TSV, selected_ids=all_libs)
RNASEQ_LIBRARY_IDS = samples.library_ids

# --- Helper: count matrix path by alignment method ---
def get_count_matrix(wc, de_map):
    if wc.align_method == "star":
        return f"{D_RNASEQ}/counts/{wc.experiment}.star.featurecounts.tsv"
    elif wc.align_method == "salmon":
        return f"{D_RNASEQ}/counts/{wc.experiment}.salmon.txi.rds"

# --- Rule all ---
rule all:
    input:
        # FASTQ QC (raw + trimmed)
        expand(
            f"{D_RNASEQ}/fastqs/{{library_id}}.trimmed_{{read}}.fastq.gz",
            library_id=RNASEQ_LIBRARY_IDS,
            read=["R1", "R2"],
        ),
        expand(
            f"{D_RNASEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.zip",
            library_id=RNASEQ_LIBRARY_IDS,
            processing=["raw", "trimmed"],
            read=["R1", "R2"],
        ),
        # EDA (per experiment × align_method)
        [f"{D_RNASEQ}/eda/{exp}.{am}.pca.pdf"
         for exp in de_map
         for am in de_map[exp]["align_method"]],
        # DE results (per experiment × align_method × tool)
        [f"{D_RNASEQ}/de/{exp}.{am}.{tool}_results.tsv"
         for exp in de_map
         for am in de_map[exp]["align_method"]
         for tool in de_map[exp]["tools"]],
        # Enrichment (per experiment × align_method × tool)
        [f"{D_RNASEQ}/enrichment/{exp}.{am}.{tool}.go.tsv"
         for exp in de_map
         for am in de_map[exp]["align_method"]
         for tool in de_map[exp]["tools"]],
        # BAM QC (STAR path only)
        [f"{D_RNASEQ}/qc/rseqc/{lib}.{de_map[exp]['ref_name']}.read_distribution.txt"
         for exp in de_map
         if "star" in de_map[exp]["align_method"]
         for lib in de_map[exp]["libs"]],
        [f"{D_RNASEQ}/qc/qualimap/{lib}.{de_map[exp]['ref_name']}/qualimapReport.html"
         for exp in de_map
         if "star" in de_map[exp]["align_method"]
         for lib in de_map[exp]["libs"]],
        # MultiQC
        f"{D_RNASEQ}/qc/multiqc.html",

# --- Symlink input FASTQs ---
rule symlink_input_fastqs:
    message: "Create symlinks for raw input FASTQs"
    input:
        r1 = lambda wc: f"{D_INPUTS}/{samples.r1_map[wc.library_id]}",
        r2 = lambda wc: f"{D_INPUTS}/{samples.r2_map[wc.library_id]}",
    output:
        r1 = f"{D_RNASEQ}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_RNASEQ}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    shell:
        """
        mkdir -p "$(dirname "{output.r1}")"
        ln -sfr "{input.r1}" "{output.r1}"
        ln -sfr "{input.r2}" "{output.r2}"
        """

# --- Include modular pipeline ---
include: "rnaseq.smk"
```

- [ ] **Step 2: Dry-run to validate DAG**

```bash
conda run -n basecamp snakemake -s workflows/test.smk --dry-run 2>&1 | head -50
```

Expected: Either successful dry-run showing rule execution order, or missing input file errors (since test data isn't committed yet). No syntax/DAG errors.

- [ ] **Step 3: Commit**

```bash
git add workflows/test.smk
git commit -m "Add test wrapper snakefile with rule all, SampleTable, and de_map dispatch"
```

## Chunk 7: CI/CD + CLAUDE.md + Final Cleanup

### Task 17: Write GitHub Actions workflows

**Files:**
- Create: `.github/workflows/generate-data.yaml`
- Create: `.github/workflows/test.yaml`

- [ ] **Step 1: Write generate-data.yaml**

```yaml
name: Generate Test Data
on:
  workflow_dispatch:
jobs:
  generate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: conda-incubator/setup-miniconda@v3
        with:
          auto-activate-base: true
      - name: Install SRA toolkit
        run: conda install -c bioconda sra-tools samtools
      - name: Run get_test_data.sh
        run: bash tools/get_test_data.sh
      - name: Validate outputs
        run: |
          test -f tests/full/inputs/mm10_chr19.fa.gz
          test -f tests/full/inputs/mm10_chr19.gtf.gz
          test -f tests/full/inputs/hg38_chr22.fa.gz
          test -f tests/full/inputs/hg38_chr22.gtf.gz
          for i in $(seq 1 8); do
            lib=$(printf "lib%03d" $i)
            test -f "tests/full/inputs/${lib}_R1.fastq.gz"
            test -f "tests/full/inputs/${lib}_R2.fastq.gz"
          done
          echo "All test data files validated"
```

- [ ] **Step 2: Write test.yaml**

```yaml
name: Pipeline Test
on:
  push:
    branches: [master]
  pull_request:
    branches: [master]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: conda-incubator/setup-miniconda@v3
        with:
          auto-activate-base: true
      - name: Install snakemake
        run: conda install -c conda-forge -c bioconda snakemake
      - name: Dry run
        run: snakemake -s workflows/test.smk --dry-run
      - name: Full pipeline test
        run: snakemake -s workflows/test.smk --use-conda --cores 4
```

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/generate-data.yaml .github/workflows/test.yaml
git commit -m "Add GitHub Actions CI/CD workflows"
```

### Task 18: Write CLAUDE.md

**Files:**
- Create: `CLAUDE.md`

- [ ] **Step 1: Write CLAUDE.md**

```markdown
# RNA-seq Biopipe Module

## Project Structure
- `rna-seq.org` — literate programming source (tangle to produce all code files)
- `workflows/rnaseq.smk` — modular pipeline (no rule all), included by wrappers
- `workflows/test.smk` — test wrapper with config, rule all, symlinks
- `config/test.yaml` — test configuration with de-map dispatch pattern

## Conventions
- All Snakemake rules prefixed `rnaseq_`
- All R scripts use argparse for CLI args
- Conda environments: rnaseq-align (tools), rnaseq-qc (BAM QC), rnaseq-r (R/Bioconductor)
- Test data committed to `tests/full/inputs/`, all other test outputs git-ignored
- Follows emseq/frag biopipe patterns — see those repos for reference

## Running
- Dry run: `conda run -n basecamp snakemake -s workflows/test.smk --dry-run`
- Full test: `conda run -n basecamp snakemake -s workflows/test.smk --use-conda --cores 4`

## de-map Pattern
The `de-map` config key dispatches differential expression experiments (like emseq's `meth-map`).
Each experiment specifies: libs, ref_name, align_method, design, contrast, tools, orgdb.
```

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "Add CLAUDE.md with project conventions"
```

### Task 19: Run get_test_data.sh and commit test inputs

**Note:** This task requires network access and SRA toolkit. The actual GEO accessions must be chosen at implementation time.

- [ ] **Step 1: Run get_test_data.sh**

```bash
bash tools/get_test_data.sh
```

- [ ] **Step 2: Verify outputs**

```bash
ls -la tests/full/inputs/
# Should see: mm10_chr19.fa.gz, mm10_chr19.gtf.gz, hg38_chr22.fa.gz, hg38_chr22.gtf.gz
# Plus 16 FASTQ files (8 libraries × R1/R2)
```

- [ ] **Step 3: Commit test data**

```bash
git add tests/full/inputs/ data/test-samples.tsv
git commit -m "Add committed test data (mouse chr19 + human chr22 from GEO)"
```

### Task 20: Final dry-run validation

- [ ] **Step 1: Run snakemake dry-run**

```bash
conda run -n basecamp snakemake -s workflows/test.smk --dry-run
```

Expected: Successful dry-run showing all 16+ rules in execution order, no errors.

- [ ] **Step 2: If errors, fix and re-commit**

Common issues: missing variables in rnaseq.smk (they're defined in test.smk and propagated via include), wildcard mismatches, lambda input resolution errors.

- [ ] **Step 3: Final commit**

```bash
git add -A
git commit -m "Pipeline dry-run validated — RNA-seq biopipe module v1.0.0"
```
