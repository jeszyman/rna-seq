# RNA-seq Biopipe Module

## Project Structure
- `rna-seq.org` -- literate programming source (tangle to produce all code files)
- `workflows/rnaseq.smk` -- modular pipeline (no rule all), included by wrappers
- `workflows/test.smk` -- test wrapper with config, rule all, symlinks
- `config/test.yaml` -- test configuration with de-map dispatch pattern

## Conventions
- All Snakemake rules prefixed `rnaseq_`
- All R scripts use argparse for CLI args
- Conda environments: rnaseq-align (tools), rnaseq-qc (BAM QC), rnaseq-r (R/Bioconductor)
- Test data committed to `tests/full/inputs/`, all other test outputs git-ignored
- Follows emseq/frag biopipe patterns -- see those repos for reference

## Running
- Dry run: `conda run -n basecamp snakemake -s workflows/test.smk --dry-run`
- Full test: `conda run -n basecamp snakemake -s workflows/test.smk --use-conda --cores 4`

## de-map Pattern
The `de-map` config key dispatches differential expression experiments (like emseq's `meth-map`).
Each experiment specifies: libs, ref_name, align_method, design, contrast, tools, orgdb.
