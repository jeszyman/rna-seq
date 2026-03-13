# workflows/test.smk
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)
# ============================================================
# RNA-SEQ FULL PIPELINE TEST WRAPPER SNAKEFILE

import os
import pandas as pd

# ------------------------------------------------------------------------------
# Load YAML Configuration
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

D_DATA      = config["main-data-dir"]
D_RNASEQ    = f"{D_DATA}/rnaseq"
D_LOGS      = f"{D_DATA}/logs"
D_BENCHMARK = f"{D_DATA}/benchmark"
D_INPUTS    = f"{D_DATA}/inputs"

ENV_ALIGN = config["envs"]["align"]
ENV_QC    = config["envs"]["qc"]
ENV_R     = config["envs"]["r"]

R_RNASEQ   = config["repos"]["rnaseq"]
SAMPLE_TSV = config["sample-tsv-path"]

# ------------------------------------------------------------------------------
# Sample Table
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Rule all
# ------------------------------------------------------------------------------

rule all:
    input:
        # FASTQ QC (trimmed outputs)
        expand(
            f"{D_RNASEQ}/fastqs/{{library_id}}.trimmed_{{read}}.fastq.gz",
            library_id=RNASEQ_LIBRARY_IDS,
            read=["R1", "R2"],
        ),
        # FastQC (raw + trimmed)
        expand(
            f"{D_RNASEQ}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.zip",
            library_id=RNASEQ_LIBRARY_IDS,
            processing=["raw", "trimmed"],
            read=["R1", "R2"],
        ),
        # EDA (per experiment x align_method)
        [f"{D_RNASEQ}/eda/{exp}.{am}.pca.pdf"
         for exp in de_map
         for am in de_map[exp]["align_method"]],
        # DE results (per experiment x align_method x tool)
        [f"{D_RNASEQ}/de/{exp}.{am}.{tool}_results.tsv"
         for exp in de_map
         for am in de_map[exp]["align_method"]
         for tool in de_map[exp]["tools"]],
        # Enrichment (per experiment x align_method x tool)
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

# ------------------------------------------------------------------------------
# Symlink input FASTQs
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Include modular pipeline
# ------------------------------------------------------------------------------

include: "rnaseq.smk"
