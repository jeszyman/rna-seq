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

# ==============================================================================
# FASTQ QC
# ==============================================================================

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

# ==============================================================================
# STAR Alignment
# ==============================================================================

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
        genome_dir     = f"{D_RNASEQ}/ref/star/{{ref_name}}",
        sa_index_bases = lambda wc: config["refs"][wc.ref_name].get("star_sa_bases", 14),
    threads: 8
    output:
        sa = f"{D_RNASEQ}/ref/star/{{ref_name}}/SA",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[STAR index] $(date) ref={wildcards.ref_name} threads={threads}"

        mkdir -p "{params.genome_dir}"
        tmpdir=$(mktemp -d)

        # Decompress if gzipped
        if file -b "{input.fasta}" | grep -qi gzip; then
            zcat "{input.fasta}" > "$tmpdir/genome.fa"
        else
            cp "{input.fasta}" "$tmpdir/genome.fa"
        fi

        if file -b "{input.gtf}" | grep -qi gzip; then
            zcat "{input.gtf}" > "$tmpdir/genes.gtf"
        else
            cp "{input.gtf}" "$tmpdir/genes.gtf"
        fi

        STAR --runMode genomeGenerate \
          --genomeDir "{params.genome_dir}" \
          --genomeFastaFiles "$tmpdir/genome.fa" \
          --sjdbGTFfile "$tmpdir/genes.gtf" \
          --genomeSAindexNbases {params.sa_index_bases} \
          --runThreadN {threads}

        rm -rf "$tmpdir"
        """

rule rnaseq_star_align:
    message: "STAR alignment of RNA-seq reads"
    conda: ENV_ALIGN
    input:
        r1       = f"{D_RNASEQ}/fastqs/{{library_id}}.trimmed_R1.fastq.gz",
        r2       = f"{D_RNASEQ}/fastqs/{{library_id}}.trimmed_R2.fastq.gz",
        star_idx = f"{D_RNASEQ}/ref/star/{{ref_name}}/SA",
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

# ==============================================================================
# Salmon Alignment
# ==============================================================================

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

        # Decompress inputs if gzipped
        if file -b "{input.fasta}" | grep -qi gzip; then
            zcat "{input.fasta}" > "$tmpdir/genome.fa"
        else
            cp "{input.fasta}" "$tmpdir/genome.fa"
        fi

        if file -b "{input.gtf}" | grep -qi gzip; then
            zcat "{input.gtf}" > "$tmpdir/genes.gtf"
        else
            cp "{input.gtf}" "$tmpdir/genes.gtf"
        fi

        # Extract transcriptome from GTF + genome FASTA using gffread
        # For gentrome: concatenate transcriptome + genome, use genome as decoy
        grep "^>" "$tmpdir/genome.fa" \
          | cut -d " " -f 1 | sed 's/>//' > "$tmpdir/decoys.txt"

        # Build gentrome (transcripts + genome)
        gffread -w "$tmpdir/transcripts.fa" -g "$tmpdir/genome.fa" "$tmpdir/genes.gtf"
        cat "$tmpdir/transcripts.fa" "$tmpdir/genome.fa" \
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

# ==============================================================================
# tximport (Salmon → gene-level counts)
# ==============================================================================

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

# ==============================================================================
# Merge featureCounts (per-library → experiment matrix)
# ==============================================================================

rule rnaseq_merge_featurecounts:
    message: "Merge per-library featureCounts into experiment count matrix"
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

# ==============================================================================
# Exploratory Data Analysis
# ==============================================================================

def get_count_matrix(wc, de_map):
    """Resolve count matrix path based on alignment method."""
    if wc.align_method == "star":
        return f"{D_RNASEQ}/counts/{wc.experiment}.star.featurecounts.tsv"
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

# ==============================================================================
# BAM QC (STAR path only)
# ==============================================================================

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

# ==============================================================================
# MultiQC
# ==============================================================================

rule rnaseq_multiqc:
    message: "MultiQC aggregation of QC reports"
    conda: ENV_ALIGN
    input:
        # MultiQC does NOT prompt inputs to run — list explicit trigger files
        fastp = expand(
            f"{D_RNASEQ}/qc/{{library_id}}_fastp.json",
            library_id=RNASEQ_LIBRARY_IDS,
        ),
    log:
        cmd = f"{D_LOGS}/rnaseq_multiqc.log",
    benchmark:
        f"{D_BENCHMARK}/rnaseq_multiqc.tsv"
    threads: 1
    output:
        html = f"{D_RNASEQ}/qc/multiqc.html",
    params:
        qc_dir  = f"{D_RNASEQ}/qc",
        out_dir = f"{D_RNASEQ}/qc",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[MultiQC] $(date)"

        multiqc "{params.qc_dir}" \
          --outdir "{params.out_dir}" \
          --filename multiqc.html \
          --force
        """

# ==============================================================================
# Differential Expression
# ==============================================================================

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

# ==============================================================================
# Functional Enrichment
# ==============================================================================

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
