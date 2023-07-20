rule pe_rna_seq_fastp:
    input:
        read1 = f"{rna_dir}/fastqs/pe/{{library}}_raw_R1.fastq.gz",
        read2 = f"{rna_dir}/fastqs/pe/{{library}}_raw_R2.fastq.gz",
    log: html = f"{log_dir}/{{library}}_pe_rna_seq_fastp.html",
    output:
        read1 = f"{rna_dir}/fastqs/pe/{{library}}_proc_R1.fastq.gz",
        read2 = f"{rna_dir}/fastqs/pe/{{library}}_proc_R2.fastq.gz",
        failed = f"{rna_dir}/fastqs/pe/{{library}}_failed_fastp.fastq.gz",
        unpaired1 = f"{rna_dir}/fastqs/pe/{{library}}_unpaired_R1.fastq.gz",
        unpaired2 = f"{rna_dir}/fastqs/pe/{{library}}_unpaired_R2.fastq.gz",
        json = f"{qc_dir}/{{library}}_fastp.json",
        cmd = f"{qc_dir}/{{library}}_fastp.log",
    params:
        script = f"{rna_script_dir}/pe_rna_seq_fastp.sh",
        threads = 4
    resources:
        mem_mb = 500
    shell:
        """
        {params.script} \
        {input.read1} \
        {input.read2} \
        {log.html} \
        {output.json} \
        {output.read1} \
        {output.read2} \
        {output.failed} \
        {output.unpaired1} \
        {output.unpaired2} \
        {params.threads} &> {output.cmd}
        """

rule pe_rna_seq_fastqc:
    input: f"{rna_dir}/fastqs/pe/{{library}}_{{processing}}_{{read}}.fastq.gz",
    log: f"{log_dir}/{{library}}_{{processing}}_{{read}}_rna_seq_fastqc.log",
    output: f"{qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.zip",
    params:
        out_dir = qc_dir,
        script = f"{rna_script_dir}/rna_seq_fastqc.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.out_dir} {params.threads} &> {log}
        """

rule pe_quant_with_salmon:
    input:
        index = f"{ref_dir}/{{build}}_salmon",
        read1 = f"{rna_dir}/fastqs/pe/{{library}}_proc_R1.fastq.gz",
        read2 = f"{rna_dir}/fastqs/pe/{{library}}_proc_R2.fastq.gz",
    log: f"{log_dir}/{{library}}_{{build}}_pe_quant_with_salmon.log",
    output: f"{rna_dir}/salmon/{{library}}_{{build}}/quant.sf",
    params:
        out_dir = f"{rna_dir}/salmon/{{library}}_{{build}}",
        script = f"{rna_script_dir}/pe_quant_with_salmon.sh",
        threads = 4,
    shell:
        """
        {params.script} \
        {input.index} \
        {input.read1} \
        {input.read2} \
        {params.out_dir} \
        {params.threads} > {log} 2>&1 &&
        [[ -s {output[0]} ]] || (echo "Output file is empty: {output[0]}" && exit 1)
        """

rule make_ensembl_de_gtf:
    conda: "rna",
    input:  f"{ref_dir}/{{build}}.gtf.gz",
    log:    f"{log_dir}/{{build}}_make_ensembl_de_gtf.log",
    output: f"{ref_dir}/{{build}}_protein_coding.gtf",
    shell:
        """
        zcat {input} | grep "protein_coding" > {output} 2> {log}
        """

rule make_txdb_from_gtf:
    input: f"{ref_dir}/{{build}}_protein_coding.gtf",
    log: f"{log_dir}/{{build}}_make_txdb_from_gtf.log",
    output: f"{ref_dir}/{{build}}_protein.txdb",
    params: script = f"{rna_script_dir}/make_txdb_from_gtf.R",
    shell:
        """
        Rscript {params.script} {input} {wildcards.build} > {log} 2>&1
        cp /tmp/{wildcards.build}_protein.txdb {output}
        """

rule make_dge_design:
    input:
        libraries_full = libraries_full_rds,
    log: f"{log_dir}/{{experiment}}_make_dge_design.log",
    output: f"{rna_dir}/{{experiment}}/design.rds",
    params:
        formula = lambda wildcards: rna_map[wildcards.experiment]['formula'],
        libs = lambda wildcards: rna_map[wildcards.experiment]['libs'],
        script = f"{rna_script_dir}/make_dge_design.R",
    shell:
        """
        Rscript {params.script} \
        {input.libraries_full} \
        "{params.formula}" \
        "{params.libs}" \
        {output} \
        > {log} 2>&1
        """

rule make_salmon_txi:
    input:
        salmon = lambda wildcards: expand(f"{rna_dir}/salmon/{{library}}_{{build}}/quant.sf",
                                          library = rna_map[wildcards.experiment]['libs'],
                                          build = rna_map[wildcards.experiment]['build']),
        txdb = lambda wildcards: rna_map[wildcards.experiment]['txdb'],
    log: f"{log_dir}/{{experiment}}_make_salmon_txi.log",
    output: f"{rna_dir}/{{experiment}}/txi.rds",
    params:
        script = rna_script_dir + "/make_salmon_txi.R",
    shell:
        """
        Rscript {params.script} \
        "{input.salmon}" \
        {input.txdb} \
        {output} > {log} 2>&1
        """

rule norm_txi:
    input:
        design = f"{rna_dir}/{{experiment}}/design.rds",
        txi = f"{rna_dir}/{{experiment}}/txi.rds",
    log: f"{log_dir}/{{experiment}}_norm_txi.log",
    output:
        dge = f"{rna_dir}/{{experiment}}/dge.rds",
        glm = f"{rna_dir}/{{experiment}}/fit.rds",
        cpm = f"{rna_dir}/{{experiment}}/cpm.tsv",
    params: script = f"{rna_script_dir}/norm_txi.R",
    shell:
        """
        Rscript {params.script} \
        {input.design} \
        {input.txi} \
        {output.dge} \
        {output.glm} \
        {output.cpm} \
        > {log} 2>&1
        """

rule make_cpm_pca:
    input:
        cpm = f"{rna_dir}/{{experiment}}/cpm.tsv",
        libraries_full = libraries_full_rds,
    log: f"{log_dir}/{{experiment}}_make_cpm_pca.log",
    output:
        f"{rna_dir}/{{experiment}}/pca.png",
        f"{rna_dir}/{{experiment}}/pca.svg",
    params:
        formula = lambda wildcards: rna_map[wildcards.experiment]['formula'],
        script = f"{rna_script_dir}/make_cpm_pca.R",
    shell:
        """
        Rscript {params.script} \
        {input.cpm} \
        "{params.formula}" \
        {input.libraries_full} \
        {output} > {log} 2>&1
        """

rule make_edger_contrast_de:
    input:
        design = lambda wildcards: dge_map[wildcards.contrast]['design'],
        fit = lambda wildcards: dge_map[wildcards.contrast]['fit'],
    log: f"{log_dir}/{{contrast}}_make_edger_contrast_de.log",
    output: f"{rna_dir}/dge/{{contrast}}.tsv",
    params:
        cohorts_str = lambda wildcards: dge_map[wildcards.contrast]['cohorts_str'],
        script = f"{rna_script_dir}/make_edger_contrast_de.R",
    shell:
        """
        Rscript {params.script} \
        "{params.cohorts_str}" \
        {input.design} \
        {input.fit} \
        {output} > {log} 2>&1
        """
