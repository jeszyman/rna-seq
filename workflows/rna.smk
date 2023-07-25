rule make_wtrans_filtered_gtf:
    input: f"{ref_dir}/{{build}}.gtf.gz",
    log: f"{log_dir}/{{build}}_make_wtrans_filtered_gtf.log",
    output: f"{ref_dir}/{{build}}_wtrans.gtf.gz",
    params: script = f"{rna_script_dir}/make_wtrans_filtered_gtf.sh",
    shell:
        """
        {params.script} {input} {output} > {log} 2>&1
        """

rule make_annotation_from_gtf:
    input: f"{ref_dir}/{{build}}_wtrans.gtf.gz",
    log: f"{log_dir}/{{build}}_make_annotation_from_gtf.log",
    output: f"{ref_dir}/{{build}}_wtrans_annotation.tsv",
    params:
        bmart_data =  lambda wildcards: build_map[wildcards.build]['bmart_data'],
        script = f"{rna_script_dir}/make_annotation_from_gtf.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {params.bmart_data} \
        {output} \
        > {log} 2>&1
        """

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

rule make_edger_lrt:
    input:
        design = lambda wildcards: lrt_map[wildcards.contrast]['design'],
        fit = lambda wildcards: lrt_map[wildcards.contrast]['fit'],
    log: f"{log_dir}/{{contrast}}_make_edger_lrt.log",
    output: f"{rna_dir}/contrasts/lrt_{{contrast}}/lrt_{{contrast}}.tsv",
    params: script = f"{rna_script_dir}/make_edger_lrt.R",
    shell:
        """
        Rscript {params.script} {input} {output} > {log} 2>&1
        """

rule make_edger_contrast_de:
    input:
        design = lambda wildcards: dge_map[wildcards.contrast]['design'],
        fit = lambda wildcards: dge_map[wildcards.contrast]['fit'],
    log: f"{log_dir}/{{contrast}}_make_edger_contrast_de.log",
    output: f"{rna_dir}/contrasts/{{contrast}}/{{contrast}}.tsv",
    params:
        annotation_tsv = lambda wildcards: dge_map[wildcards.contrast]['annotation_tsv'],
        cohorts_str = lambda wildcards: dge_map[wildcards.contrast]['cohorts_str'],
        script = f"{rna_script_dir}/make_edger_contrast_de.R",
    shell:
        """
        Rscript {params.script} \
        {input.design} \
        {input.fit} \
        {params.annotation_tsv} \
        "{params.cohorts_str}" \
        {output} > {log} 2>&1
        """

rule rna_volcano:
    input: f"{rna_dir}/contrasts/{{contrast}}/{{contrast}}.tsv",
    log: f"{log_dir}/{{contrast}}_rna_volcano.log",
    output: f"{rna_dir}/contrasts/{{contrast}}/{{contrast}}_volcano.pdf",
    params: script = f"{rna_script_dir}/rna_volcano.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        > {log} 2>&1
        """

rule gsea_from_edger:
    input: f"{rna_dir}/contrasts/{{contrast}}/{{contrast}}.tsv",
    log: f"{log_dir}/{{contrast}}_{{pathset}}_gsea_from_edger.log",
    output:
       f"{rna_dir}/contrasts/{{contrast}}/{{contrast}}_gsea_{{pathset}}.tsv",
       f"{rna_dir}/contrasts/{{contrast}}/{{contrast}}_gsea_{{pathset}}.xlsx",
    params: script = f"{rna_script_dir}/gsea_from_edger.R",
    shell:
        """
        Rscript {params.script} {input} {wildcards.pathset} {output} > {log} 2>&1
        """

rule make_dge_design:
    input:
        libraries_full = libraries_full_rds,
    log: f"{log_dir}/{{experiment}}_make_dge_design.log",
    output: f"{rna_dir}/models/{{experiment}}/{{experiment}}_design.rds",
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
        gtf = lambda wildcards: f"{ref_dir}/{rna_map[wildcards.experiment]['build']}_wtrans.gtf.gz",
    log: f"{log_dir}/{{experiment}}_make_salmon_txi.log",
    output: f"{rna_dir}/models/{{experiment}}/{{experiment}}_txi.rds",
    params:
        script = rna_script_dir + "/make_salmon_txi.R",
    shell:
        """
        Rscript {params.script} \
        {input.gtf} \
        "{input.salmon}" \
        {output} > {log} 2>&1
        """

rule norm_txi_edger:
    input:
        design = f"{rna_dir}/models/{{experiment}}/{{experiment}}_design.rds",
        txi = f"{rna_dir}/models/{{experiment}}/{{experiment}}_txi.rds",
    log: f"{log_dir}/{{experiment}}_norm_txi_edger.log",
    output:
        dge = f"{rna_dir}/models/{{experiment}}_edger/{{experiment}}dge.rds",
        glm = f"{rna_dir}/models/{{experiment}}_edger/{{experiment}}_fit.rds",
        cpm = f"{rna_dir}/models/{{experiment}}_edger/{{experiment}}_cpm.tsv",
    params: script = f"{rna_script_dir}/norm_txi_edger.R",
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
        cpm = f"{rna_dir}/models/{{experiment}}_edger/{{experiment}}_cpm.tsv",
        libraries_full = libraries_full_rds,
    log: f"{log_dir}/{{experiment}}_make_cpm_pca.log",
    output:
        f"{rna_dir}/models/{{experiment}}/{{experiment}}_pca.png",
        f"{rna_dir}/models/{{experiment}}/{{experiment}}_pca.svg",
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
