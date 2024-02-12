######################################
###   RNA-seq Initial Processing   ###
######################################


# - https://www.biostars.org/p/106590/

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
        json = f"{rna_qc_dir}/{{library}}_fastp.json",
        cmd = f"{rna_qc_dir}/{{library}}_fastp.log",
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
    output: f"{rna_qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.zip",
    params:
        out_dir = rna_qc_dir,
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



# Make an experimental design for the list of libraries.


rule make_dge_design:
    input:
        libraries_full = libraries_full_rds,
    log: f"{log_dir}/{{experiment}}_make_dge_design.log",
    output: f"{rna_dir}/models/unadjusted/{{experiment}}/design.rds",
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



# Annotate and summarize counts for salmon files


rule make_salmon_txi:
    input:
        salmon = lambda wildcards: expand(f"{rna_dir}/salmon/{{library}}_{{build}}/quant.sf",
                                          library = rna_map[wildcards.experiment]['libs'],
                                          build = rna_map[wildcards.experiment]['build']),
        gtf = lambda wildcards: f"{ref_dir}/{rna_map[wildcards.experiment]['build']}_wtrans.gtf.gz",
    log: f"{log_dir}/{{experiment}}_make_salmon_txi.log",
    output: f"{rna_dir}/models/unadjusted/{{experiment}}/txi.rds",
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
        design = f"{rna_dir}/models/unadjusted/{{experiment}}/design.rds",
        txi = f"{rna_dir}/models/unadjusted/{{experiment}}/txi.rds",
    log: f"{log_dir}/{{experiment}}_unadjusted_norm_txi_edger.log",
    output:
        dge = f"{rna_dir}/models/unadjusted/{{experiment}}/edger_dge.rds",
        glm = f"{rna_dir}/models/unadjusted/{{experiment}}/edger_fit.rds",
        cpm = f"{rna_dir}/models/unadjusted/{{experiment}}/edger_cpm.tsv",
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



# Makes a logCPM-based PCA plot.


rule make_cpm_pca:
    input:
        cpm = f"{rna_dir}/models/unadjusted/{{experiment}}/edger_cpm.tsv",
        libraries_full = libraries_full_rds,
    log: f"{log_dir}/{{experiment}}_make_cpm_pca.log",
    output:
        f"{rna_dir}/models/unadjusted/{{experiment}}/pca.png",
        f"{rna_dir}/models/unadjusted/{{experiment}}/pca.pdf",
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

rule make_rna_batch_corrections:
    input:
        design = f"{rna_dir}/models/unadjusted/{{experiment}}/design.rds",
        libs = libraries_full_rds,
        txi = f"{rna_dir}/models/unadjusted/{{experiment}}/txi.rds",
    log: f"{log_dir}/{{experiment}}_make_rna_batch_corrections.log",
    output:
        design = f"{rna_dir}/models/combat/{{experiment}}/design.rds",
        dge = f"{rna_dir}/models/combat/{{experiment}}/edger_dge.rds",
        glm = f"{rna_dir}/models/combat/{{experiment}}/edger_fit.rds",
        pca = f"{rna_dir}/models/combat/{{experiment}}/pca.pdf",
    params:
        batch_var = lambda wildcards: rna_map[wildcards.experiment]['batch_var'],
        covars = lambda wildcards: rna_map[wildcards.experiment]['covars'],
        script = f"{rna_script_dir}/make_rna_batch_corrections.R",
    shell:
        """
	cp {input.design} {output.design}
	Rscript {params.script} \
	--batch_var {params.batch_var} \
	--covars {params.covars} \
	--design_rds {output.design} \
	--dge_rds {output.dge} \
	--glm_rds {output.glm} \
	--libraries_full_rds {input.libs} \
	--pdf {output.pca} \
	--txi_rds {input.txi} > {log} 2>&1
	"""

rule make_rlog_matrix:
    input:
        dge = f"{rna_dir}/models/{{adjustment}}/{{experiment}}/edger_dge.rds",
        libs = libraries_full_rds,
    log: f"{log_dir}/{{experiment}}_{{adjustment}}_make_rlog_matrix.log",
    output:
        f"{rna_dir}/models/{{adjustment}}/{{experiment}}/rlog_mat.rds"
    params:
        formula = lambda wildcards: rna_map[wildcards.experiment]['formula'],
        script = f"{rna_script_dir}/make_rlog_matrix.R",
    shell:
        """
        Rscript {params.script} \
        --dge_rds {input.dge} \
        --formula_string "{params.formula}" \
        --libraries_full_rds {input.libs} \
        --out_mat_file {output} \
        > {log} 2>&1
        """

rule model_heatmap_unadjusted:
    input:
        dge = f"{rna_dir}/models/unadjusted/{{experiment}}/edger_dge.rds",
        libs = libraries_full_rds,
    log: f"{log_dir}/{{experiment}}_model_heatmap_unadjusted.log",
    output: f"{rna_dir}/models/unadjusted/{{experiment}}/heat.pdf",
    params:
        fct_str= lambda wildcards: rna_map[wildcards.experiment]['fct_str'],
        nrow = 100,
        script = f"{rna_script_dir}/model_heatmap.R",
    shell:
        """
        Rscript {params.script} \
        --dge_rds {input.dge} \
        --fct_str "{params.fct_str}" \
        --libs_rds {input.libs} \
        --nrow {params.nrow} \
        --pdf_out {output} > {log} 2>&1
        """

rule model_heatmap_combat:
    input:
        dge = f"{rna_dir}/models/combat/{{experiment}}/edger_dge.rds",
        libs = libraries_full_rds,
    log: f"{log_dir}/{{experiment}}_model_heatmap_combat.log",
    output: f"{rna_dir}/models/combat/{{experiment}}/heat.pdf",
    params:
        fct_str= lambda wildcards: rna_map[wildcards.experiment]['fct_str'],
        nrow = 100,
        script = f"{rna_script_dir}/model_heatmap.R",
    shell:
        """
        Rscript {params.script} \
        --dge_rds {input.dge} \
        --fct_str "{params.fct_str}" \
        --libs_rds {input.libs} \
        --nrow {params.nrow} \
        --pdf_out {output} > {log} 2>&1
        """

rule make_edger_contrast_de:
    input:
        design = lambda wildcards: f"{rna_dir}/models/{rna_dge_map[wildcards.contrast]['correction']}/{rna_dge_map[wildcards.contrast]['model']}/design.rds",
        fit = lambda wildcards: f"{rna_dir}/models/{rna_dge_map[wildcards.contrast]['correction']}/{rna_dge_map[wildcards.contrast]['model']}/edger_fit.rds",
        annotation_tsv = lambda wildcards: f"{ref_dir}/{rna_dge_map[wildcards.contrast]['build']}_wtrans_annotation.tsv",
    log: f"{log_dir}/{{contrast}}_make_edger_contrast_de.log",
    output: f"{rna_dir}/contrasts/{{contrast}}/edger_dge.tsv",
    params:
        cohorts_str = lambda wildcards: rna_dge_map[wildcards.contrast]['cohorts_str'],
        script = f"{rna_script_dir}/make_edger_contrast_de.R",
    shell:
        """
        Rscript {params.script} \
        --annotation_tsv {input.annotation_tsv} \
        --cohorts_str "{params.cohorts_str}" \
        --design_rds {input.design} \
        --fit_rds {input.fit} \
        --res_tsv {output} > {log} 2>&1
        """

rule rna_volcano:
    input: f"{rna_dir}/contrasts/{{contrast}}/edger_dge.tsv",
    log: f"{log_dir}/{{contrast}}_rna_volcano.log",
    output: f"{rna_dir}/contrasts/{{contrast}}/volcano.pdf",
    params: script = f"{rna_script_dir}/rna_volcano.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        > {log} 2>&1
        """


# solve s2n

# https://www.google.com/search?q=fgsea+signal+to+noise&oq=fgsea+signal+to+noise&gs_lcrp=EgZjaHJvbWUyBggAEEUYOTIGCAEQLhhA0gEIMzU5NGowajGoAgCwAgA&sourceid=chrome&ie=UTF-8
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1674-0#Sec8

rule gsea_from_edger:
    input: f"{rna_dir}/contrasts/{{contrast}}/edger_dge.tsv",
    log: f"{log_dir}/{{contrast}}_{{pathset}}_gsea_from_edger.log",
    output:
       f"{rna_dir}/contrasts/{{contrast}}/gsea_{{pathset}}.tsv",
       f"{rna_dir}/contrasts/{{contrast}}/gsea_{{pathset}}.xlsx",
    params: script = f"{rna_script_dir}/gsea_from_edger.R",
    shell:
        """
        Rscript {params.script} {input} {wildcards.pathset} {output} > {log} 2>&1
        """
