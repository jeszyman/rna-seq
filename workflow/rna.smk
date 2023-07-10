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
    input: f"{ref_dir}/{{build}}_protein.gtf",
    log: f"{log_dir}/{{build}}_make_txdb_from_gtf.log",
    output: f"{ref_dir}/{{build}}_protein_txdb",
    params: script = f"{rna_scriptdir}/make_txdb_from_gtf.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        > {log} 2>&1
        """

rule make_salmon_txi:
    input: expand(salmon_dir + "/{library}.sf", library = RNA_LIBS),
    log: logdir + "/{experiment}_make_salmon_txi.log",
    output: rna_dir + "/{experiment}_txi.rdata",
    params:
        script = rna_script_dir + "/make_salmon_txi.R",
        txdb = txdb,
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output} \
        {params.txdb} \
        > {log} 2>&1
        """

rule all_rna_eda:
    container: "/home/jeszyman/sing_containers/atac.1.1.0.sif",
    input: analysis + "/{experiment}_txi.rdata",
    log: logdir + "/{experiment}_rna_eda.log",
    output:
        pca = results + "/figures/{experiment}_pca.pdf",
        rdata = analysis + "/{experiment}_eda.rdata",
    params:
        factor_str = factor_str,
        library_tsv = library_tsv,
        script = rna_scriptdir + "/all_rna_eda.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output.pca} \
        {output.rdata} \
        "{params.factor_str}" \
        {params.library_tsv} \
        > {log} 2>&1
        """
