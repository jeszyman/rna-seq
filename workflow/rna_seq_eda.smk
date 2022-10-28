

rule make_salmon_txi:
    container: rna_container,
    input: expand(salmon + "/{library}.quant.sf", library = BULK_RNA_LIBS),
    log: logdir + "/{experiment}_make_salmon_txi.log",
    output: analysis + "/{experiment}_txi.rdata",
    params:
        script = rna_scriptdir + "/make_salmon_txi.R",
        txdb = "TxDb.Mmusculus.UCSC.mm10.ensGene",
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
        library_tsv = inputs + "/libraries.tsv",
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
