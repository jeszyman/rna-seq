print("Integration testing snakefile for bulk RNA-seq\n")

# Import common packages
import pandas as pd
import re
import numpy as np





rule all:
    input:
        expand(salmon + "/{library}.quant.sf", library = BULK_RNA_LIBS)

rule symlink_salmon:
    container: rna_container
    input: lambda wildcards: bulk_rna_lib_dict[wildcards.librfary],
    log: logdir + "/{library}_symlink_salmon.log",
    output: salmon = "/{library}.quant.sf",
    params:
        script = rna_scriptdir + "/symlink_salmon.sh"
    shell:
        """
        {params.script} {input} {ouput} &> {log}
        """

include: " <INCLUDE FILE LOCATION (VIA CONFIG PARAM)>"
