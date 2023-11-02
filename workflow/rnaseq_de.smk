# Preamble

print("Integration testing snakefile for Post-QC RNA-seq Differential Expression\n")

# Import common packages
import pandas as pd
import re
import numpy as np

# Variable naming



# Functions, miscellaneous



# All rule

rule all:
    input:
        design
        tmm
        ebayes
        dds

# Symlink inputs

rule symlink_rnaseq_de_inputs:
    input:

# Include statements

#include: " <INCLUDE FILE LOCATION (VIA CONFIG PARAM)>"
