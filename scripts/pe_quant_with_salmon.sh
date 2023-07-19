#!/usr/bin/env bash

index="${1}"
read1="${2}"
read2="${3}"
out_dir="${4}"
threads="${5}"

salmon quant \
       --index $index \
       --libType A \
       --mates1 $read1 \
       --mates2 $read2 \
       --output $out_dir \
       --threads $threads \
       --validateMappings
