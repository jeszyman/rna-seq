input="${1}"
outdir="${2}"
threads="${3}"

fastqc  --outdir $outdir \
        --quiet \
        --threads $threads $input
