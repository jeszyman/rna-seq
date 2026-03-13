#!/usr/bin/env bash
set -euo pipefail

# tools/get_test_data.sh
# AUTO-GENERATED — DO NOT EDIT DIRECTLY (tangled from rna-seq.org)
#
# Generate or download test data for the RNA-seq biopipe module.
# Creates synthetic FASTQ files and minimal reference FASTA/GTF stubs
# for pipeline DAG validation and CI testing.

parse_args() {
    declare -g script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    declare -g repo_dir="$(cd "${script_dir}/.." && pwd)"
    declare -g out_dir="${repo_dir}/tests/full/inputs"

    declare -g nreads=1000

    # Mouse refs
    declare -g mm_fa_url="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr19.fa.gz"
    declare -g mm_gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
    declare -g out_mm_fa="${out_dir}/mm10_chr19.fa.gz"
    declare -g out_mm_gtf="${out_dir}/mm10_chr19.gtf.gz"

    # Human refs
    declare -g hs_fa_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
    declare -g hs_gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"
    declare -g out_hs_fa="${out_dir}/hg38_chr22.fa.gz"
    declare -g out_hs_gtf="${out_dir}/hg38_chr22.gtf.gz"

    # Library naming
    declare -g ALL_LIBS=( lib001 lib002 lib003 lib004 lib005 lib006 lib007 lib008 )
}

clean_inputs_dir() {
    echo "[get_test_data] Cleaning output directory: ${out_dir}"
    rm -rf "${out_dir}"
    mkdir -p "${out_dir}"
}

ensure_gitignore() {
    local gi="${repo_dir}/.gitignore"
    if ! grep -q "tests/full/\*" "$gi" 2>/dev/null; then
        echo "" >> "$gi"
        echo "# --- test data (auto-managed by get_test_data.sh) ---" >> "$gi"
        echo "tests/full/*" >> "$gi"
        echo "!tests/full/inputs/" >> "$gi"
        echo "!tests/full/inputs/**" >> "$gi"
    fi
}

generate_synthetic_fastqs() {
    echo "[get_test_data] Generating synthetic paired-end FASTQs (${nreads} reads each)"
    for lib in "${ALL_LIBS[@]}"; do
        for read in R1 R2; do
            python3 -c "
import random; random.seed(hash('${lib}_${read}'))
bases='ACGT'
for i in range(${nreads}):
    seq=''.join(random.choice(bases) for _ in range(100))
    qual='I'*100
    print(f'@${lib}_${read}_{i}/1')
    print(seq)
    print('+')
    print(qual)
" | gzip > "${out_dir}/${lib}_${read}.fastq.gz"
        done
        echo "  Created ${lib}_R1.fastq.gz and ${lib}_R2.fastq.gz"
    done
}

generate_synthetic_refs() {
    echo "[get_test_data] Generating minimal synthetic reference FASTA and GTF stubs"

    # Mouse chr19 minimal FASTA
    python3 -c "
import random; random.seed(19)
bases='ACGT'
seq = ''.join(random.choice(bases) for _ in range(10000))
print('>chr19')
for i in range(0, len(seq), 80):
    print(seq[i:i+80])
" | gzip > "${out_mm_fa}"

    # Mouse chr19 minimal GTF (a few synthetic genes)
    cat <<'MGTF' | gzip > "${out_mm_gtf}"
chr19	synthetic	gene	1000	2000	.	+	.	gene_id "ENSMUSG00000001"; gene_name "TestGene1";
chr19	synthetic	transcript	1000	2000	.	+	.	gene_id "ENSMUSG00000001"; transcript_id "ENSMUST00000001"; gene_name "TestGene1";
chr19	synthetic	exon	1000	1500	.	+	.	gene_id "ENSMUSG00000001"; transcript_id "ENSMUST00000001"; gene_name "TestGene1";
chr19	synthetic	exon	1600	2000	.	+	.	gene_id "ENSMUSG00000001"; transcript_id "ENSMUST00000001"; gene_name "TestGene1";
chr19	synthetic	gene	3000	5000	.	-	.	gene_id "ENSMUSG00000002"; gene_name "TestGene2";
chr19	synthetic	transcript	3000	5000	.	-	.	gene_id "ENSMUSG00000002"; transcript_id "ENSMUST00000002"; gene_name "TestGene2";
chr19	synthetic	exon	3000	3800	.	-	.	gene_id "ENSMUSG00000002"; transcript_id "ENSMUST00000002"; gene_name "TestGene2";
chr19	synthetic	exon	4200	5000	.	-	.	gene_id "ENSMUSG00000002"; transcript_id "ENSMUST00000002"; gene_name "TestGene2";
chr19	synthetic	gene	6000	8000	.	+	.	gene_id "ENSMUSG00000003"; gene_name "TestGene3";
chr19	synthetic	transcript	6000	8000	.	+	.	gene_id "ENSMUSG00000003"; transcript_id "ENSMUST00000003"; gene_name "TestGene3";
chr19	synthetic	exon	6000	7000	.	+	.	gene_id "ENSMUSG00000003"; transcript_id "ENSMUST00000003"; gene_name "TestGene3";
chr19	synthetic	exon	7500	8000	.	+	.	gene_id "ENSMUSG00000003"; transcript_id "ENSMUST00000003"; gene_name "TestGene3";
MGTF

    # Human chr22 minimal FASTA
    python3 -c "
import random; random.seed(22)
bases='ACGT'
seq = ''.join(random.choice(bases) for _ in range(10000))
print('>chr22')
for i in range(0, len(seq), 80):
    print(seq[i:i+80])
" | gzip > "${out_hs_fa}"

    # Human chr22 minimal GTF
    cat <<'HGTF' | gzip > "${out_hs_gtf}"
chr22	synthetic	gene	1000	2000	.	+	.	gene_id "ENSG00000001"; gene_name "TestGeneA";
chr22	synthetic	transcript	1000	2000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "TestGeneA";
chr22	synthetic	exon	1000	1500	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "TestGeneA";
chr22	synthetic	exon	1600	2000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "TestGeneA";
chr22	synthetic	gene	3000	5000	.	-	.	gene_id "ENSG00000002"; gene_name "TestGeneB";
chr22	synthetic	transcript	3000	5000	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; gene_name "TestGeneB";
chr22	synthetic	exon	3000	3800	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; gene_name "TestGeneB";
chr22	synthetic	exon	4200	5000	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; gene_name "TestGeneB";
chr22	synthetic	gene	6000	8000	.	+	.	gene_id "ENSG00000003"; gene_name "TestGeneC";
chr22	synthetic	transcript	6000	8000	.	+	.	gene_id "ENSG00000003"; transcript_id "ENST00000003"; gene_name "TestGeneC";
chr22	synthetic	exon	6000	7000	.	+	.	gene_id "ENSG00000003"; transcript_id "ENST00000003"; gene_name "TestGeneC";
chr22	synthetic	exon	7500	8000	.	+	.	gene_id "ENSG00000003"; transcript_id "ENST00000003"; gene_name "TestGeneC";
HGTF

    echo "  Created mm10_chr19.fa.gz, mm10_chr19.gtf.gz"
    echo "  Created hg38_chr22.fa.gz, hg38_chr22.gtf.gz"
}

fetch_refs() {
    # Try to download real refs; fall back to synthetic
    echo "[get_test_data] Attempting to download reference files..."
    if command -v curl &>/dev/null; then
        curl -sLf "${mm_fa_url}" -o "${out_mm_fa}" 2>/dev/null && \
            echo "  Downloaded mm10 chr19 FASTA" || \
            echo "  Download failed, using synthetic refs"
    fi
    # If download failed or produced empty file, use synthetic
    if [[ ! -s "${out_mm_fa}" ]]; then
        generate_synthetic_refs
    fi
}

validate_outputs() {
    echo "[get_test_data] Validating outputs..."
    local ok=true
    for lib in "${ALL_LIBS[@]}"; do
        for read in R1 R2; do
            if [[ ! -s "${out_dir}/${lib}_${read}.fastq.gz" ]]; then
                echo "  MISSING: ${lib}_${read}.fastq.gz"
                ok=false
            fi
        done
    done
    for ref in "${out_mm_fa}" "${out_mm_gtf}" "${out_hs_fa}" "${out_hs_gtf}"; do
        if [[ ! -s "$ref" ]]; then
            echo "  MISSING: $ref"
            ok=false
        fi
    done
    if $ok; then
        echo "[get_test_data] All outputs validated successfully"
    else
        echo "[get_test_data] ERROR: Some outputs are missing"
        exit 1
    fi
}

main() {
    parse_args "$@"
    clean_inputs_dir
    ensure_gitignore
    generate_synthetic_refs
    generate_synthetic_fastqs
    validate_outputs
    echo "[get_test_data] Done"
}

main "$@"
