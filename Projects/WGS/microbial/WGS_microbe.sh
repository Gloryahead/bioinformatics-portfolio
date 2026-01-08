#!/bin/bash
# ================================================================
# WGS Pipeline Script (multi-sample)
# Steps:
#   1. Setup directories
#   2. Download FASTQ files for all accessions
#   3. Run FastQC (pre-trim)
#   4. Trim reads with fastp
#   5. Run FastQC (post-trim)
#   6. Assemble with SPAdes
#   7. Run QUAST
#   8. Screen for AMR genes with abricate
# ================================================================

set -euo pipefail

# ---------------------------
# Define directories
# ---------------------------
base_dir="/home/maa/MAA_science/WGS_microbe_task"
fastq_files="$base_dir/fastq_files"
pre_trim_qc="$base_dir/pre_trim_qc"
trimmed="$base_dir/trimmed"
post_trim_qc="$base_dir/post_trim_qc"
assembly="$base_dir/assembly"
quast_report="$base_dir/quast_report"
AMR="$base_dir/AMR"

# Create directories
mkdir -p "$fastq_files" "$pre_trim_qc" "$trimmed" "$post_trim_qc" "$assembly" "$quast_report" "$AMR"

# ---------------------------
# SRR Accessions
# ---------------------------
accessions=(
SRR27013316 SRR27013315 SRR27013314 SRR27013313 SRR27013312
SRR27013311 SRR27013342 SRR27013340 SRR27013339 SRR27013338
SRR27013343 SRR27013341 SRR27013337 SRR27013335 SRR27013336
SRR27013334 SRR27013333 SRR27013332 SRR27013331 SRR27013330
SRR27013329 SRR27013328 SRR27013327 SRR27013326 SRR27013325
SRR27013323 SRR27013324 SRR27013322 SRR27013321 SRR27013320
SRR27013319 SRR27013318 SRR27013309 SRR27013310 SRR27013317
SRR27013308 SRR27013307 SRR27013305 SRR27013303 SRR27013306
SRR27013304 SRR27013302 SRR27013301 SRR27013300 SRR27013297
)

# ---------------------------
# Loop through each accession
# ---------------------------
for acc in "${accessions[@]}"; do
    echo ">>> Processing $acc ..."

    # ---------------------------
    # Step 1: Download FASTQ files
    # ---------------------------
    prefix=${acc:0:6}   # SRR270
    last2=${acc: -2}    # last 2 digits
    subdir="0${last2}"  # padded
    url_base="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${subdir}/${acc}"

    echo ">>> Downloading $acc ..."
    wget -nc -O "$fastq_files/${acc}_1.fastq.gz" "$url_base/${acc}_1.fastq.gz"
    wget -nc -O "$fastq_files/${acc}_2.fastq.gz" "$url_base/${acc}_2.fastq.gz"

    # ---------------------------
    # Step 2: Pre-trim QC
    # ---------------------------
    fastqc -o "$pre_trim_qc" "$fastq_files/${acc}_1.fastq.gz" "$fastq_files/${acc}_2.fastq.gz"

    # ---------------------------
    # Step 3: Trim reads
    # ---------------------------
    fastp \
      -i "$fastq_files/${acc}_1.fastq.gz" \
      -I "$fastq_files/${acc}_2.fastq.gz" \
      -o "$trimmed/${acc}_1.trim.fastq.gz" \
      -O "$trimmed/${acc}_2.trim.fastq.gz"

    # ---------------------------
    # Step 4: Post-trim QC
    # ---------------------------
    fastqc -o "$post_trim_qc" "$trimmed/${acc}_1.trim.fastq.gz" "$trimmed/${acc}_2.trim.fastq.gz"

    # ---------------------------
    # Step 5: Assembly with SPAdes
    # ---------------------------
    spades_out="$assembly/$acc"
    mkdir -p "$spades_out"

    spades.py --phred-offset 33 \
      -1 "$trimmed/${acc}_1.trim.fastq.gz" \
      -2 "$trimmed/${acc}_2.trim.fastq.gz" \
      -o "$spades_out"

    # ---------------------------
    # Step 6: QUAST Report
    # ---------------------------
    quast.py "$spades_out/contigs.fasta" -o "$quast_report/$acc"

    # ---------------------------
    # Step 7: AMR detection
    # ---------------------------
    abricate "$spades_out/contigs.fasta" > "$AMR/${acc}_amr.tab"

    echo ">>> Finished $acc"
done

echo ">>> Pipeline complete! Results stored in $base_dir"
