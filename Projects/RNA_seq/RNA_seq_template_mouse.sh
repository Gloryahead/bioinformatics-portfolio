#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Directory Structure Visual Representation ---
# The project will be organized under the 'base_dir' as follows:
#
# MAA_science/
# └── RNA_seq_human/                  # base_dir
#     ├── fastq_files/                # Raw FASTQ input files
#     │   └── SRRxxxxxxx_1.fastq.gz
#     │   └── SRRxxxxxxx_2.fastq.gz
#     │   └── ...
#     ├── hg38_dir/                   # Genome reference and annotation
#     │   ├── Homo_sapiens_assembly38.fasta
#     │   ├── h_sapiens.gff3
#     │   └── genomeIndex/            # STAR genome index files
#     │       └── SA.idx
#     │       └── ...
#     ├── qc_raw_reads/               # FastQC reports for raw reads
#     │   └── SRRxxxxxxx_1_fastqc.html
#     │   └── SRRxxxxxxx_2_fastqc.html
#     │   └── ...
#     ├── multiqc_raw/                # MultiQC report aggregating raw QC and fastp stats
#     │   ├── multiqc_report_raw_reads.html
#     │   └── multiqc_report_fastp.html
#     ├── trimmed_fastq/              # Trimmed FASTQ files and fastp reports
#     │   ├── SRRxxxxxxx_trimmed_1.fastq.gz
#     │   ├── SRRxxxxxxx_trimmed_2.fastq.gz
#     │   ├── SRRxxxxxxx_fastp.html
#     │   └── SRRxxxxxxx_fastp.json
#     │   └── ...
#     ├── qc_trimmed_reads/           # FastQC reports for trimmed reads
#     │   └── SRRxxxxxxx_trimmed_1_fastqc.html
#     │   └── SRRxxxxxxx_trimmed_2_fastqc.html
#     │   └── ...
#     ├── multiqc_trimmed/            # MultiQC report aggregating trimmed QC
#     │   └── multiqc_report_trimmed_reads.html
#     ├── star_mapped/                # STAR alignment output (BAMs, logs, etc.)
#     │   ├── SRRxxxxxxx_Aligned.sortedByCoord.out.bam
#     │   └── SRRxxxxxxx_Log.final.out
#     │   └── ...
#     ├── IGV/                        # Indexed BAMs ready for IGV visualization
#     │   ├── SRRxxxxxxx_Aligned.sortedByCoord.out.bam
#     │   ├── SRRxxxxxxx_Aligned.sortedByCoord.out.bam.bai
#     │   ├── SRRxxxxxxx_head.txt
#     │   └── SRRxxxxxxx_flagstat.txt
#     │   └── ...
#     └── counts/                     # FeatureCounts gene quantification output
#         └── SRRxxxxxxx_counts.txt
#         └── ...
#
# -----------------------------------------------


# --- Configuration Variables ---
# Define base directory for the project
base_dir="/xdisk/clsmith1/maarowosegbe/RNA-seq/Test_120525"

# Genome and Annotation references
mm10_genome_dir="/xdisk/clsmith1/maarowosegbe/RNA-seq/STAR_mm10_dir"
STAR_index="$mm10_genome_dir/genomeIndex"
mm10_fasta="$mm10_genome_dir/GRCm38.primary_assembly.genome.fa"
gene_annotation="$mm10_genome_dir/GRCm38_annot.gff3" # Renamed for clarity and consistency

# Input/Output directories
fastq_files="$base_dir/fastq_files"
trimmed_reads="$base_dir/trimmed_fastq"
mapped_dir="$base_dir/star_mapped"
counts_dir="$base_dir/counts"
IGV_dir="$base_dir/IGV" # For BAMs ready for IGV visualization

# Quality Control directories
pre_trim_qc="$base_dir/qc_raw_reads"       # Renamed for clarity
post_trim_qc="$base_dir/qc_trimmed_reads"  # Renamed for clarity
multiqc_pre_trim="$base_dir/multiqc_raw"   # Renamed for clarity
multiqc_post_trim="$base_dir/multiqc_trimmed" # Renamed for clarity

# Setting threads and RAM for STAR
THREADS=12   # Leave a few cores for the system
RAM_LIMIT=200000000000 # 200GB is usually enough for Human/Mouse

# --- Create Directories ---
echo "--- Setting up project directories ---"
mkdir -p \
    "$pre_trim_qc" "$multiqc_pre_trim" \
    "$trimmed_reads" "$post_trim_qc" "$multiqc_post_trim" \
    "$STAR_index" "$mapped_dir" "$IGV_dir" "$counts_dir" \
    "$hg38_genome_dir" # Ensure genome directory exists
echo "Project directories created/verified."

# --- Sample List ---
# Declare an array of sample IDs
declare -a samples=(
    "SRR28119110" "SRR28119111"
    "SRR28119112" "SRR28119113"
)


# --- Check for Reference Files ---
echo "--- Checking for genome reference and annotation files ---"
if [ ! -f "$hg38_fasta" ]; then
    echo "ERROR: Genome FASTA file not found at $hg38_fasta. Please place it there."
    exit 1
fi
echo "Genome FASTA file found: $hg38_fasta."


# --- Step 0: Download Genome Annotation (if not already present) ---
if [ ! -f "$gene_annotation" ]; then
    echo "--- Downloading genome annotation (Homo_sapiens.GRCh38.114.gff3.gz)... ---"
    # Using curl for potentially better progress feedback and error handling
    curl -o "$mm10_genome_dir/gencode.vM25.annotation.gff3.gz" \
         "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gff3.gz" \
         --silent --show-error --fail
    echo "Unzipping annotation file..."
    gunzip "$mm10_genome_dir/gencode.vM25.annotation.gff3.gz"
    mv "$mm10_genome_dir/gencode.vM25.annotation.gff3" "$gene_annotation"
    echo "Genome annotation prepared: $gene_annotation"
else
    echo "Genome annotation file already exists: $gene_annotation. Skipping download."
fi


# --- Step 1: Build STAR Genome Index (if not already built) ---
# This step only needs to be run once per genome.
# It's better to put it outside the loop to avoid redundant operations.
if [ ! -d "$STAR_index" ] || [ ! -f "$STAR_index/SA.idx" ]; then
    echo "--- Building STAR genome index (this may take a while)... ---"
    echo "Deleting existing (potentially incomplete) index directory to ensure clean build..."
    rm -rf "$STAR_index" # Ensure a clean slate for index generation
    mkdir -p "$STAR_index" # Recreate the directory

    STAR --runMode genomeGenerate \
         --runThreadN $THREADS \
         --genomeDir "$STAR_index" \
         --genomeFastaFiles "$mm10_fasta" \
         --sjdbGTFfile "$gene_annotation" \
         --sjdbOverhang 100 \
         --limitGenomeGenerateRAM $RAM_LIMIT
    echo "STAR genome index built successfully in $STAR_index."
else
    echo "STAR genome index already exists and appears complete in $STAR_index. Skipping generation."
fi


# --- Main Processing Loop ---
for sample_id in "${samples[@]}"; do
    echo -e "\n=== Processing sample: $sample_id ==="

    # Define sample-specific input/output files
    fastq_r1="$fastq_files/${sample_id}_1.fastq.gz"
    fastq_r2="$fastq_files/${sample_id}_2.fastq.gz"
    trimmed_r1="$trimmed_reads/${sample_id}_trimmed_1.fastq.gz"
    trimmed_r2="$trimmed_reads/${sample_id}_trimmed_2.fastq.gz"

    # New variables for the uncompressed files
    #trimmed_r1_unzipped="$trimmed_reads/${sample_id}_trimmed_1.fastq"
    #trimmed_r2_unzipped="$trimmed_reads/${sample_id}_trimmed_2.fastq"

    star_prefix="$mapped_dir/${sample_id}_"
    mapped_bam="${star_prefix}Aligned.sortedByCoord.out.bam"
    indexed_bam_path="$IGV_dir/${sample_id}_Aligned.sortedByCoord.out.bam" # Final BAM path after indexing
    counts_file="$counts_dir/${sample_id}_counts.txt"

    # --- Pre-check: Ensure raw FASTQ files exist ---
    if [ ! -f "$fastq_r1" ] || [ ! -f "$fastq_r2" ]; then
        echo "WARNING: Raw FASTQ files for $sample_id not found ($fastq_r1, $fastq_r2). Skipping sample."
        continue # Skip to the next sample
    fi
    echo "Raw FASTQ files found for $sample_id."

    # --- Step 2: Quality Control before Trimming (FastQC) ---
    echo "--- Running FastQC on raw reads for $sample_id... ---"
    fastqc "$fastq_r1" "$fastq_r2" -o "$pre_trim_qc"
    echo "FastQC for raw reads completed."

    # --- Step 3: Adapter Trimming (fastp) ---
    echo "--- Running adapter trimming with fastp for $sample_id... ---"
    fastp -i "$fastq_r1" -o "$trimmed_r1" \
          -I "$fastq_r2" -O "$trimmed_r2" \
          --detect_adapter_for_pe \
          --json "$trimmed_reads/${sample_id}_fastp.json" \
          --html "$trimmed_reads/${sample_id}_fastp.html" \
          --thread 8 # Adjust threads as appropriate for your system
    echo "Trimming for $sample_id completed."

    # --- NEW STEP: Decompress trimmed FASTQ files for STAR ---
    echo "--- Decompressing trimmed FASTQ files for $sample_id... ---"
    #gunzip "$trimmed_r1"
    #gunzip "$trimmed_r2"
    echo "Decompression for $sample_id completed."

    # --- Step 4: Quality Control after Trimming (FastQC) ---
    echo "--- Running FastQC on trimmed reads for $sample_id... ---"
    fastqc "$trimmed_r1_unzipped" "$trimmed_r2_unzipped" -o "$post_trim_qc"
    echo "FastQC for trimmed reads completed."

    # --- Step 5: Align trimmed FASTQ files (STAR) ---
    echo "--- Mapping trimmed FASTQ files with STAR for $sample_id... ---"
    STAR --genomeDir "$STAR_index" \
         --runThreadN $THREADS \
         --readFilesIn "$trimmed_r1" "$trimmed_r2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$star_prefix" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --quantMode GeneCounts
    echo "Mapping for $sample_id completed. BAM saved to $mapped_bam."

    # --- Step 6: Prepare BAM for IGV (Copy and Index) ---
    echo "--- Preparing BAM for IGV for $sample_id... ---"
    # Copy the sorted BAM to the IGV directory
    cp "$mapped_bam" "$indexed_bam_path"
    # Index the BAM file for IGV visualization
    samtools index "$indexed_bam_path"
    echo "BAM for IGV prepared and indexed: $indexed_bam_path and ${indexed_bam_path}.bai"

    # --- Step 7: Quick BAM check (optional, but good for verification) ---
    echo "--- Performing quick BAM checks for $sample_id... ---"
    echo "Viewing first 10 lines of BAM (saved to $IGV_dir/${sample_id}_head.txt):"
    samtools view "$indexed_bam_path" | head > "$IGV_dir/${sample_id}_head.txt"

    echo "Counting reads in BAM (saved to $IGV_dir/${sample_id}_count.txt):"
    samtools view -c "$indexed_bam_path" > "$IGV_dir/${sample_id}_count.txt"

    echo "Flagstat summary (saved to $IGV_dir/${sample_id}_flagstat.txt):"
    samtools flagstat "$indexed_bam_path" > "$IGV_dir/${sample_id}_flagstat.txt"
    echo "BAM checks completed."

done

# --- FINAL STEP: Combined Gene Quantification (for all samples) ---

# Define the final output file and capture all indexed BAMs
final_counts_file="$base_dir/final_gene_counts_matrix.txt"
all_bams="$IGV_dir/*_Aligned.sortedByCoord.out.bam"

# Use a checkpoint to ensure this computationally intensive step is only run once
if [ ! -f "$final_counts_file" ]; then
    echo -e "\n--- Running final featureCounts for ALL SAMPLES to create matrix ---"

    # NOTE: The -T 8 flag is added for speed (using 8 threads)
    featureCounts \
        -T 8 \
        -p \
        -t exon \
        -g gene_id \
        -a "$gene_annotation" \
        -o "$final_counts_file" \
        $all_bams # Wildcard expands to include all BAM files

    echo "Combined gene count matrix generated: $final_counts_file"
else
    echo -e "\nCombined gene count matrix already exists. Skipping final quantification."
fi


# --- Final MultiQC Reports ---
echo -e "\n--- Generating MultiQC reports ---"

echo "Running MultiQC on pre-trim FastQC results and fastp reports..."
multiqc "$pre_trim_qc" "$trimmed_reads" -o "$multiqc_pre_trim" -n "multiqc_report_raw_fastp.html"
# The above command now combines pre-trim FastQC and fastp reports into one MultiQC report.

echo "Running MultiQC on post-trim FastQC results..."
multiqc "$post_trim_qc" -o "$multiqc_post_trim" -n "multiqc_report_trimmed_reads.html"

echo "--- All samples processed. Transfer counts files for DE analysis. ---"
echo "--- Don't forget to check the MultiQC reports for overall quality assessment! ---"
echo "Reports are located in: $multiqc_pre_trim/multiqc_report_raw_fastp.html and $multiqc_post_trim/multiqc_report_trimmed_reads.html"
