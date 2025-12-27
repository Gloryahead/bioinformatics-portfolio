#!/bin/bash

# ==============================================================================
# CUT&RUN Deeptools Track Generation (FINAL VERSION)
# ==============================================================================

# --- Configuration ---
core=8
cores=8 

# Define base directory for the project
base_dir="/your/cluster/path/project_name"
echo "Project Base Directory: $base_dir"

# --- Define Required Paths (CRITICAL: These were missing from your snippet) ---
# Assuming the file structure from previous scripts in /CR12/
bigwig_file="$base_dir/bigwig"
filtered_bam_files="$base_dir/bam/filtered_blacklist"


# BigWig output paths
bigwig_RPKM_10bp="$bigwig_file/bigwig_RPKM_10bp" 
bigwig_RPKM_20bp="$bigwig_file/bigwig_RPKM_20bp"
bigwig_RPKM_25bp="$bigwig_file/bigwig_RPKM_25bp"
bigwig_RPKM_50bp="$bigwig_file/bigwig_RPKM_50bp"

# --- Create Directories ---
echo -e "\n--- Setting up project directories ---"
mkdir -p \
    "$bigwig_RPKM_10bp" "$bigwig_RPKM_20bp" \
    "$bigwig_RPKM_25bp" "$bigwig_RPKM_50bp"
    
echo "Project directories created/verified."
echo "----------------------------------------------------"

# --- Sample List ---
declare -a samples=(
    "CR12CGR" "CR12CI" "CR12CS" "CR12DGR"
    "CR12DI" "CR12DS" "CR13CGR" "CR13CI"
    "CR13CPF" "CR13CS" "CR13DGR" "CR13DI"
    "CR13DPF" "CR13DS" "CR15CI" "CR15CK4me3"
    "CR15CPF" "CR15CS" "CR15DI" "CR15DK4me3"
    "CR15DPF" "CR15DS"
)


# --- Main Processing Loop ---
for sample_id in "${samples[@]}"; do
    echo -e "\n=== Processing sample: $sample_id ==="
    
    # Input BAM file
    sorted_bam_output="$filtered_bam_files/${sample_id}_mm10.filtered.sorted.bam" 

    # Output BigWig files
    bigwig_RPKM_output_10="$bigwig_RPKM_10bp/${sample_id}_RPKM_10bp.bw"
    bigwig_RPKM_output_20="$bigwig_RPKM_20bp/${sample_id}_RPKM_20bp.bw"
    bigwig_RPKM_output_25="$bigwig_RPKM_25bp/${sample_id}_RPKM_25bp.bw"
    bigwig_RPKM_output_50="$bigwig_RPKM_50bp/${sample_id}_RPKM_50bp.bw"

    # --- 1. Generate RPKM BigWig Tracks at Multiple Resolutions ---

    # 1a. Generate RPKM BigWig with bin 10 bp
    echo "Generating 10 bp BigWig..."    

    bamCoverage \
      -b "$sorted_bam_output" \
      -o "$bigwig_RPKM_output_10" \
      --binSize 10 \
      --normalizeUsing RPKM \
      --smoothLength 0 \
      --centerReads \
      --ignoreForNormalization chrM \
      -p "$cores"
    echo "Scaled RPKM (10 bp) complete for $sample_id."

    # 1b. Generate RPKM BigWig with bin 20 bp
    echo "Generating 20 bp BigWig..."

    bamCoverage \
      -b "$sorted_bam_output" \
      -o "$bigwig_RPKM_output_20" \
      --binSize 20 \
      --normalizeUsing RPKM \
      --smoothLength 0 \
      --centerReads \
      --ignoreForNormalization chrM \
      -p "$cores"
    echo "Scaled RPKM (20 bp) complete for $sample_id."

    # 1c. Generate RPKM BigWig with bin 25 bp
    echo "Generating 25 bp BigWig..."
    
    bamCoverage \
      -b "$sorted_bam_output" \
      -o "$bigwig_RPKM_output_25" \
      --binSize 25 \
      --normalizeUsing RPKM \
      --smoothLength 0 \
      --centerReads \
      --ignoreForNormalization chrM \
      -p "$cores"
    echo "Scaled RPKM (25 bp) complete for $sample_id."

    # 1d. Generate RPKM BigWig with bin 50 bp
    echo "Generating 50 bp BigWig..."

    bamCoverage \
      -b "$sorted_bam_output" \
      -o "$bigwig_RPKM_output_50" \
      --binSize 50 \
      --normalizeUsing RPKM \
      --smoothLength 0 \
      --centerReads \
      --ignoreForNormalization chrM \
      -p "$cores"
    echo "Scaled RPKM (50 bp) complete for $sample_id."

done

echo "Pipeline finished."
