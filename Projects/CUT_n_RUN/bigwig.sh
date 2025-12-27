#!/bin/bash

# ==============================================================================
# CUT&RUN Deeptools Track Generation (FINAL VERSION)
# Calculation: SF = 1 / BW% (as per EpiCypher user request)
# ==============================================================================

# --- Configuration ---
cores=8 

# Define base directory for the project
base_dir="/your/cluster/path/project_name"
echo "Project Base Directory: $base_dir"

# --- Define Required Paths (CRITICAL: These were missing from your snippet) ---
# Assuming the file structure from previous scripts in /CR12/
bigwig_file="$base_dir/bigwig"
sam_ecoli_sum="$base_dir/sam_ecoli/sam_ecoli_summary" 
sam_mm10_sum="$base_dir/sam_mm10/sam_mm10_summary"
filtered_bam_files="$base_dir/bam/filtered_blacklist"
genome_size="/path/to/your/genome_size/mm10.chrom.sizes" # Assuming static path

# BigWig output paths
bigwig_bin_10="$bigwig_file/bigwig_bin_10" 
bigwig_bin_25="$bigwig_file/bigwig_bin_25"
bigwig_bin_50="$bigwig_file/bigwig_bin_50"

# --- Create Directories ---
echo -e "\n--- Setting up project directories ---"
mkdir -p \
    "$bigwig_bin_10" "$bigwig_bin_25" \
    "$bigwig_bin_50" 
    
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

    # Define sample-specific files
    sam_ecoli_sum_output="$sam_ecoli_sum/${sample_id}_ecoli.txt"
    sam_mm10_sum_output="$sam_mm10_sum/${sample_id}_mm10.txt"
    
    # Input BAM file
    sorted_bam_output="$filtered_bam_files/${sample_id}_mm10.filtered.sorted.bam" 

    # Output BigWig files
    bigwig_scaled_output_10="$bigwig_bin_10/${sample_id}_scaled_10.bw"
    bigwig_scaled_output_25="$bigwig_bin_25/${sample_id}_scaled_25.bw"
    bigwig_scaled_output_50="$bigwig_bin_50/${sample_id}_scaled_50.bw"

    # --- 1. Spike-in Normalization Calculation ---
    echo "--- Calculating Spike-in Bandwidth and Scaling Factor... ---"

    if [ ! -f "$sam_ecoli_sum_output" ] || [ ! -f "$sam_mm10_sum_output" ]; then
        echo "ERROR: Bowtie2 summary files missing. Cannot calculate scaling factor. Using SF=1.0."
        scale_factor=1.0
    else
        # 1. ROBUSTLY extract Unique Paired Reads
        ecoli_unique_pairs=$(grep 'aligned concordantly exactly 1 time' "$sam_ecoli_sum_output" | sed 's/^[ \t]*//' | awk '{print $1}')
        mm10_unique_pairs=$(grep 'aligned concordantly exactly 1 time' "$sam_mm10_sum_output" | sed 's/^[ \t]*//' | awk '{print $1}')
        
        ecoli_unique_pairs=${ecoli_unique_pairs//[^0-9]/} 
        mm10_unique_pairs=${mm10_unique_pairs//[^0-9]/}
        ecoli_unique_pairs=${ecoli_unique_pairs:-0}
        mm10_unique_pairs=${mm10_unique_pairs:-0}
        
        echo "N_ecoli (Unique Pairs): $ecoli_unique_pairs"
        echo "N_mm10 (Unique Pairs): $mm10_unique_pairs"

        total_unique_pairs=$((ecoli_unique_pairs + mm10_unique_pairs))
        
        # 2. Calculate Spike-in Bandwidth (BW) AS PERCENTAGE
        if [ "$total_unique_pairs" -eq 0 ]; then
            echo "ERROR: Total unique pairs is zero. Using SF=1.0."
            spike_bandwidth_percent=0.0
            scale_factor=1.0
        else
            # Calculate BW as a percentage: (N_ecoli / N_total) * 100
            spike_bandwidth_percent=$(awk "BEGIN { printf \"%.6f\", ($ecoli_unique_pairs / $total_unique_pairs) * 100 }")
            
            # 3. Calculate Scaling Factor (SF = 1 / BW%) - AS REQUESTED
            if (( $(echo "$spike_bandwidth_percent == 0" | bc -l) )); then
                scale_factor=1.0
            else
                # SF = 1 / BW(%) - Using the exact, requested formula.
                scale_factor=$(awk "BEGIN { printf \"%.6f\", 1.0 / $spike_bandwidth_percent }")
            fi
        fi
    fi 

    echo "Spike-in Bandwidth (BW%): $spike_bandwidth_percent%"
    echo "Calculated Scaling Factor (SF): $scale_factor"
    
    # Check if BAM file exists before proceeding with bamCoverage
    if [ ! -f "$sorted_bam_output" ]; then
        echo "FATAL ERROR: Sorted BAM file not found at $sorted_bam_output. Skipping track generation."
        continue
    fi

    # --- 2. Generate scaled BigWig Tracks at Multiple Resolutions ---

    # 2a. Generate scaled BigWig with bin 10 bp
    echo "Generating 10 bp BigWig..."
    # The crucial command is run here, setting $TMPDIR before the command if needed for I/O safety
    bamCoverage -b "$sorted_bam_output" --scaleFactor "$scale_factor" --binSize 10 -o "$bigwig_scaled_output_10" -p "$cores"
    echo "Scaled BigWig (10 bp) complete for $sample_id."
    
    # 2b. Generate scaled BigWig with bin 25 bp
    echo "Generating 25 bp BigWig..."
    bamCoverage -b "$sorted_bam_output" --scaleFactor "$scale_factor" --binSize 25 -o "$bigwig_scaled_output_25" -p "$cores"
    echo "Scaled BigWig (25 bp) complete for $sample_id."
    
    # 2c. Generate scaled BigWig with bin 50 bp
    echo "Generating 50 bp BigWig..."
    bamCoverage -b "$sorted_bam_output" --scaleFactor "$scale_factor" --binSize 50 -o "$bigwig_scaled_output_50" -p "$cores"
    echo "Scaled BigWig (50 bp) complete for $sample_id."
    
done
