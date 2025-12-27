#!/bin/bash

# ==============================================================================
# MACS2 Peak Calling Script (Sharp and Broad for applicable factors)
# Supports: CR12, CR13, CR15 samples (TFs/Histones/Scaffolds)
# Method: Sharp call for ALL samples; Broad call for S and PF factors.
# ==============================================================================

# --- Configuration & Paths ---
cores=8

# --- Project Paths (Based on user input) ---
base_dir="/path/to/project" 
# NOTE: Ensure all BAM files for all samples are located here:
filtered_bam_files="$base_dir/bam/filtered_blacklist"

# Output Directory
peak_base="$base_dir/peak_calling_combined"
macs2_out="$peak_base/MACS2"

# Global Configs
MACS2_GSIZE="2.65e9" # Mus musculus effective genome size (mm10)

mkdir -p "$macs2_out"
echo "Output directory created: $macs2_out"

# --- Sample Definitions ---
# All unique factor/histone samples from your list
declare -a TARGET_SAMPLES=(
    "CR12CGR" "CR12CS" "CR12DGR" "CR12DS" 
    "CR13CGR" "CR13CS" "CR13CPF" "CR13DGR" "CR13DS" "CR13DPF" 
    "CR15CK4me3" "CR15CS" "CR15CPF" "CR15DK4me3" "CR15DS" "CR15DPF"
)

# --- Helper Function to get Control ID and Factor Type ---
get_control_info() {
    local sample_id=$1
    local control_id=""
    local is_broad="NO"

    # 1. Determine Control ID based on prefix (CRXX + C/D)
    if [[ "$sample_id" =~ ^(CR12|CR13)C ]]; then
        control_id="${BASH_REMATCH[1]}CI" # e.g., CR12CI or CR13CI
    elif [[ "$sample_id" =~ ^(CR12|CR13)D ]]; then
        control_id="${BASH_REMATCH[1]}DI" # e.g., CR12DI or CR13DI
    elif [[ "$sample_id" =~ ^CR15C ]]; then
        control_id="CR15CI"
    elif [[ "$sample_id" =~ ^CR15D ]]; then
        control_id="CR15DI"
    else
        echo "UNKNOWN_CONTROL NO"
        return
    fi

    # 2. Determine if the factor is a BROAD type (S=Sin3A, PF=PF1)
    if [[ "$sample_id" =~ (S|PF)$ ]]; then
        is_broad="YES"
    fi

    echo "$control_id $is_broad"
}

# --- FUNCTION: MACS2 Peak Calling ---
run_macs2() {
    local sample_id=$1
    local control_id=$2
    local is_broad=$3 # "YES" or "NO"
    
    local treatment_bam="$filtered_bam_files/${sample_id}_mm10.filtered.sorted.bam"
    local control_bam="$filtered_bam_files/${control_id}_mm10.filtered.sorted.bam"
    
    echo -e "\n--- Running MACS2 for $sample_id (Control: $control_id) ---"

    # Expert Check: Ensure input files exist before starting MACS2
    if [ ! -f "$treatment_bam" ] || [ ! -f "$control_bam" ]; then
        echo "ERROR: Missing BAM file(s). Treatment: $treatment_bam or Control: $control_bam. Skipping MACS2."
        return 1
    fi
    
    # 1. SHARP Peak Call (Mandatory for ALL factors)
    echo "  -> 1. MACS2 SHARP Call (Standard Peak Detection)"
    # Default behavior is sharp, no special broad flags needed.
    macs2 callpeak -t "$treatment_bam" -c "$control_bam" \
        -f BAMPE -g "$MACS2_GSIZE" -n "${sample_id}_MACS2_SHARP" \
        --outdir "$macs2_out" -q 0.1 --keep-dup all 

    # 2. BROAD Peak Call (Only if factor type is determined to be BROAD)
    if [ "$is_broad" == "YES" ]; then
        echo "  -> 2. MACS2 BROAD Call (For Scaffold/Corepressor)"
        # Use --broad and --broad-cutoff flags
        macs2 callpeak -t "$treatment_bam" -c "$control_bam" \
            -f BAMPE -g "$MACS2_GSIZE" -n "${sample_id}_MACS2_BROAD" \
            --outdir "$macs2_out" -q 0.1 --broad --broad-cutoff 0.1 --keep-dup all
    fi
}


# ==============================================================================
# EXECUTION
# ==============================================================================

for sample in "${TARGET_SAMPLES[@]}"; do
    # Get the correct control ID and whether it's a broad factor
    INFO=($(get_control_info "$sample"))
    CONTROL="${INFO[0]}"
    IS_BROAD="${INFO[1]}"

    if [ "$CONTROL" == "UNKNOWN_CONTROL" ]; then
        echo "Skipping $sample: Could not determine control."
        continue
    fi
    
    echo -e "\n======================================================="
    echo "Processing Sample: **$sample** (Control: **$CONTROL**)"
    echo "Requires Broad Call? **$IS_BROAD**"
    echo "======================================================="
    
    # Run the MACS2 function
    run_macs2 "$sample" "$CONTROL" "$IS_BROAD"
done

echo -e "\n\n--- MACS2 Peak Calling Complete for All Samples (Sharp & Conditional Broad) ---"
echo "Results saved in: $macs2_out"
