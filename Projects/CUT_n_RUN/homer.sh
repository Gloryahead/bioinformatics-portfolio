#!/bin/bash

# ==============================================================================
# HOMER CUT&RUN/ChIP-seq Analysis Master Script (FINAL VERSION)
# Comprehensive analysis: Tag Dir, Peaks (Sharp/Broad), Motifs, SE, Annotation, 
# Overlap, Quantification, BigWig, and MultiWig Hub Generation (Phases 1-9).
# ==============================================================================

# --- Configuration & Paths ---
cores=8
GENOME_VERSION="/groups/clsmith1/maarowosegbe/micromamba/envs/CUT_n_RUN_env/bin/data/genomes/mm10" # Mus musculus genome
NORMALIZATION_FACTOR="1e7" # 10 million tags for normalization

# --- Project Paths ---
base_dir="/xdisk/clsmith1/maarowosegbe/CR12_13_15"
filtered_bam_files="$base_dir/bam/filtered_blacklist"

# Output Directories
homer_base="$base_dir/peak_calling_combined/HOMER"
tag_dirs="$homer_base/tag_directories"
peaks_out="$homer_base/peaks"
motif_out="$homer_base/motifs"
se_out="$homer_base/super_enhancers"
annotation_out="$homer_base/annotation" 
overlap_out="$homer_base/overlap_quantification" 
multi_overlap_out="$homer_base/multi_factor_overlap" 

# PHASE 9: UCSC Hub Configuration
# SETTING USER-PROVIDED S3 ENDPOINT URL:
HUB_WEB_DIR="$homer_base/UCSC_Hub_Files"
HUB_URL="http://maa-genomics-hub-2025.s3-website.us-east-2.amazonaws.com" 

mkdir -p "$tag_dirs" "$peaks_out" "$motif_out" "$se_out" "$annotation_out" "$overlap_out" "$multi_overlap_out" "$HUB_WEB_DIR"
echo "Output directories created. Hub web directory placeholder: $HUB_WEB_DIR"

# --- Sample Definitions ---
declare -a TARGET_SAMPLES=(
    "CR12CGR" "CR12CS" "CR12DGR" "CR12DS"
    "CR13CGR" "CR13CS" "CR13CPF" "CR13DGR" "CR13DS" "CR13DPF"
    "CR15CK4me3" "CR15CS" "CR15CPF" "CR15DK4me3" "CR15DS" "CR15DPF"
)
declare -a CONTROL_SAMPLES=("CR12CI" "CR12DI" "CR13CI" "CR13DI" "CR15CI" "CR15DI")

# --- Helper Function: Get Control ID and Broad Factor Status ---
get_control_info() {
    local sample_id=$1
    local control_id=""
    local is_broad="NO"

    if [[ "$sample_id" =~ ^(CR12|CR13)C ]]; then
        control_id="${BASH_REMATCH[1]}CI"
    elif [[ "$sample_id" =~ ^(CR12|CR13)D ]]; then
        control_id="${BASH_REMATCH[1]}DI"
    elif [[ "$sample_id" =~ ^CR15C ]]; then
        control_id="CR15CI"
    elif [[ "$sample_id" =~ ^CR15D ]]; then
        control_id="CR15DI"
    else
        echo "UNKNOWN_CONTROL NO"
        return
    fi

    if [[ "$sample_id" =~ (S|PF)$ ]]; then
        is_broad="YES"
    fi
    echo "$control_id $is_broad"
}

# ==============================================================================
# PHASE 1: TAG DIRECTORY CREATION
# ==============================================================================

echo -e "\n======================================================="
echo "PHASE 1: CREATING HOMER TAG DIRECTORIES"
echo "======================================================="

ALL_SAMPLES=("${TARGET_SAMPLES[@]}" "${CONTROL_SAMPLES[@]}")
printf "%s\n" "${ALL_SAMPLES[@]}" | sort -u > temp_samples.txt

while read sample_id; do
    local_bam="$filtered_bam_files/${sample_id}_mm10.filtered.sorted.bam"
    local_tag_dir="$tag_dirs/${sample_id}_tag"

    if [ ! -f "$local_bam" ]; then
        echo "WARNING: BAM not found for $sample_id. Skipping."
        continue
    fi

    if [ ! -d "$local_tag_dir" ]; then
        echo "  -> Creating Tag Directory for: $sample_id"
        makeTagDirectory "$local_tag_dir" "$local_bam" \
            -format sam -sspe -genome "$GENOME_VERSION"
    else
        echo "  -> Tag Dir exists for $sample_id. Skipping."
    fi
done < temp_samples.txt

rm temp_samples.txt

# ==============================================================================
# PHASE 2: PEAK CALLING (Sharp and Conditional Broad)
# ==============================================================================

echo -e "\n======================================================="
echo "PHASE 2: HOMER PEAK CALLING"
echo "======================================================="

run_homer_peaks() {
    local sample_id=$1
    local control_id=$2
    local is_broad=$3

    local treatment_tag="$tag_dirs/${sample_id}_tag"
    local control_tag="$tag_dirs/${control_id}_tag"

    if [ ! -d "$treatment_tag" ] || [ ! -d "$control_tag" ]; then
        echo "ERROR: Missing Tag Dir. Skipping $sample_id."
        return
    fi

    # ---- SHARP PEAKS (TFs/K4me3/Baseline - style factor defaults) ----
    local sharp_output="$peaks_out/${sample_id}_HOMER_SHARP.txt"
    echo "  -> HOMER SHARP Peaks (-style factor)"
    
    findPeaks "$treatment_tag" -i "$control_tag" -o "$sharp_output" \
        -style factor -p $cores -norm $NORMALIZATION_FACTOR

    # ---- BROAD PEAKS (Selective for S/PF factors - style histone) ----
    if [ "$is_broad" == "YES" ]; then
        local broad_output="$peaks_out/${sample_id}_HOMER_BROAD.txt"
        echo "  -> HOMER BROAD Peaks (-style histone)"
        
        findPeaks "$treatment_tag" -i "$control_tag" -o "$broad_output" \
            -style histone -size 1000 -p $cores -norm $NORMALIZATION_FACTOR
    fi
}

for sample in "${TARGET_SAMPLES[@]}"; do
    INFO=($(get_control_info "$sample"))
    CONTROL="${INFO[0]}"
    IS_BROAD="${INFO[1]}"

    if [ "$CONTROL" == "UNKNOWN_CONTROL" ]; then
        echo "Skipping $sample: Unknown control"
        continue
    fi

    echo -e "\n======================================================="
    echo "Processing Sample: $sample (Control: $CONTROL)"
    echo "Broad peaks required? $IS_BROAD"
    echo "======================================================="

    run_homer_peaks "$sample" "$CONTROL" "$IS_BROAD"
done

# ==============================================================================
# PHASE 3: MOTIF ENRICHMENT (ON SHARP PEAKS)
# ==============================================================================

echo -e "\n======================================================="
echo "PHASE 3: MOTIF ENRICHMENT"
echo "======================================================="

run_motif_analysis() {
    local sample_id=$1
    local peak_file="$peaks_out/${sample_id}_HOMER_SHARP.txt"
    local motif_dir="$motif_out/${sample_id}_Motif_Results"

    if [ ! -f "$peak_file" ]; then
        echo "ERROR: No sharp peak file for $sample_id. Skipping motifs."
        return
    fi

    echo "  -> Running Motif Enrichment for $sample_id"

    findMotifsGenome.pl "$peak_file" "$GENOME_VERSION" "$motif_dir" \
        -size 200 -mask -preparsed -p $cores
}

for sample in "${TARGET_SAMPLES[@]}"; do
    run_motif_analysis "$sample"
done

# ==============================================================================
# PHASE 4: ENHANCER / SUPER ENHANCER ANALYSIS
# ==============================================================================

echo -e "\n======================================================="
echo "PHASE 4: SUPER ENHANCER (SE) ANALYSIS (ROSE Implementation)"
echo "======================================================="

declare -a SE_TARGETS=("CR12CGR" "CR12DGR" "CR13CGR" "CR13DGR" "CR15CK4me3" "CR15DK4me3")

run_super_enhancers() {
    local sample_id=$1
    local control_id=$2
    
    local treatment_tag="$tag_dirs/${sample_id}_tag"
    local control_tag="$tag_dirs/${control_id}_tag"
    local se_results_dir="$se_out/${sample_id}_SE_Results"

    mkdir -p "$se_results_dir"

    echo -e "\n--- Processing Super Enhancers for $sample_id ---"
    
    # 1. Identify all candidate regulatory regions (CSRs - Enhancers)
    local csr_output="$se_results_dir/${sample_id}_CandidateRegions.txt"
    echo "  -> 1. Identifying Candidate Regulatory Regions (CSRs)"
    
    findPeaks "$treatment_tag" -i "$control_tag" -o "$csr_output" \
        -style super -L 0 -fdr 0.001 -minDist 500 -p $cores

    if [ ! -f "$csr_output" ]; then
        echo "WARNING: No Candidate Regulatory Regions found for $sample_id. Skipping ROSE."
        return
    fi

    # 2. Run findCSREs.pl (HOMER's implementation of ROSE)
    echo "  -> 2. Ranking and Defining Super Enhancers (ROSE method)"
    
    findCSREs.pl "$treatment_tag" -r "$csr_output" -o "$se_results_dir/${sample_id}_SE" -norm -p $cores
    
    echo "Super Enhancer results saved in: $se_results_dir"
}

for sample in "${SE_TARGETS[@]}"; do
    INFO=($(get_control_info "$sample"))
    CONTROL="${INFO[0]}"
    
    if [ "$CONTROL" == "UNKNOWN_CONTROL" ]; then
        echo "Skipping $sample: Unknown control for SE analysis."
        continue
    fi
    
    run_super_enhancers "$sample" "$CONTROL"
done

# ==============================================================================
# PHASE 5: PEAK ANNOTATION AND PAIRWISE OVERLAP
# ==============================================================================

echo -e "\n======================================================="
echo "PHASE 5: PEAK ANNOTATION AND PAIRWISE OVERLAP"
echo "======================================================="

## SECTION 5A: Individual Peak Annotation
echo -e "\n--- 5A: Annotating Individual Sharp Peak Sets ---"

run_annotation() {
    local sample_id=$1
    local peak_file="$peaks_out/${sample_id}_HOMER_SHARP.txt"
    local annotation_file="$annotation_out/${sample_id}_Annotation.txt"

    if [ ! -f "$peak_file" ]; then
        echo "WARNING: Sharp peak file not found for $sample_id. Skipping annotation."
        return
    fi

    echo "  -> Annotating $sample_id..."

    annotatePeaks.pl "$peak_file" "$GENOME_VERSION" > "$annotation_file" -p $cores
}

for sample in "${TARGET_SAMPLES[@]}"; do
    run_annotation "$sample"
done


## SECTION 5B: Pairwise Overlap Analysis (C vs. D)
echo -e "\n--- 5B: Pairwise Overlap Analysis (Control vs. Treated) ---"

# --- GR Overlap (CR12 and CR13) ---
for batch in CR12 CR13; do
    local C_sample="${batch}CGR"
    local D_sample="${batch}DGR"

    local C_peak_file="$peaks_out/${C_sample}_HOMER_SHARP.txt"
    local D_peak_file="$peaks_out/${D_sample}_HOMER_SHARP.txt"
    local output_file="$annotation_out/${batch}_GR_Overlap_Venn.txt"

    if [ ! -f "$C_peak_file" ] || [ ! -f "$D_peak_file" ]; then
        echo "WARNING: Skipping GR Overlap for $batch (Missing peak files)."
        continue
    fi

    echo "  -> Analyzing Overlap for ${C_sample} vs. ${D_sample}..."
    mergePeaks "$C_peak_file" "$D_peak_file" -d given -delim "," -venn "$output_file"
done

# --- H3K4me3 Overlap (CR15) ---
local C_sample="CR15CK4me3"
local D_sample="CR15DK4me3"

local C_peak_file="$peaks_out/${C_sample}_HOMER_SHARP.txt"
local D_peak_file="$peaks_out/${D_sample}_HOMER_SHARP.txt"
local output_file="$annotation_out/CR15_K4me3_Overlap_Venn.txt"

if [ ! -f "$C_peak_file" ] || [ ! -f "$D_peak_file" ]; then
    echo "WARNING: Skipping K4me3 Overlap for CR15 (Missing peak files)."
else
    echo "  -> Analyzing Overlap for ${C_sample} vs. ${D_sample}..."
    mergePeaks "$C_peak_file" "$D_peak_file" -d given -delim "," -venn "$output_file"
fi


# ==============================================================================
# PHASE 6: COREPRESSOR OVERLAP QUANTIFICATION (S/PF on GR Enhancers/SEs)
# ==============================================================================

echo -e "\n======================================================="
echo "PHASE 6: COREPRESSOR OVERLAP QUANTIFICATION (S/PF on GR Enhancer/SEs)"
echo "======================================================="

declare -a ANALYSIS_PAIRS=(
    "CR12 C GR S"   
    "CR12 D GR S"   
    "CR13 C GR S"   
    "CR13 D GR S"   
    "CR13 C GR PF"  
    "CR13 D GR PF"  
)

run_quantification() {
    local batch=$1
    local condition=$2
    local se_factor=$3
    local target_factor=$4
    local region_type=$5

    local se_sample="${batch}${condition}${se_factor}"
    local target_sample="${batch}${condition}${target_factor}"
    local control_id="${batch}${condition}I"

    if [ "$region_type" == "SE" ]; then
        local regions_file="$se_out/${se_sample}_SE_Results/${se_sample}_SE_SuperEnhancers.txt"
        local output_suffix="SE_Quant"
    else
        local regions_file="$se_out/${se_sample}_SE_Results/${se_sample}_CandidateRegions.txt"
        local output_suffix="Enhancer_Quant"
    fi

    local target_tag_dir="$tag_dirs/${target_sample}_tag"
    local control_tag_dir="$tag_dirs/${control_id}_tag"
    local output_file="$overlap_out/${se_sample}_vs_${target_factor}_${output_suffix}.txt"

    echo "  -> Quantifying ${target_sample} over ${se_sample} ${region_type}s..."

    if [ ! -f "$regions_file" ] || [ ! -d "$target_tag_dir" ] || [ ! -d "$control_tag_dir" ]; then
        echo "ERROR: Missing required input file/directory. Skipping ${region_type} quantification."
        return 1
    fi

    analyzeRegions.pl "$regions_file" "$GENOME_VERSION" \
        -d "$target_tag_dir" "$control_tag_dir" -norm $NORMALIZATION_FACTOR -raw \
        > "$output_file"
}

for pair_info in "${ANALYSIS_PAIRS[@]}"; do
    read -r batch condition se_factor target_factor <<< "$pair_info"

    run_quantification "$batch" "$condition" "$se_factor" "$target_factor" "SE"
    run_quantification "$batch" "$condition" "$se_factor" "$target_factor" "Enhancer"
done

# ==============================================================================
# PHASE 7: MULTI-FACTOR OVERLAP (S vs PF vs GR vs K4me3)
# ==============================================================================

echo -e "\n======================================================="
echo "PHASE 7: MULTI-FACTOR OVERLAP (S, PF, GR, K4me3)"
echo "======================================================="

run_multi_overlap() {
    local condition=$1 # C or D
    local output_name="${condition}_All_Factors_CR13_Overlap_Venn.txt"
    local output_file="$multi_overlap_out/$output_name"

    local gr_file="$peaks_out/CR13${condition}GR_HOMER_SHARP.txt"
    local s_file="$peaks_out/CR13${condition}S_HOMER_SHARP.txt"
    local pf_file="$peaks_out/CR13${condition}PF_HOMER_SHARP.txt"
    local k4me3_file="$peaks_out/CR15${condition}K4me3_HOMER_SHARP.txt"

    echo "  -> Analyzing Multi-Factor Overlap for Condition $condition (GR, S, PF, K4me3)..."

    if [ ! -f "$gr_file" ] || [ ! -f "$s_file" ] || [ ! -f "$pf_file" ] || [ ! -f "$k4me3_file" ]; then
        echo "WARNING: Skipping multi-factor overlap for Condition $condition (Missing one or more peak files)."
        return 1
    fi

    mergePeaks "$gr_file" "$s_file" "$pf_file" "$k4me3_file" \
        -d 200 -delim "," -venn "$output_file"

    echo "  -> Results saved to: $output_file"
}

run_multi_overlap "C"
run_multi_overlap "D"

# ==============================================================================
# PHASE 8: BIGWIG TRACK GENERATION FOR UCSC
# ==============================================================================

echo -e "\n======================================================="
echo "PHASE 8: CREATING BIGWIG FILES FOR GENOME BROWSER 📊"
echo "======================================================="

declare -a ALL_VIS_SAMPLES=("${TARGET_SAMPLES[@]}" "${CONTROL_SAMPLES[@]}")

run_makeUCSCfile() {
    local sample_id=$1
    local tag_dir="$tag_dirs/${sample_id}_tag"

    if [ ! -d "$tag_dir" ]; then
        echo "WARNING: Tag Dir not found for $sample_id. Skipping BigWig creation."
        return
    fi

    echo "  -> Generating BigWig for: $sample_id"

    makeUCSCfile "$tag_dir" -o auto -style factor -norm $NORMALIZATION_FACTOR -p $cores -fragLength auto

    echo "  -> File created in: ${tag_dirs}/${sample_id}_tag/${sample_id}_tag_bigwig.bw"
}

for sample in "${ALL_VIS_SAMPLES[@]}"; do
    run_makeUCSCfile "$sample"
done

# ==============================================================================
# PHASE 9: MULTIWIG TRACK HUB GENERATION
# ==============================================================================

echo -e "\n======================================================="
echo "PHASE 9: CREATING UCSC MULTIWIG TRACK HUB 🌐"
echo "======================================================="

run_multi_wig_hub() {
    # 1. Collect all generated BigWig files
    local bw_files=()
    for sample in "${ALL_VIS_SAMPLES[@]}"; do
        # Note: BigWig files are saved directly inside their respective tag directories
        local bw_file="${tag_dirs}/${sample}_tag/${sample}_tag_bigwig.bw" 
        if [ -f "$bw_file" ]; then
            bw_files+=("$bw_file")
        fi
    done

    if [ ${#bw_files[@]} -eq 0 ]; then
        echo "ERROR: No BigWig files found. Skipping Track Hub creation."
        return 1
    fi

    echo "  -> Found ${#bw_files[@]} BigWig files. Generating Track Hub..."

    # 2. Run makeMultiWigHub.pl
    # The output files are written to the local HUB_WEB_DIR
    makeMultiWigHub.pl "${bw_files[@]}" "$GENOME_VERSION" "Mouse" "CR Multi-Factor Analysis Hub" \
        -webDir "$HUB_WEB_DIR" \
        -url "$HUB_URL" \
        -norm $NORMALIZATION_FACTOR \
        -trackType "signal"

    echo -e "\n*** TRACK HUB CREATION COMPLETE ***"
    echo "The hub files have been written to the local directory: $HUB_WEB_DIR"
    echo "Remember to copy the contents of $HUB_WEB_DIR to your S3 bucket."
    echo -e "\nURL to share with UCSC: $HUB_URL/hub.txt\n"
}

run_multi_wig_hub

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

echo -e "\n\n======================================================="
echo "--- ALL HOMER ANALYSIS PHASES COMPLETE (1 through 9) ---"
echo "======================================================="
echo "Final Action Required: Copy the contents of the $HUB_WEB_DIR folder "
echo "to the root of your S3 bucket: $HUB_URL"
