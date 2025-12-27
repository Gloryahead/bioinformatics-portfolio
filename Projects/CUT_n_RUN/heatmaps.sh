#!/bin/bash

# ==============================================================================
# CUT&RUN Deeptools Visualization Script (FINAL, CORRECTED, HPC SAFE)
# ==============================================================================

set -euo pipefail

cores=8
BW_SUFFIX="_RPKM_10bp.bw"

# --- Configuration & Paths ---
base_dir="/xdisk/clsmith1/maarowosegbe/CR12_13_15"
peak_dir="$base_dir/peak_calling_combined/MACS2"
transcript_file="/xdisk/clsmith1/maarowosegbe/CUTnRUN_6n7/mm10"
heatmap_file="$base_dir/heatmaps"
bigwig_dir="$base_dir/bigwig/bigwig_RPKM_10bp"

mkdir -p "$heatmap_file"

# --- Sample Definitions ---
declare -a ALL_SAMPLES=(
    "CR12CGR" "CR12CI" "CR12CS" "CR12DGR" "CR12DI" "CR12DS"
    "CR13CGR" "CR13CI" "CR13CPF" "CR13CS" "CR13DGR" "CR13DI"
    "CR13DPF" "CR13DS"
    "CR15CI" "CR15CK4me3" "CR15CPF" "CR15CS"
    "CR15DI" "CR15DK4me3" "CR15DPF" "CR15DS"
)

declare -a MACS2_Q_LEVELS=("Q01" "Q005" "Q001")

# Peak-containing samples
declare -a TARGET_PEAK_SAMPLES=(
    "CR12CGR" "CR12CS" "CR12DGR" "CR12DS"
    "CR13CGR" "CR13CPF" "CR13CS" "CR13DGR" "CR13DPF" "CR13DS"
    "CR15CK4me3" "CR15CPF" "CR15CS" "CR15DK4me3" "CR15DPF" "CR15DS"
)

# -------------------------
# MAP ALL BIGWIG FILES
# -------------------------
declare -A ALL_BW_FILES=()
declare -a BW_ARGS=()

for S in "${ALL_SAMPLES[@]}"; do
    BW="$bigwig_dir/${S}${BW_SUFFIX}"
    if [[ ! -f "$BW" ]]; then
        echo "ERROR: Missing BigWig file: $BW" >&2
        exit 1
    fi
    ALL_BW_FILES[$S]="$BW"
    BW_ARGS+=("$BW")
done


# -------------------------
# DIFFERENTIAL PEAK PAIRS
# -------------------------
declare -A DIFFERENTIAL_PAIRS=()

for S in "${TARGET_PEAK_SAMPLES[@]}"; do
    COND="${S:3:1}"   # C or D
    if [[ "$COND" == "C" ]]; then
        D="${S/C/D}"
        [[ " ${TARGET_PEAK_SAMPLES[*]} " =~ " $D " ]] && DIFFERENTIAL_PAIRS[$S]="$D"
    fi
done

# -------------------------
# CROSS-FACTOR PAIRS
# -------------------------
declare -A CROSS_FACTOR_PAIRS=(
    ["CR12CS"]="CR12CGR" ["CR12DS"]="CR12DGR"
    ["CR13CS"]="CR13CGR" ["CR13DS"]="CR13DGR"
    ["CR15CS"]="CR15CK4me3" ["CR15DS"]="CR15DK4me3"
    ["CR13CPF"]="CR13CS" ["CR13DPF"]="CR13DS"
    ["CR15CPF"]="CR15CS" ["CR15DPF"]="CR15DS"
    ["CR15CK4me3"]="CR15CPF" ["CR15DK4me3"]="CR15DPF"
)

REGIONS_GTF="$transcript_file/filtered_transcripts.gtf"
MATRIX_BODY="$heatmap_file/body_body.matrix.gz"
MATRIX_TSS="$heatmap_file/tss_ref_point.matrix.gz"


# ==============================================================================
# PART 1 – Convert MACS2 peaks → BED
# ==============================================================================
echo "Converting MACS2 peaks → BED"

for S in "${TARGET_PEAK_SAMPLES[@]}"; do
    for Q in "${MACS2_Q_LEVELS[@]}"; do
        NP="$peak_dir/${S}_MACS2_${Q}_SHARP_peaks.narrowPeak"
        BED="$peak_dir/${S}_MACS2_${Q}.bed"
        if [[ -f "$NP" && ! -f "$BED" ]]; then
            awk 'OFS="\t" {print $1,$2,$3,"'$S'_Peak_"NR,$8}' "$NP" > "$BED"
        fi
    done
done


# ==============================================================================
# PART 2 – MASTER MATRICES (CORRECT ARRAY PASSING)
# ==============================================================================
echo "Computing Master matrices"

if [[ ! -f "$MATRIX_BODY" ]]; then
    echo "Recomputing MATRIX_BODY..."
    computeMatrix scale-regions \
        -S "${BW_ARGS[@]}" \
        -R "$REGIONS_GTF" \
        --beforeRegionStartLength 3000 \
        --regionBodyLength 5000 \
        --afterRegionStartLength 3000 \
        --skipZeros \
        -p $cores \
        -o "$MATRIX_BODY"
fi

if [[ ! -f "$MATRIX_TSS" ]]; then
    echo "Recomputing MATRIX_TSS..."
    computeMatrix reference-point \
        -S "${BW_ARGS[@]}" \
        -R "$REGIONS_GTF" \
        --referencePoint TSS \
        --beforeRegionStartLength 1500 \
        --afterRegionStartLength 1500 \
        --skipZeros \
        -p $cores \
        -o "$MATRIX_TSS"
fi


# ==============================================================================
# PART 3 – DIFFERENTIAL PEAK MATRICES
# ==============================================================================
echo "Computing Differential Peak Matrices"

for C in "${!DIFFERENTIAL_PAIRS[@]}"; do
    D="${DIFFERENTIAL_PAIRS[$C]}"
    Q="Q005"

    C_BED="$peak_dir/${C}_MACS2_${Q}.bed"
    D_BED="$peak_dir/${D}_MACS2_${Q}.bed"
    INT="$peak_dir/${C}_${D}_${Q}_Intersection.bed"
    OUT="$heatmap_file/${C}_vs_${D}_${Q}_IntersectionPeak.matrix.gz"

    if [[ -f "$C_BED" && -f "$D_BED" && ! -f "$OUT" ]]; then

        if [[ ! -f "$INT" ]]; then
            echo "Creating Intersection Peak Set for $C vs $D"
            cat "$C_BED" "$D_BED" \
                | sort -k1,1 -k2,2n \
                | bedtools intersect -a stdin -b "$D_BED" -u -wa \
                | cut -f1-3 \
                | bedtools merge -i stdin \
                > "$INT"
        fi

        computeMatrix reference-point \
            --referencePoint center \
            --beforeRegionStartLength 2000 \
            --afterRegionStartLength 2000 \
            --binSize 10 \
            --regionsFileName "$INT" \
            --scoreFileName "${ALL_BW_FILES[$C]}" "${ALL_BW_FILES[$D]}" \
            -p $cores \
            -o "$OUT"
    fi
done


# ==============================================================================
# PLOTTING FUNCTION (SAFE)
# ==============================================================================
plot_subset() {
    local matrix="$1"
    local title="$2"
    local A="$3"
    local B="$4"
    local prefix="$5"

    computeMatrixOperations subset \
        -m "$matrix" \
        --samples "${ALL_BW_FILES[$A]}" "${ALL_BW_FILES[$B]}" \
        -o "$heatmap_file/${prefix}.subset.matrix.gz"

    plotHeatmap \
        -m "$heatmap_file/${prefix}.subset.matrix.gz" \
        -out "$heatmap_file/${prefix}.png" \
        --samplesLabel "$A" "$B" \
        --sortUsing sum \
        --sortRegions keep \
        --colorMap viridis \
        --plotTitle "$title"

    plotProfile \
        -m "$heatmap_file/${prefix}.subset.matrix.gz" \
        -out "$heatmap_file/${prefix}.profile.png" \
        --samplesLabel "$A" "$B" \
        --plotTitle "$title"
}


# ==============================================================================
# PART 4 – GENE BODY / TSS COMPARISONS
# ==============================================================================
echo "Generating Gene/TSS plots"

for A in "${!CROSS_FACTOR_PAIRS[@]}"; do
    B="${CROSS_FACTOR_PAIRS[$A]}"

    keyA="${A:4}"
    keyB="${B:4}"
    batch="${A:0:4}"

    plot_subset "$MATRIX_BODY" "Gene Body" "$A" "$B" \
        "${batch}_${keyA}_vs_${keyB}_Body"

    plot_subset "$MATRIX_TSS" "TSS" "$A" "$B" \
        "${batch}_${keyA}_vs_${keyB}_TSS"
done


# ==============================================================================
# PART 5 – QC ON INDIVIDUAL PEAK FILES
# ==============================================================================
echo "Generating QC Plots"

for S in "${TARGET_PEAK_SAMPLES[@]}"; do
    for Q in "${MACS2_Q_LEVELS[@]}"; do

        batch="${S:0:4}"
        cond="${S:3:1}"

        if [[ "$cond" == "C" ]]; then
            CONTROL="${batch}CI"
        else
            CONTROL="${batch}DI"
        fi

        BED="$peak_dir/${S}_MACS2_${Q}.bed"
        QC="$heatmap_file/${S}_vs_${CONTROL}_${Q}.matrix.gz"

        if [[ -f "$BED" && ! -f "$QC" ]]; then
            echo "QC: $S vs $CONTROL ($Q)"
            computeMatrix reference-point \
                --referencePoint center \
                --beforeRegionStartLength 2000 \
                --afterRegionStartLength 2000 \
                --binSize 10 \
                --regionsFileName "$BED" \
                --scoreFileName "${ALL_BW_FILES[$S]}" "${ALL_BW_FILES[$CONTROL]}" \
                -p $cores \
                -o "$QC"
        fi

        [[ -f "$QC" ]] && \
            plot_subset "$QC" "Peak QC" "$S" "$CONTROL" \
                "${S}_vs_${CONTROL}_${Q}_QC"
    done
done

echo "All matrices and plots generated successfully."

