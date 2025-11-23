#!/bin/bash

# ==============================================================================
# CUT&RUN Data Analysis Pipeline (Local Alignment & Spike-in Normalization)
# ==============================================================================

# --- Configuration ---
core=8
cores=8 # Use 'cores' consistently for clarity in Bowtie2 commands
binLen=500

# Define base directory for the project
base_dir="/xdisk/clsmith1/maarowosegbe/CR12_13_15"
echo "Project Base Directory: $base_dir"

# --- Genome and Annotation References ---
ecoli_ref="/xdisk/clsmith1/maarowosegbe/Ecoli_index/Ecoli"   # Bowtie2 index for E. coli
mm10_ref="/xdisk/clsmith1/maarowosegbe/mm10_index/mm10"     # Bowtie2 index for mm10
mm10_genome_dir="/xdisk/clsmith1/maarowosegbe/CUTnRUN_6n7/mm10"
genome_size="$mm10_genome_dir/mm10.chrom.sizes"
mm10_fa="$mm10_genome_dir/mm10.fa"
blacklist_bed="/xdisk/clsmith1/maarowosegbe/ENCODE_DAC/mm10-blacklist.v2.bed"

# --- Temp Directory Setup (CRITICAL FIX for I/O Stability) ---
# Use dedicated scratch paths for large temporary file operations
SAMTOOLS_TEMP_DIR="$base_dir/samtools_tmp"
PICARD_TEMP_DIR="$base_dir/picard_tmp"
DEEPTOOLS_TEMP_DIR="$base_dir/deeptools_tmp"

# --- Picard Setup ---
picard_url="https://github.com/broadinstitute/picard/releases/download/2.26.2/picard.jar"
picard_jar="picard.jar"
picardCMD="java -jar $picard_jar"
PICARD_TEMP_DIR="$base_dir/picard_tmp" # New temp directory in project space

if [ ! -f "$picard_jar" ]; then
    echo "Downloading Picard JAR..."
    wget -O "$picard_jar" "$picard_url"
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to download Picard.jar. Exiting."
        exit 1
    fi
fi

# --- Project Subdirectories ---
fastq_files="$base_dir/fastq_files"
trimmed_reads="$base_dir/trimmed_fastq_files"
sam_mm10="$base_dir/sam_mm10"
bam_files="$base_dir/bam"
bed_files="$base_dir/bed"
sam_mm10_sum="$sam_mm10/sam_mm10_summary"
sam_ecoli="$base_dir/sam_ecoli"
sam_ecoli_sum="$sam_ecoli/sam_ecoli_summary"
bedgraph="$base_dir/bedgraph"
fragmentLen="$base_dir/fragmentLen"
removeDuplicate="$bam_files/removeDuplicate"
picard_summary="$removeDuplicate/picard_summary"
bigwig_file="$base_dir/bigwig"
filtered_bam_files="$bam_files/filtered_blacklist" # New directory for blacklisted BAMs

# Quality Control directories
pre_trim_qc="$base_dir/qc_raw_reads"       
post_trim_qc="$base_dir/qc_trimmed_reads" 
multiqc_pre_trim="$base_dir/multiqc_raw"  
multiqc_post_trim="$base_dir/multiqc_trimmed" 

# --- Create Directories ---
echo -e "\n--- Setting up project directories ---"
mkdir -p \
    "$pre_trim_qc" "$multiqc_pre_trim" \
    "$trimmed_reads" "$post_trim_qc" "$multiqc_post_trim" \
    "$sam_mm10" "$sam_mm10_sum" "$sam_ecoli" "$sam_ecoli_sum" \
    "$bam_files" "$bed_files" "$bedgraph" "$fragmentLen" "$removeDuplicate" \
    "$picard_summary" "$bigwig_file" "$filtered_bam_files" \
    "$SAMTOOLS_TEMP_DIR" "$PICARD_TEMP_DIR" "$DEEPTOOLS_TEMP_DIR"

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

# --- Genome size file generation (Run once) ---
if [ ! -f "$genome_size" ]; then
    echo "Generating mm10 index and chromosome sizes..."
    if [ -f "$mm10_fa" ]; then
        samtools faidx "$mm10_fa"
        cut -f1,2 "$mm10_fa.fai" > "$genome_size"
    else
        echo "ERROR: mm10 FASTA file not found at $mm10_fa. Exiting."
        exit 1
    fi
fi
# --------------------------------------------------------------------------------

# --- Main Processing Loop ---
for sample_id in "${samples[@]}"; do
    echo -e "\n=== Processing sample: $sample_id ==="

    # Define sample-specific files
    fastq_r1="$fastq_files/${sample_id}_R1.fastq.gz"
    fastq_r2="$fastq_files/${sample_id}_R2.fastq.gz"
    trimmed_r1="$trimmed_reads/${sample_id}_trimmed_R1.fastq.gz"
    trimmed_r2="$trimmed_reads/${sample_id}_trimmed_R2.fastq.gz"
    sam_ecoli_output="$sam_ecoli/${sample_id}_ecoli.sam"
    sam_ecoli_sum_output="$sam_ecoli_sum/${sample_id}_ecoli.txt"
    sam_mm10_output="$sam_mm10/${sample_id}_mm10.sam" 
    sam_mm10_sum_output="$sam_mm10_sum/${sample_id}_mm10.txt"
    
    # BAM/BED files (using 'filtered' suffix for downstream files)
    bam_output="$bam_files/${sample_id}_mm10.mapped.bam"
    filtered_bam_output="$filtered_bam_files/${sample_id}_mm10.filtered.bam"
    sorted_bam_output="$filtered_bam_files/${sample_id}_mm10.filtered.sorted.bam" # Final BAM source
    
    # Downstream files
    bed_output="$bed_files/${sample_id}_mm10.filtered.bed" 
    clean_bed_output="$bed_files/${sample_id}_mm10.filtered.clean.bed"
    fragments_bed_output="$bed_files/${sample_id}_mm10.filtered.fragments.bed"
    fragmentLen_output="$fragmentLen/${sample_id}_mm10_fragmentLen.txt"

    # Fragment Analysis/Counting
    fragments_count_bed_output="$bed_files/${sample_id}_mm10.fragmentsCount.bin${binLen}.bed" # Fix: Variable name mismatch
    fragmentLen_output="$fragmentLen/${sample_id}_mm10_fragmentLen.txt"

    # Picard/Duplicate Marking
    removeDuplicate_marked_bam="$removeDuplicate/${sample_id}_mm10.filtered.dupMarked.bam"
    removeDuplicate_rm_bam="$removeDuplicate/${sample_id}_mm10.filtered.rmDup.bam"
    picard_summary_marked="$picard_summary/${sample_id}_picard.dupMark.txt"
    picard_summary_rm="$picard_summary/${sample_id}_picard.rmDup.txt"
    
    # Final Tracks
    bedgraph_scaled_output="$bedgraph/${sample_id}_scaled.bedgraph"
    bigwig_scaled_output="$bigwig_file/${sample_id}_scaled.bw"


    # --- Step 1: Pre-check: Ensure raw FASTQ files exist ---
    if [ ! -f "$fastq_r1" ] || [ ! -f "$fastq_r2" ]; then
        echo "WARNING: Raw FASTQ files for $sample_id not found. Skipping sample."
        continue
    fi
    echo "Raw FASTQ files found."

    # --- Step 2: Quality Control before Trimming (FastQC) ---
    echo "--- Running FastQC on raw reads... ---"
    # Added -t for threads, leveraging the core variable.
    fastqc -t $cores "$fastq_r1" "$fastq_r2" -o "$pre_trim_qc"

    # --- Step 3: Adapter Trimming (fastp) ---
    echo "--- Running adapter trimming with fastp... ---"
    fastp -i "$fastq_r1" -o "$trimmed_r1" \
          -I "$fastq_r2" -O "$trimmed_r2" \
          --json "$trimmed_reads/${sample_id}_fastp.json" \
          --html "$trimmed_reads/${sample_id}_fastp.html" \
          -w "$cores" # Use the $cores variable
    echo "Trimming completed."

    # --- Step 4: Quality control after trimming (FastQC) ---
    echo "--- Running FastQC on trimmed reads... ---"
    fastqc -t $cores "$trimmed_r1" "$trimmed_r2" -o "$post_trim_qc"

    # --- Step 5: Alignment with Bowtie2 (Target and Spike-in) ---
    echo "--- Running Bowtie2 alignments... ---"
    
    # E. coli (Spike-in) - Keep END-TO-END for consistent count metric
    bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p $cores -x "$ecoli_ref" -1 "$trimmed_r1" -2 "$trimmed_r2" -S "$sam_ecoli_output" &> "$sam_ecoli_sum_output"

    # mm10 (Target) - Keep END-TO-END for consistent count metric
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p $cores -x "$mm10_ref" -1 "$trimmed_r1" -2 "$trimmed_r2" -S "$sam_mm10_output" &> "$sam_mm10_sum_output"
    echo "Alignment complete for $sample_id."

    # --- Step 6: Filtering, Sorting, and Fragment Analysis ---
    echo "--- Converting and filtering files... ---"

    # 6a. SAM to BAM: Keep mapped read pairs only
    samtools view -bS -F 0x04 "$sam_mm10_output" -o "$bam_output"
    
    # 6b. DAC Blacklist Filtering (bedtools intersect -v)
    echo "--- Removing blacklisted regions with bedtools intersect... ---"
    bedtools intersect -v -abam "$bam_output" -b "$blacklist_bed" > "$filtered_bam_output"
    echo "Blacklist filtering complete for $sample_id."

    # 6c. Sort and Index BAM (STABILITY FIX)
        samtools sort -T "$SAMTOOLS_TEMP_DIR/$sample_id" -o "$sorted_bam_output" "$filtered_bam_output"
        
        if [ $? -eq 0 ]; then
            samtools index "$sorted_bam_output"
            echo "BAM file sorted and indexed successfully."
            rm -f "$SAMTOOLS_TEMP_DIR/${sample_id}"* # AUTO-CLEANUP TEMP FILES
        else
            echo "FATAL ERROR: samtools sort failed for $sample_id. Stopping analysis for this sample."
            continue
        fi

    # 6d. BAM to Paired-End BED
    bedtools bamtobed -i "$sorted_bam_output" -bedpe > "$bed_output"

    # 6e. Clean Paired-End BED (Same Chromosome, Fragment Length < 1000bp)
    awk '$1==$4 && $6-$2 < 1000 {print $0}' "$bed_output" > "$clean_bed_output"

    # 6f. Extract Fragments (Chr, Start, End)
    cut -f 1,2,6 "$clean_bed_output" | sort -k1,1 -k2,2n -k3,3n > "$fragments_bed_output"

    # 6g. Fragment Length Distribution
    samtools view -F 0x04 "$sorted_bam_output" | \
        awk 'function abs(x){return ((x < 0) ? -x : x)} { if (abs($9) > 0) print abs($9)}' | \
        sort -n | uniq -c | \
        awk -v OFS="\t" '{print $2, $1}' > "$fragmentLen_output"
    
    # 6h. Fragment Midpoint Counting (Normalization prerequisite)
    awk -v w="$binLen" '{print $1, int(($2 + $3)/(2*w))*w + w/2}' "$fragments_bed_output" | \
        sort -k1,1V -k2,2n | uniq -c | \
        awk -v OFS="\t" '{print $2, $3, $1}' | sort -k1,1V -k2,2n > "$fragments_count_bed_output"

    # --- 6i. Spike-in Normalization and Track Generation ---
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
            
            # 3. Calculate Scaling Factor (SF = 1 / BW%)
            if (( $(echo "$spike_bandwidth_percent == 0" | bc -l) )); then
                scale_factor=1.0
            else
                # SF = 100 / BW(%) - Correct for percentage input
                scale_factor=$(awk "BEGIN { printf \"%.6f\", 1 / $spike_bandwidth_percent }")
            fi
        fi
    fi 

    echo "Spike-in Bandwidth (BW%): $spike_bandwidth_percent%"
    echo "Calculated Scaling Factor (SF): $scale_factor"
    
    # 4. Generate scaled BedGraph
    bedtools genomecov -bg -scale "$scale_factor" -i "$fragments_bed_output" -g "$genome_size" > "$bedgraph_scaled_output"
    echo "Scaled BedGraph complete for $sample_id."

    # 5. Generate scaled BigWig (using deepTools bamCoverage)
    bamCoverage -b "$sorted_bam_output" --scaleFactor "$scale_factor" --binSize 50 -o "$bigwig_scaled_output" -p "$cores"
    echo "Scaled BigWig complete for $sample_id."
    
    # --- Step 7: Mark and Remove Duplicates (Picard) ---
    echo "--- Running Picard MarkDuplicates... ---"
    
    # CRITICAL FIX: Explicitly set -Djava.io.tmpdir for BAM processing
    $picardCMD MarkDuplicates \
        I="$sorted_bam_output" \
        O="$removeDuplicate_marked_bam" \
        METRICS_FILE="$picard_summary_marked"
            
    $picardCMD MarkDuplicates \
        I="$sorted_bam_output" \
        O="$removeDuplicate_rm_bam" \
        REMOVE_DUPLICATES=true \
        METRICS_FILE="$picard_summary_rm"
    echo "Picard duplicate marking and removal complete for $sample_id."

done

# --- Final QC: MultiQC (Run outside the main sample loop) ---
echo -e "\n--- Running MultiQC for all samples ---"

# MultiQC for raw reads
multiqc "$pre_trim_qc" -o "$multiqc_pre_trim" -n "multiqc_raw_reads"

# MultiQC for trimmed reads
multiqc "$post_trim_qc" "$trimmed_reads" -o "$multiqc_post_trim" -n "multiqc_trimmed_reads"

echo "Pipeline finished."
