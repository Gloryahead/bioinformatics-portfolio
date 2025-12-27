#!/bin/bash
# --- Configuration ---
cores=12

# Define base directory for the project
base_dir="/your/cluster/path/project_name"
echo "Project Base Directory: $base_dir"

# --- Genome and Annotation References ---
mm10_ref="/path/to/your/bowtie2/genome_index"      # Bowtie2 index
mm10_genome_dir="/path/to/your/genome"
genome_size="$mm10_genome_dir/mm10.chrom.sizes"
mm10_fa="$mm10_genome_dir/mm10.fa"
blacklist_bed="/path/to/your/blacklist_file/mm10-blacklist.v2.bed" # change mm10-blacklist.v2.bed to your preferred genome blacklist
GENOME_VERSION="mm10" # Mus musculus genome

# --- Temp Directory Setup ---
SAMTOOLS_TEMP_DIR="$base_dir/samtools_tmp"
PICARD_TEMP_DIR="$base_dir/picard_tmp"
DEEPTOOLS_TEMP_DIR="$base_dir/deeptools_tmp"

# --- Picard Setup ---
picard_url="https://github.com/broadinstitute/picard/releases/download/2.26.2/picard.jar"
picard_jar="$base_dir/picard.jar" 
picardCMD="java -Djava.io.tmpdir=$PICARD_TEMP_DIR -jar $picard_jar" 

if [ ! -f "$picard_jar" ]; then
    echo "Downloading Picard JAR..."
    wget -O "$picard_jar" "$picard_url" || { echo "ERROR: Failed to download Picard.jar"; exit 1; }
fi

# --- Project Subdirectories ---
fastq_files="$base_dir/fastq_files"
trimmed_reads="$base_dir/trimmed_fastq_files"
sam_mm10="$base_dir/sam_mm10"
bam_files="$base_dir/bam"
bed_files="$base_dir/bed"
sam_mm10_sum="$sam_mm10/sam_mm10_summary"
bedgraph="$base_dir/bedgraph"
fragmentLen="$base_dir/fragmentLen"
removeDuplicate="$bam_files/removeDuplicate"
picard_summary="$removeDuplicate/picard_summary"
bigwig_file="$base_dir/bigwig"
filtered_bam_files="$bam_files/filtered_blacklist"
macs2_out="$base_dir/MACS2"
homer_out="$base_dir/HOMER"
motif_out="$base_dir/Motifs"
visualization_dir="$base_dir/visualization"  # <--- NEW DIRECTORY

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
    "$sam_mm10" "$sam_mm10_sum" \
    "$bam_files" "$bed_files" "$bedgraph" "$fragmentLen" "$removeDuplicate" \
    "$picard_summary" "$bigwig_file" "$filtered_bam_files" \
    "$macs2_out" "$homer_out" "$motif_out" "$visualization_dir" \
    "$SAMTOOLS_TEMP_DIR" "$PICARD_TEMP_DIR" "$DEEPTOOLS_TEMP_DIR"

echo "Project directories verified."

# --- Sample List ---
declare -a samples=(
    "SRR29684800"
    "SRR29684801"
    "SRR29684802"
    "SRR29684803"
    "SRR29684804"
    "SRR29684805"
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
# --------------------------------------------------------------------------------
for sample_id in "${samples[@]}"; do
    echo -e "\n=== Processing sample: $sample_id ==="

    # Define inputs/outputs for this iteration
    fastq_r1="$fastq_files/${sample_id}_1.fastq.gz"
    fastq_r2="$fastq_files/${sample_id}_2.fastq.gz"
    trimmed_r1="$trimmed_reads/${sample_id}_trimmed_R1.fastq.gz"
    trimmed_r2="$trimmed_reads/${sample_id}_trimmed_R2.fastq.gz"
    
    # Alignment Outputs
    sam_output="$sam_mm10/${sample_id}_mm10.sam" 
    bam_sorted="$bam_files/${sample_id}_mm10.sorted.bam"
    
    # Deduplication Outputs
    dedup_bam="$removeDuplicate/${sample_id}_mm10.dedup.bam"
    dedup_metrics="$picard_summary/${sample_id}_picard.dupMark.txt"
    
    # Filtering Outputs
    filtered_bam="$filtered_bam_files/${sample_id}_mm10.final.bam"
    
    # BigWig Output (DeepTools)
    bw_output="$bigwig_file/${sample_id}.bw"

    # --- Step 1: Pre-check ---
    if [ ! -f "$fastq_r1" ] || [ ! -f "$fastq_r2" ]; then
        echo "WARNING: Raw FASTQ files for $sample_id not found. Skipping."
        continue
    fi

    # --- Step 2: Quality Control (FastQC) ---
    echo "--- Step 2: FastQC on raw reads ---"
    # fastqc -t $cores "$fastq_r1" "$fastq_r2" -o "$pre_trim_qc"

    # --- Step 3: Adapter Trimming (fastp) ---
    echo "--- Step 3: Adapter trimming (fastp) ---"
    if [ ! -f "$trimmed_r1" ]; then
        fastp -i "$fastq_r1" -o "$trimmed_r1" \
              -I "$fastq_r2" -O "$trimmed_r2" \
              --json "$trimmed_reads/${sample_id}_fastp.json" \
              --html "$trimmed_reads/${sample_id}_fastp.html" \
              -w "$cores"
    else
        echo "Trimmed files exist, skipping..."
    fi

    # --- Step 4: Alignment (Bowtie2) ---
    echo "--- Step 4: Aligning to mm10 ---"
    if [ ! -f "$sam_output" ] && [ ! -f "$bam_sorted" ]; then
        bowtie2 -q -x "$mm10_ref" -p "$cores" \
            -1 "$trimmed_r1" -2 "$trimmed_r2" \
            --sensitive --no-mixed --no-discordant --no-unal \
            -S "$sam_output"
    fi

    # --- Step 5: SAM to Sorted BAM ---
    echo "--- Step 5: Converting and Sorting BAM ---"
    if [ ! -f "$bam_sorted" ]; then
        samtools view -bS -@ "$cores" "$sam_output" | \
        samtools sort -@ "$cores" -T "$SAMTOOLS_TEMP_DIR/${sample_id}" -o "$bam_sorted"
        samtools index "$bam_sorted"
        rm "$sam_output" # Save space
    fi

    # --- Step 6: Mark Duplicates (Picard) ---
    echo "--- Step 6: Marking Duplicates ---"
    if [ ! -f "$dedup_bam" ]; then
        $picardCMD MarkDuplicates \
            I="$bam_sorted" \
            O="$dedup_bam" \
            M="$dedup_metrics" \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT
        samtools index "$dedup_bam"
    fi

    # --- Step 7: Filtering (Blacklist + Quality + Mito) ---
    echo "--- Step 7: Filtering alignments ---"
    if [ ! -f "$filtered_bam" ]; then
        temp_filtered="$filtered_bam_files/${sample_id}_temp.bam"
        
        samtools view -h -q 30 -f 2 -F 1804 "$dedup_bam" | \
        grep -v "chrM" | samtools view -b - > "$temp_filtered"

        bedtools intersect -v -abam "$temp_filtered" -b "$blacklist_bed" > "$filtered_bam"
        samtools index "$filtered_bam"
        
        echo "Alignment Stats for $sample_id:" > "$sam_mm10_sum/${sample_id}_stats.txt"
        samtools flagstat "$filtered_bam" >> "$sam_mm10_sum/${sample_id}_stats.txt"
        
        rm "$temp_filtered"
    fi

    # --- Step 8: BigWig Generation (deepTools) ---
    echo "--- Step 8: Generating BigWig tracks (deepTools) ---"
    if [ ! -f "$bw_output" ]; then
        bamCoverage --bam "$filtered_bam" \
            --outFileName "$bw_output" \
            --outFileFormat bigwig \
            --binSize 10 \
            --normalizeUsing RPKM \
            --extendReads \
            --ignoreDuplicates \
            --numberOfProcessors "$cores"
    fi

    # --- Step 9: Peak Calling ---
    
    # 9a. MACS2
    echo "--- Step 9a: Calling Peaks (MACS2) ---"
    macs2_sample_dir="$macs2_out/${sample_id}"
    mkdir -p "$macs2_sample_dir"
    
    if [ ! -f "$macs2_sample_dir/${sample_id}_peaks.narrowPeak" ]; then
        macs2 callpeak \
            -t "$filtered_bam" \
            -f BAMPE \
            -g mm \
            -n "${sample_id}" \
            --outdir "$macs2_sample_dir" \
            -q 0.01 \
            --nomodel \
            --keep-dup all \
            --call-summits
    fi

    # 9b. HOMER Peak Calling & Annotation
    echo "--- Step 9b: Calling & Annotating Peaks (HOMER) ---"
    homer_sample_dir="$homer_out/${sample_id}"
    homer_tag_dir="$homer_sample_dir/tag_dir"
    mkdir -p "$homer_sample_dir"
    mkdir -p "$homer_sample_dir/go"
    mkdir -p "$homer_sample_dir/genomeOntology"
    
    if [ ! -f "$homer_sample_dir/peaks_${sample_id}_fdr1e-05_annotated.txt" ]; then
        # Create Tag Directory
        if [ ! -d "$homer_tag_dir" ]; then
            makeTagDirectory "$homer_tag_dir" \
                "$filtered_bam" \
                -genome "$GENOME_VERSION" \
                -checkGC \
                -tbp 1 \
                -format sam 
        fi
        
        # Find Peaks
        findPeaks "$homer_tag_dir" \
            -style dnase \
            -o "$homer_sample_dir/peaks_${sample_id}_fdr1e-05.txt" \
            -fdr 1e-05 \
            -region \
            -minDist 150 \
            -size 150
            
        # Annotate Peaks
        annotatePeaks.pl "$homer_sample_dir/peaks_${sample_id}_fdr1e-05.txt" "$GENOME_VERSION" \
            -d "$homer_tag_dir" \
            -go "$homer_sample_dir/go" \
            -genomeOntology "$homer_sample_dir/genomeOntology" \
            > "$homer_sample_dir/peaks_${sample_id}_fdr1e-05_annotated.txt"
    fi

    # --- Step 10: Motif Finding (Using MACS2 peaks) ---
    echo "--- Step 10: Motif Finding (HOMER using MACS2 peaks) ---"
    motif_sample_dir="$motif_out/${sample_id}_motifs"
    input_peak_file="$macs2_sample_dir/${sample_id}_peaks.narrowPeak"
    
    if [ -f "$input_peak_file" ]; then
        mkdir -p "$motif_sample_dir"
        # Only run if directory is empty to avoid overwriting/re-running
        if [ -z "$(ls -A $motif_sample_dir)" ]; then
            findMotifsGenome.pl "$input_peak_file" "$GENOME_VERSION" \
                "$motif_sample_dir" \
                -size 200 \
                -len 8,10,12,15 \
                -mask \
                -p "$cores"
        fi
    fi

    # --- Step 11: Create Visualization Files (Tracks & BEDs) ---
    echo "--- Step 11: Creating Visualization Files ---"
    
    # 11a. HOMER BigWig (Alternative to deepTools)
    # Note: Using $genome_size which we generated from mm10.fa
    # Outputting to the HOMER sample directory by default, then we can link/copy if needed
    if [ -d "$homer_tag_dir" ]; then
         makeUCSCfile "$homer_tag_dir" \
             -o auto \
             -bigWig "$genome_size" \
             -style dnase
    fi

    # 11b. Convert HOMER peaks to standard BED for sharing
    homer_peak_file="$homer_sample_dir/peaks_${sample_id}_fdr1e-05.txt"
    if [ -f "$homer_peak_file" ]; then
        pos2bed.pl "$homer_peak_file" > "$visualization_dir/${sample_id}_homer_peaks.bed"
    fi

    # 11c. Copy MACS2 peaks to visualization folder
    macs2_peak_file="$macs2_sample_dir/${sample_id}_peaks.narrowPeak"
    if [ -f "$macs2_peak_file" ]; then
        cp "$macs2_peak_file" "$visualization_dir/${sample_id}_macs2_peaks.bed"
    fi

    echo "=== Finished processing $sample_id ==="
done

echo "All samples processed successfully."
