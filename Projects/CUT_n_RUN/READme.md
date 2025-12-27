# CUT&RUN Analysis Pipeline: Spike-in Normalization & Advanced Downstream Analysis

![Bash](https://img.shields.io/badge/Language-Bash-4EAA25)
![Environment](https://img.shields.io/badge/Manager-Conda%20%7C%20Bioconda-green)
![Normalization](https://img.shields.io/badge/Method-Spike--in%20Normalized-blue)

## 📌 Overview
This repository contains a modular, high-performance computing (HPC) pipeline for processing **CUT&RUN** chromatin profiling data. The workflow is designed for reproducibility and quantitative accuracy, featuring **E. coli spike-in normalization** to account for global shifts in chromatin accessibility.

Beyond standard alignment and peak calling, this pipeline includes advanced downstream modules for **Super Enhancer identification (ROSE logic)**, **Motif enrichment**, and **automated UCSC Track Hub generation**.



## 🧬 Pipeline Architecture
The analysis is split into modular Bash scripts optimized for SLURM job schedulers:

1.  **Preprocessing & Alignment (`run_bowtie...sh`)**
    * **QC:** `FastQC` & `MultiQC`.
    * **Trimming:** `fastp` (adapter & quality trimming).
    * **Alignment:** Dual alignment to **Target Genome (mm10)** and **Spike-in Genome (E. coli)** using `Bowtie2`.
    * **Filtering:** Removal of mitochondrial reads, unmapped reads, and **ENCODE Blacklist regions**.
    * **Normalization:** Calculation of scaling factors based on *E. coli* read depth (Formula: `SF = 1 / Spike_Bandwidth%`).

2.  **Peak Calling**
    * **`run_macs2_peaks.sh`:** Adaptive calling strategies.
        * **Sharp Mode:** For Transcription Factors (TFs).
        * **Broad Mode:** For Histone Marks (e.g., H3K27me3) and Corepressors.
    * **`homer.sh`:** Secondary peak calling using `findPeaks` (Factor vs. Histone styles) for validation.

3.  **Visualization (`bigwig.sh`)**
    * Generates `bigWig` tracks at multiple resolutions (10bp, 25bp, 50bp).
    * Produces both **RPKM-normalized** and **Spike-in scaled** tracks.

4.  **Downstream Analysis (`homer.sh`)**
    * **Motif Discovery:** `findMotifsGenome.pl` on sharp peaks.
    * **Super Enhancers:** Identification of large regulatory elements using ROSE-like algorithms.
    * **Overlap Analysis:** Quantification of factor binding at defined genomic regions.
    * **UCSC Hub:** Auto-generates structure for S3-hosted Track Hubs.

## 🛠 Installation & Dependencies

This project is part of a larger portfolio. To download and set up this specific pipeline:

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/Gloryahead/bioinformatics-portfolio.git](https://github.com/Gloryahead/bioinformatics-portfolio.git)
   cd bioinformatics-portfolio/Projects/CUT_n_RUN

## 📂 Directory Output Structure
The pipeline organizes inputs and outputs into a strict hierarchy to ensure reproducibility:

Project_Root/
├── bam/
│   ├── filtered_blacklist/  # Final analysis-ready BAMs
├── bigwig/                  # Scaled .bw tracks
├── MACS2/                   # .narrowPeak and .broadPeak files
├── HOMER/
│   ├── motifs/              # De novo motif results
│   ├── super_enhancers/     # ROSE/SE output
│   └── UCSC_Hub_Files/      # Ready for S3 upload
└── qc_reports/              # MultiQC HTML reports
