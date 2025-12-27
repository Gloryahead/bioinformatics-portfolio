# Automated RNA-seq Analysis Pipeline (STAR & FeatureCounts)

![Bash](https://img.shields.io/badge/Language-Bash-4EAA25)
![Aligner](https://img.shields.io/badge/Aligner-STAR-blue)
![Quantification](https://img.shields.io/badge/Quant-FeatureCounts-orange)

## 📌 Overview
This repository contains a robust, end-to-end bioinformatics pipeline for **bulk RNA-seq analysis**. Designed for High-Performance Computing (HPC) environments, it automates the workflow from raw FASTQ files to gene-level quantification, producing a final count matrix ready for differential expression analysis (e.g., DESeq2).

The pipeline features **automated reference management** (downloading annotations/building indices if missing) and **IGV-ready outputs** to streamline downstream visualization.



## 🧬 Pipeline Workflow
1.  **Reference Setup:** Checks for and builds STAR genome indices and downloads Gencode annotations if not present.
2.  **Quality Control:** * Raw Reads: `FastQC` + `MultiQC` aggregation.
    * Trimmed Reads: Post-trimming verification.
3.  **Trimming:** Adapter detection and removal using `fastp`.
4.  **Alignment:** Splicing-aware mapping to the reference genome (mm10/GRCm38) using `STAR`.
5.  **Visualization Prep:** Automatically sorts, indexes, and organizes BAM files for the Integrative Genomics Viewer (IGV).
6.  **Quantification:** Generates a combined gene counts matrix using `Subread featureCounts`.

## 📂 Directory Output Structure
The pipeline organizes inputs and outputs into a strict hierarchy to ensure reproducibility:

```text
MAA_science/
└── RNA_seq_human/                  # Project Base Directory
    ├── fastq_files/                # Raw FASTQ inputs
    ├── hg38_dir/                   # Genome references & STAR Index
    ├── qc_raw_reads/               # Initial FastQC reports
    ├── multiqc_raw/                # Aggregated Pre-trim QC
    ├── trimmed_fastq/              # Cleaned reads + fastp JSONs
    ├── qc_trimmed_reads/           # Post-trim FastQC
    ├── multiqc_trimmed/            # Aggregated Post-trim QC
    ├── star_mapped/                # Raw STAR alignment outputs
    ├── IGV/                        # Indexed BAMs (.bai) ready for visualization
    └── counts/                     # Final FeatureCounts matrix
