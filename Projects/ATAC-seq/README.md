# ATAC-seq Analysis Pipeline (HPC/SLURM)

![Bash](https://img.shields.io/badge/Language-Bash-4EAA25)
![Platform](https://img.shields.io/badge/Platform-SLURM%20%7C%20HPC-blue)
![Bioinformatics](https://img.shields.io/badge/Focus-Genomics-red)

## 📌 Overview
This repository contains a robust, end-to-end automated pipeline for **ATAC-seq (Assay for Transposase-Accessible Chromatin)** data analysis. 

Designed for High-Performance Computing (HPC) environments using **SLURM**, this pipeline processes raw sequencing data (FASTQ) through to peak calling and motif discovery. It features built-in idempotency (check-pointing), allowing the script to resume analysis for specific samples without re-running completed steps—critical for large genomic datasets.

## 🧬 Workflow Architecture
The pipeline executes the following biological and computational steps:

1.  **Quality Control:** Raw read analysis using `FastQC` and `MultiQC`.
2.  **Trimming:** Adapter and low-quality base removal using `fastp`.
3.  **Alignment:** Mapping to reference genome (mm10/hg38) using `Bowtie2` (optimized for sensitive local alignment).
4.  **Post-Processing:**
    * Sorting and Indexing (`Samtools`)
    * Duplicate Removal (`Picard MarkDuplicates`)
    * **Blacklist Filtering:** Removal of mitochondrial reads and ENCODE blacklisted regions.
5.  **Signal Generation:** Creation of normalized BigWig tracks using `deepTools` (RPKM normalized).
6.  **Peak Calling:**
    * Narrow peak calling using `MACS2`.
    * Peak annotation and motif discovery using `HOMER`.


## 🛠 Prerequisites & Dependencies
The pipeline assumes the following tools are available in the system path or Conda environment:
* `fastqc` / `multiqc`
* `fastp`
* `bowtie2`
* `samtools`
* `picard` (Java)
* `bedtools`
* `deeptools`
* `macs2`
* `homer`


## 🚀 Usage

### Configuration
Edit the `Configuration` section at the top of the script to set your file paths:

```bash


## Define your project base directory
base_dir="/path/to/your/project"
cores=12  # Adjust based on SLURM allocation

## 📂 Directory Output Structure
The pipeline organizes inputs and outputs into a strict hierarchy to ensure reproducibility:

```text
project_root/
├── qc_raw_reads/       # FastQC reports
├── trimmed_fastq/      # Cleaned reads
├── bam/                # Aligned & Filtered BAMs
├── bigwig/             # Visualization tracks (.bw)
├── MACS2/              # Peak calls
├── HOMER/              # Annotated peaks & Motifs
└── visualization/      # Ready-to-load BED files for IGV/UCSC
