# ===================================================================
# SCRIPT 1: PROCESS mm10 ALIGNMENTS
# - Self-contained script for use in a SLURM job.
# - Loads packages, processes mm10 files, and saves results.
# ===================================================================

# --- 1. PACKAGE MANAGEMENT ---
# This block is required because this script runs in its own session.

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define list of required packages
cran_packages <- c("stringr", "dplyr", "magrittr")
bioconductor_packages <- c("Biobase", "SummarizedExperiment", "MatrixGenerics", "matrixStats", "DESeq2", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges", "chromVAR")

# Function to install missing packages
install_if_missing <- function(pkgs, install_func) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install_func(pkg)
    }
  }
}

# Install any missing packages
install_if_missing(cran_packages, install.packages)
install_if_missing(bioconductor_packages, BiocManager::install)

# Load required libraries
suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(magrittr)
})


# --- 2. PROCESS PRIMARY GENOME (mm10) ALIGNMENTS ---

# Define paths
projPath_mm10 <- "/path/to/your/project/bowtie2_summary_files"
output_file <- file.path(projPath_mm10, "alignResult_mm10.csv") # Define output file

# Automatically find all .txt files
file_list_mm10 <- list.files(path = projPath_mm10, pattern = "\\.txt$")

# Error Handling
if (length(file_list_mm10) == 0) {
  stop("Execution halted: No '.txt' files found in directory: '", projPath_mm10, "'")
}

# Dynamically create lists
sampleList_mm10 <- tools::file_path_sans_ext(file_list_mm10)

# Initialize data frame
alignResult <- data.frame() 

# Loop and process files
for(hist in sampleList_mm10){
  file_path <- file.path(projPath_mm10, paste0(hist, ".txt"))
  alignRes <- read.table(file_path, header = FALSE, fill = TRUE)
  alignRate <- substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo <- strsplit(hist, "_")[[1]]
  
  tempResult <- data.frame(
    Targets = histInfo[1], 
    Alignment = histInfo[2], 
    SequencingDepth = as.numeric(as.character(alignRes$V1[1])), 
    MappedFragNum_mm10 = as.numeric(as.character(alignRes$V1[4])) + as.numeric(as.character(alignRes$V1[5])), 
    AlignmentRate_mm10 = as.numeric(alignRate)
  )
  
  alignResult <- rbind(alignResult, tempResult)
}

# --- 3. SAVE OUTPUT FOR NEXT SCRIPT ---

# Write the clean, numerical data to a CSV file.
# The next script in the pipeline will read this file.
write.csv(alignResult, file = output_file, row.names = FALSE)

print(paste("mm10 alignment results successfully saved to:", output_file))
print("--- Script 1 Complete ---")
