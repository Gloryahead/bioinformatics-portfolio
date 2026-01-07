# ===================================================================
# SCRIPT 6: GENERATE REPLICATE REPRODUCIBILITY HEATMAP (Corrected for pipe syntax)
# - Self-contained script for use in a SLURM job.
# - Automatically discovers and merges binned fragment count files.
# - Calculates a pairwise correlation matrix to assess reproducibility.
# - Saves the final correlation plot as a high-resolution PNG.
# ===================================================================


# --- 1. PACKAGE MANAGEMENT ---

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define list of required packages
cran_packages <- c("dplyr", "purrr", "corrplot", "RColorBrewer", "magrittr", "tidyr")

# Function to install missing packages
install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
}

# Install any missing packages
install_if_missing(cran_packages)

# Load all the required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(corrplot)
  library(RColorBrewer)
  library(magrittr)
  library(tidyr)
})


# --- 2. AUTOMATED FILE DISCOVERY AND EFFICIENT DATA LOADING ---

# Define the path to your binned bed files
projPath <- "/path/to/your/project"

# Automatically find all relevant fragment count files
count_files <- list.files(path = projPath, pattern = "_mm10\\.fragmentsCount\\.bin500\\.bed$")

# Error Handling: Stop if no files are found
if (length(count_files) == 0) {
  stop("Execution halted: No '...fragmentsCount.bin500.bed' files found in directory: '", projPath, "'")
}

print("Found files to process:")
print(count_files)

# Use a highly efficient pipeline to read, rename, and merge all files
fragCount <- map(count_files, function(file) {
  full_path <- file.path(projPath, file)
  sample_name <- gsub("_mm10\\.fragmentsCount\\.bin500\\.bed$", "", file)
  read.table(full_path, header = FALSE, col.names = c("chrom", "bin", sample_name))
}) %>%
  reduce(full_join, by = c("chrom", "bin")) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))


# --- 3. CALCULATE CORRELATION MATRIX ---

print("Calculating correlation matrix...")

# Calculate the correlation matrix on log2-transformed data.
M <- fragCount %>%
  select(-c(chrom, bin)) %>%
  mutate(across(everything(), ~ .x / sum(.x) * 1e6)) %>%  # normalize to counts per million
  {log2(. + 1)} %>%
  cor(use = "complete.obs")

# --- 4. GENERATE AND SAVE CORRELATION PLOT ---

# Define the output path for the final figure
output_figure_file <- file.path(projPath, "Figure_Replicate_Correlation.png")
corr_colors <- colorRampPalette(c("midnightblue", "white", "darkred"))(100)

# Dynamically determine the number of clusters (groups) to box.
num_clusters <- M %>% colnames() %>% sub("_[^_]*$", "", .) %>% unique() %>% length()
print(paste("Automatically detected", num_clusters, "histone groups to cluster."))

print(paste("Saving plot to:", output_figure_file))

# Open a PNG device to save the plot
png(output_figure_file, width = 10, height = 10, units = "in", res = 300)

# Create the correlation plot
corrplot(M, 
         method = "color", 
         col = corr_colors,
         outline = TRUE, 
         addgrid.col = "darkgray", 
         order = "hclust", 
         addrect = num_clusters,
         rect.col = "black", 
         rect.lwd = 3,
         cl.pos = "b", 
         tl.col = "indianred4", 
         tl.cex = 1.2, 
         cl.cex = 1, 
         addCoef.col = "black", 
         number.digits = 2, 
         number.cex = 1)

# Close the PNG device
dev.off()

print("--- Script 6 Complete ---")
