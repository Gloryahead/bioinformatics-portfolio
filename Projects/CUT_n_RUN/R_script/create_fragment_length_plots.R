# ===================================================================
# SCRIPT 4: GENERATE FRAGMENT LENGTH PLOTS (Corrected for ggplot2 v3.4.0+)
# - Self-contained script for use in a SLURM job.
# - Automatically discovers fragment length files, processes them,
#   and generates a two-panel summary figure.
# ===================================================================


# --- 1. PACKAGE MANAGEMENT ---
# Required for every self-contained script.

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define list of required packages
cran_packages <- c("ggplot2", "dplyr", "purrr", "viridis", "ggpubr", "magrittr")

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
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(viridis)
  library(ggpubr)
  library(magrittr)
})


# --- 2. AUTOMATED FILE DISCOVERY AND DATA LOADING ---

# Define the path to your fragment length files
projPath <- "/path/to/your/project"

# Automatically find all files that end with "_fragmentLen.txt"
frag_files <- list.files(path = projPath, pattern = "_fragmentLen\\.txt$")

# Error Handling: Stop if no files are found
if (length(frag_files) == 0) {
  stop("Execution halted: No '_fragmentLen.txt' files found in directory: '", projPath, "'")
}

# Use purrr::map_dfr for a fast and efficient way to read and combine all files
fragLen <- map_dfr(frag_files, function(file) {
  full_path <- file.path(projPath, file)
  sampleInfo <- gsub("_fragmentLen\\.txt$", "", file)
  histInfo <- sub("_[^_]*$", "", sampleInfo)
  replicateInfo <- sub(".*_", "", sampleInfo)
  
  read.table(full_path, col.names = c("fragLen", "fragCount")) %>%
    mutate(
      Weight = fragCount / sum(fragCount),
      Histone = histInfo,
      Replicate = replicateInfo,
      sampleInfo = sampleInfo
    )
})


# --- 3. PREPARE DATA FOR PLOTTING (SET FACTOR LEVELS) ---

ordered_samples <- gsub("_fragmentLen\\.txt$", "", frag_files)
ordered_histones <- unique(sub("_[^_]*$", "", ordered_samples))

fragLen$sampleInfo <- factor(fragLen$sampleInfo, levels = ordered_samples)
fragLen$Histone <- factor(fragLen$Histone, levels = ordered_histones)


# --- 4. GENERATE INDIVIDUAL PLOTS ---

# A. Violin plot of fragment size distribution
fig2A <- fragLen %>% 
  ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 100), limits = c(0, 800)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Fragment Length (bp)") +
  xlab(NULL) +
  ggtitle("Fragment Length Distribution")

# B. Line plot of raw fragment counts
fig2B <- fragLen %>% 
  ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Replicate)) +
  # **FIX APPLIED HERE:** Changed 'size' to 'linewidth' to align with modern ggplot2
  geom_line(linewidth = 1) + 
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length (bp)") +
  ylab("Raw Read Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("Fragment Counts vs. Length")


# --- 5. COMBINE PLOTS AND SAVE TO FILE ---

# Arrange the two plots side-by-side
final_figure <- ggarrange(fig2A, fig2B, ncol = 2, widths = c(1.2, 1))

# Define the output path for the final figure
output_figure_file <- file.path(projPath, "Figure_Fragment_Length_Analysis.png")

# Save the combined figure to a file
ggsave(output_figure_file, final_figure, width = 18, height = 9, dpi = 300)

print(paste("Final figure successfully saved to:", output_figure_file))
print("--- Script 4 Complete ---")
