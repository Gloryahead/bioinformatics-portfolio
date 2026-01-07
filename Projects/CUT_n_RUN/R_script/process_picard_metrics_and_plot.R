# ===================================================================
# SCRIPT 5: PROCESS PICARD METRICS AND GENERATE PLOTS (Corrected for Data Type Error)
# - Self-contained script for use in a SLURM job.
# - Loads the main summary file from previous steps.
# - Automatically discovers and parses Picard duplication files.
# - Merges the data and saves a final multi-panel figure.
# ===================================================================


# --- 1. PACKAGE MANAGEMENT ---

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define list of required packages for this specific task
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


# --- 2. LOAD DEPENDENT DATA FROM PREVIOUS STEPS ---

# Define paths
summary_path <- "/path/to/your/project/csv_files"
merged_summary_file <- file.path(summary_path, "alignSummary_merged.csv")

# Error handling: Check if the main summary file exists
if (!file.exists(merged_summary_file)) {
  stop("Execution halted: Merged summary file not found: '", merged_summary_file, "'. Please run previous scripts first.")
}

# Load the main alignment summary data
alignSummary <- read.csv(merged_summary_file)

# 🔹 Extract genome name dynamically (e.g., "mm10")
genome <- sub(".*sam_(.*?)_summary.*", "\\1", summary_path)

# 🔹 Add Alignment column for joining
alignSummary$Alignment <- genome

# --- 3. AUTOMATED DISCOVERY AND LOADING OF PICARD DATA ---

# Define the path to your picard summary files
projPath <- "/home/michael/Downloads/picard_summary/CR12"

# Automatically find all Picard output files
picard_files <- list.files(path = projPath, pattern = "_picard\\.rmDup\\.txt$")

# Error Handling: Stop if no files are found
if (length(picard_files) == 0) {
  stop("Execution halted: No '_picard.rmDup.txt' files found in directory: '", projPath, "'")
}

# Use purrr::map_dfr to efficiently read and combine all Picard files
dupResult <- map_dfr(picard_files, function(file) {
  full_path <- file.path(projPath, file)
  
  sampleInfo <- gsub("_picard\\.rmDup\\.txt$", "", file)
  histInfo <- strsplit(sampleInfo, "_")[[1]][1]
  replicateInfo <- strsplit(sampleInfo, "_")[[1]][2]
  
  # Read the Picard metrics file
  # **FIX 1:** Added 'comment.char = "#"' to properly ignore Picard's header lines.
  read.table(full_path, header = TRUE, fill = TRUE, sep = "\t", comment.char = "#") %>%
    head(1) %>%
    # **FIX 2:** Use as.numeric() to guarantee the columns are numeric before math operations.
    transmute(
      Targets = histInfo,
      Alignment = replicateInfo,
      MappedFragNum_hg38 = as.numeric(READ_PAIRS_EXAMINED),
      DuplicationRate = as.numeric(PERCENT_DUPLICATION) * 100,
      EstimatedLibrarySize = as.numeric(ESTIMATED_LIBRARY_SIZE),
      UniqueFragNum = as.numeric(READ_PAIRS_EXAMINED) * (1 - as.numeric(PERCENT_DUPLICATION))
    )
})

# 🔹 Ensure same Alignment value ("mm10") for joining
dupResult$Alignment <- genome

# --- 4. PREPARE DATA AND MERGE ---

ordered_histones <- unique(map_chr(strsplit(gsub("_picard\\.rmDup\\.txt$", "", picard_files), "_"), 1))
dupResult$Targets <- factor(dupResult$Targets, levels = ordered_histones)

alignDupSummary <- left_join(alignSummary, dupResult, by = c("Targets", "Alignment"))

alignDupSummary_display <- alignDupSummary %>%
  mutate(DuplicationRate_str = paste0(round(DuplicationRate, 2), "%"))


# --- 5. GENERATE PLOTS ---

# A. Duplication Rate
fig3A <- dupResult %>% 
  ggplot(aes(x = Targets, y = DuplicationRate, fill = Targets)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(aes(color = Alignment), position = position_jitter(0.15), size = 2) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Duplication Rate (%)") +
  xlab(NULL) +
  ggtitle("A. Duplication Rate")

# B. Estimated Library Size
fig3B <- dupResult %>% 
  ggplot(aes(x = Targets, y = EstimatedLibrarySize, fill = Targets)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(aes(color = Alignment), position = position_jitter(0.15), size = 2) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Estimated Library Size") +
  xlab(NULL) +
  ggtitle("B. Estimated Library Size")

# C. Number of Unique Fragments
fig3C <- dupResult %>% 
  ggplot(aes(x = Targets, y = UniqueFragNum, fill = Targets)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(aes(color = Alignment), position = position_jitter(0.15), size = 2) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("# of Unique Fragments") +
  xlab(NULL) +
  ggtitle("C. Unique Fragments")


# --- 6. COMBINE PLOTS AND SAVE TO FILE ---

final_figure <- ggarrange(fig3A, fig3B, fig3C, 
                          ncol = 3, 
                          common.legend = TRUE, legend = "bottom")

output_figure_file <- file.path(projPath, "Figure_Picard_Duplication_Metrics.png")

ggsave(output_figure_file, final_figure, width = 20, height = 8, dpi = 300)

print(paste("Final figure successfully saved to:", output_figure_file))
print("--- Script 5 Complete ---")
