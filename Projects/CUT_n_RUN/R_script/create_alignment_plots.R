# ===================================================================
# SCRIPT 3: GENERATE ALIGNMENT QUALITY PLOTS
# - Self-contained script for use in a SLURM job.
# - Loads packages, reads processed data, and generates a multi-panel figure.
# ===================================================================


# --- 1. PACKAGE MANAGEMENT ---
# This block is required as the script runs in its own, new session.

# Ensure BiocManager is installed (though not strictly needed for these packages, it's good practice)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define list of required plotting packages
cran_packages <- c("ggplot2", "dplyr", "viridis", "ggpubr", "magrittr")

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
  library(viridis)
  library(ggpubr)
  library(magrittr)
})


# --- 2. LOAD PROCESSED DATA ---

# Define paths to the data files created by previous scripts
projPath_mm10 <- "/path/to/your/project"
projPath_spike <- "/path/to/your/project"

# It's good practice to create specific filenames for each output
mm10_results_file <- file.path(projPath_mm10, "alignResult_mm10.csv")
spike_results_file <- file.path(projPath_spike, "alignResult_spikeIn.csv") # Assumes Script 2 saves this

# Error handling: Check if the input files exist before proceeding
if (!file.exists(mm10_results_file)) {
  stop("Execution halted: Input file not found: '", mm10_results_file, "'. Please run Script 1 first.")
}
if (!file.exists(spike_results_file)) {
  stop("Execution halted: Input file not found: '", spike_results_file, "'. Please run the spike-in processing script first.")
}

# Load the data into data frames
alignResult <- read.csv(mm10_results_file)
spikeAlign <- read.csv(spike_results_file)


# --- 3. PREPARE DATA FOR PLOTTING (RE-CREATE FACTOR LEVELS) ---

# The order of samples on the x-axis is determined by factor levels. This is lost in CSV files.
# We re-create the intended order by dynamically reading the source filenames again.
# This ensures the plot order matches the file processing order.
file_list_mm10 <- list.files(path = projPath_mm10, pattern = "\\.txt$")
sampleList_mm10 <- tools::file_path_sans_ext(file_list_mm10)
histone_order <- unique(sapply(strsplit(sampleList_mm10, "_"), `[`, 1))

# Apply the correct factor ordering to both data frames for consistent plots
alignResult$Targets <- factor(alignResult$Targets, levels = histone_order)
spikeAlign$Targets <- factor(spikeAlign$Targets, levels = histone_order)


# --- 4. GENERATE INDIVIDUAL PLOTS ---

# A. Sequencing Depth
fig1A <- alignResult %>% 
  ggplot(aes(x = Targets, y = SequencingDepth/1000000, fill = Targets)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(aes(color = Alignment), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Sequencing Depth (Millions)") +
  xlab(NULL) + 
  ggtitle("A. Sequencing Depth")

# B. Alignable Fragments (mm10)
fig1B <- alignResult %>% 
  ggplot(aes(x = Targets, y = MappedFragNum_mm10/1000000, fill = Targets)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(aes(color = Alignment), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Mapped Fragments (Millions)") +
  xlab(NULL) +
  ggtitle("B. Alignable Fragments (mm10)")

# C. Alignment Rate (mm10)
fig1C <- alignResult %>% 
  ggplot(aes(x = Targets, y = AlignmentRate_mm10, fill = Targets)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(aes(color = Alignment), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("% Mapped Fragments") +
  xlab(NULL) +
  ggtitle("C. Alignment Rate (mm10)")

# D. Alignment Rate (E. coli Spike-in)
fig1D <- spikeAlign %>% 
  ggplot(aes(x = Targets, y = AlignmentRate_spikeIn, fill = Targets)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(aes(color = Alignment), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("% Spike-in Alignment") +
  xlab(NULL) +
  ggtitle("D. Alignment Rate (E.coli)")


# --- 5. COMBINE PLOTS AND SAVE TO FILE ---

# Arrange the four plots into a 2x2 grid.
# The `common.legend` argument is powerful and will create a single legend for all plots.
final_figure <- ggarrange(fig1A, fig1B, fig1C, fig1D, 
                          ncol = 2, nrow = 2, 
                          common.legend = TRUE, legend = "bottom")

# Define the output path for the final figure
output_figure_file <- file.path(projPath_mm10, "Figure_Alignment_Summary.png")

# Save the combined figure to a file
ggsave(output_figure_file, final_figure, width = 12, height = 10, dpi = 300)

print(paste("Final figure successfully saved to:", output_figure_file))
print("--- Script 3 Complete ---")
