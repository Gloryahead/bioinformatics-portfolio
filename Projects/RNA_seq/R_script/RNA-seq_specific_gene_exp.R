# ==============================================================================
# TARGETED GENE EXPRESSION PLOTTING (Box, Violin, Bar)
# ==============================================================================
# This script performs the following:
# 1. Defines input files and genes of interest.
# 2. Checks and installs necessary visualization libraries.
# 3. Loads the processed DGEList object (.rds file).
# 4. Subsets data for specific genes and reshapes it (Wide -> Long).
# 5. Generates three types of plots:
#    - Box Plot (Distribution quartiles)
#    - Violin Plot (Density distribution)
#    - Bar Plot (Mean +/- Standard Error)
# 6. Saves all plots to the same directory as the input file.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. INPUT SETTINGS & LIBRARY LOADING
# ------------------------------------------------------------------------------

# DEFINITION: Input Parameters
# Update these paths and genes for different analyses.
input_file <- "/path/to/your/dge_v.rds"
my_genes   <- c("Acan", "Sox9", "Col2a1", "Agc1") 

# CHECK: Install CRAN packages if missing
# We check for ggplot2, tidyr, and dplyr.
required_packages <- c("ggplot2", "tidyr", "dplyr")
missing_packages  <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages)) {
  message("Installing missing packages...")
  install.packages(missing_packages)
}

# LOAD: Load the libraries
library(ggplot2) # Plotting engine
library(tidyr)   # Data reshaping (pivot_longer)
library(dplyr)   # Data manipulation (group_by, summarise)
library(tools)   # File path handling


# ------------------------------------------------------------------------------
# 2. DATA IMPORT & PATH SETUP
# ------------------------------------------------------------------------------

# SETUP: Define Output Directory
# We extract the folder path from your input file. Results go here.
output_dir <- dirname(input_file)
base_name  <- file_path_sans_ext(basename(input_file))
print(paste("Output Directory set to:", output_dir))

# READ: Load the processed RNA-seq object
dge_v <- readRDS(input_file)
print("Data loaded successfully.")


# ------------------------------------------------------------------------------
# 3. DATA PROCESSING & RESHAPING
# ------------------------------------------------------------------------------

# *** CRITICAL FIX: TRANSFORM DATA TO MATCH "FIRST PLOT" SCALE ***
# This shifts all values to ensure they are positive, matching your manual method.
dge_v$E <- dge_v$E + abs(min(dge_v$E)) + 1

# MATCH: Find Ensembl IDs for Gene Symbols
found_genes <- dge_v$genes[dge_v$genes$Symbol %in% my_genes, ]

# CHECK: verify we found the genes
if (nrow(found_genes) == 0) stop("No genes matched your list!")
print(paste("Found", nrow(found_genes), "genes."))

# SUBSET: Extract Expression Matrix
# We pull the numeric data (E) for only the rows matching our genes.
subset_expression <- t(dge_v$E[rownames(found_genes), , drop=FALSE])

# FORMAT: Prepare the dataframe
plot_data <- as.data.frame(subset_expression)
colnames(plot_data) <- found_genes$Symbol 

# ANNOTATE: Add Metadata
plot_data$Treatment <- dge_v$targets$Treatment

# TRANSFORM: Reshape Wide -> Long
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = all_of(found_genes$Symbol), 
    names_to = "Gene",                  
    values_to = "Expression"            
  )

# ------------------------------------------------------------------------------
# 4. PLOTTING
# ------------------------------------------------------------------------------

# PLOT A: Box Plot
# Shows median, quartiles, and outliers.
p_box <- ggplot(plot_data_long, aes(x = Treatment, y = Expression)) +
  geom_boxplot(aes(fill = Treatment), alpha = 0.7) +   
  geom_jitter(width = 0.2, alpha = 0.6) +              
  theme_minimal() +
  labs(title = "Box Plot: Expression Levels", y = "log2 Expression") +
  theme(legend.position = "none") +
  facet_wrap(~ Gene, scales = "free_y") 

# PLOT B: Violin Plot
# Shows the probability density of the data.
p_violin <- ggplot(plot_data_long, aes(x = Treatment, y = Expression, fill = Treatment)) +
  geom_violin(alpha = 0.6, scale = "width") +          
  geom_jitter(width = 0.2, alpha = 0.8, size = 2) +    
  stat_summary(fun = mean, geom = "point", size = 3) + # Add Mean dot
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + 
  theme_minimal() +
  labs(title = "Violin Plot: Expression Distribution", y = "log2 Expression") +
  theme(legend.position = "none") +
  facet_wrap(~ Gene, scales = "free_y")

# PLOT C: Bar Plot
# Shows Mean + Standard Error (Requires summarizing data first).
summary_data <- plot_data_long %>%
  group_by(Treatment, Gene) %>%
  summarise(
    mean = mean(Expression),
    se = sd(Expression)/sqrt(n()),
    .groups = "drop"
  )

p_bar <- ggplot(summary_data, aes(x = Treatment, y = mean, fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.7) +           
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) + 
  theme_minimal() +
  labs(title = "Bar Plot: Mean Expression", y = "Mean log2 Expression") +
  theme(legend.position = "none") +
  facet_wrap(~ Gene, scales = "free_y")


# ------------------------------------------------------------------------------
# 5. SAVING OUTPUTS
# ------------------------------------------------------------------------------

# SAVE: Export plots to disk
ggsave(file.path(output_dir, paste0(base_name, "_BoxPlot.png")), p_box, width = 8, height = 6)
ggsave(file.path(output_dir, paste0(base_name, "_ViolinPlot.png")), p_violin, width = 8, height = 6)
ggsave(file.path(output_dir, paste0(base_name, "_BarPlot.png")), p_bar, width = 8, height = 6)

print(paste("All plots saved to:", output_dir))
