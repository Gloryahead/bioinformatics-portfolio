# ==============================================================================
# RNA-SEQ VISUALIZATION PIPELINE: Volcano Plots & Heatmaps
# ==============================================================================
# Purpose: Generate publication-ready visualizations from Differential Expression data.
# Inputs:  1. Differential Expression Table (TSV) from Limma/EdgeR
#          2. DGEList object (RDS) containing normalized counts
# Outputs: 1. EnhancedVolcano Plot (PNG)
#          2. Custom Labeled Volcano Plot (PDF)
#          3. Top DEGs Heatmap (PDF)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. CONFIGURATION & PATHS (Change this section for new projects)
# ------------------------------------------------------------------------------

# Set this to TRUE if you need to install packages. Set to FALSE for daily use.
RUN_INSTALL <- FALSE 

# Define the full path to your Differential Expression results file
input_deg_file <- "/path/to/your/DEG_P1_lfc0.tsv"

# Define the path to your R object file (dge_v)
input_rds_file <- "/path/to/your/dge_v.rds"

# Automatically extract the directory from the input file
# All output files will be saved here.
output_dir <- dirname(input_deg_file)

# Significance Thresholds
FC_threshold <- 2      # Fold Change (2 = log2(1))
P_threshold  <- 0.05   # Adjusted P-value (FDR)


# ------------------------------------------------------------------------------
# 2. PACKAGE INSTALLATION & LOADING
# ------------------------------------------------------------------------------
if (RUN_INSTALL) {
  # Install CRAN packages
  install.packages(c("data.table", "dplyr", "tidyr", "readr", "stringr", 
                     "ggplot2", "ggrepel", "ggfortify", "ggprism", "pheatmap"))
  
  # Install Bioconductor packages
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("EnhancedVolcano")
}

# Load Libraries
# We load individual tidyverse components to avoid system dependency issues (ragg)
library(data.table)      # Fast file reading
library(dplyr)           # Data manipulation
library(tidyr)           # Data reshaping
library(stringr)         # String manipulation
library(ggplot2)         # Plotting core
library(ggrepel)         # Smart labels that don't overlap
library(pheatmap)        # Heatmaps
library(EnhancedVolcano) # Pre-packaged Volcano plots
library(grid)            # Required for saving heatmaps

print(paste("Working Directory set to:", output_dir))


# ------------------------------------------------------------------------------
# 3. DATA IMPORT
# ------------------------------------------------------------------------------

# Load the Differential Expression Table
deg_tbl <- fread(input_deg_file)

# Load the DGEList object (contains normalized counts for heatmap)
dge_v <- readRDS(input_rds_file)

# Check if 'Symbol' column exists, if not, try to use V1 or GeneID
if (!"Symbol" %in% colnames(deg_tbl)) {
  warning("Column 'Symbol' not found. Using 'V1' or row names for labels.")
  deg_tbl$Symbol <- deg_tbl$V1 # Fallback
}


# ------------------------------------------------------------------------------
# 4. PLOT 1: ENHANCED VOLCANO (The Quick Look)
# ------------------------------------------------------------------------------

# Create custom color logic for EnhancedVolcano
# We define colors based on thresholds defined in Section 1
keyvals.colour <- ifelse(
  deg_tbl$logFC < -log2(FC_threshold) & deg_tbl$adj.P.Val < P_threshold, 'blue',
  ifelse(deg_tbl$logFC > log2(FC_threshold) & deg_tbl$adj.P.Val < P_threshold, 'red', 'grey30')
)

# Assign names to the colors for the legend
# We calculate the number of genes (n=) dynamically
names(keyvals.colour)[keyvals.colour == 'blue'] <- paste0("Downregulated (n=", sum(keyvals.colour == 'blue'), ")")
names(keyvals.colour)[keyvals.colour == 'red']  <- paste0("Upregulated (n=", sum(keyvals.colour == 'red'), ")")
names(keyvals.colour)[keyvals.colour == 'grey30'] <- 'Non-significant'

# Generate the plot
volcano_plot <- EnhancedVolcano(
  deg_tbl,
  lab = NA,                # Turn off labels here (we use custom plot for labels)
  x = 'logFC',
  y = 'adj.P.Val',
  title = 'Global Differential Expression',
  subtitle = paste0("FDR < ", P_threshold, ", FC > ", FC_threshold),
  pCutoff = P_threshold,
  FCcutoff = log2(FC_threshold),
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  ylim = c(0, max(-log10(deg_tbl$adj.P.Val), na.rm=TRUE) + 0.5), # Auto-scale Y-axis
  colCustom = keyvals.colour,
  colAlpha = 0.8,
  pointSize = 1.0,
  legendPosition = 'right'
)

# Save Plot 1
outfile_v1 <- file.path(output_dir, "Volcano_Enhanced.png")
ggsave(outfile_v1, volcano_plot, device = "png", width = 16, height = 14, units = "cm")
print(paste("Saved:", outfile_v1))


# ------------------------------------------------------------------------------
# 5. PLOT 2: CUSTOM GGPLOT VOLCANO (With Labels)
# ------------------------------------------------------------------------------

# Add columns for plotting logic
deg_tbl$neg_log10_adj.P.Val <- -log10(deg_tbl$adj.P.Val)

# Classify genes into categories
deg_tbl$significance <- ifelse(
  deg_tbl$adj.P.Val < P_threshold & abs(deg_tbl$logFC) > log2(FC_threshold),
  ifelse(deg_tbl$logFC > 0, "Upregulated", "Downregulated"),
  "Non-Significant"
)

# Pick the top 15 most significant genes to label
# We order by P-value and select the 'Symbol' column
genes_to_label <- deg_tbl %>% arrange(adj.P.Val) %>% pull(Symbol) %>% head(15)

# Build the custom plot
custom_volcano <- ggplot(deg_tbl, aes(x = logFC, y = neg_log10_adj.P.Val)) +
  # Main scatter points
  geom_point(aes(color = significance), alpha = 0.8, size = 1) +
  
  # Manual colors matching our logic
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-Significant" = "gray80")) +
  
  # Add Threshold lines
  geom_vline(xintercept = c(-log2(FC_threshold), log2(FC_threshold)), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(P_threshold), linetype = "dashed", color = "gray50") +
  
  # Smart Labeling (Only labels genes in our top 15 list)
  geom_text_repel(
    data = subset(deg_tbl, Symbol %in% genes_to_label), 
    aes(label = Symbol),
    size = 3,
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  
  # Labels and Theme
  labs(x = expression(log[2]~"Fold Change"), 
       y = expression(-log[10]~"adjusted p-value"), 
       color = "Status") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

# Save Plot 2
outfile_v2 <- file.path(output_dir, "Volcano_Labeled_Symbols.pdf")
ggsave(outfile_v2, custom_volcano, device = "pdf", width = 10, height = 8)
print(paste("Saved:", outfile_v2))


# ------------------------------------------------------------------------------
# 6. PLOT 3: HEATMAP OF SIGNIFICANT GENES
# ------------------------------------------------------------------------------

# A. Prepare Annotation (the colored bar at the top of the heatmap)
# We map samples to their treatments using the 'targets' from the dge_v object
col_annot_df <- data.frame(
  Treatment = dge_v$targets$Treatment,
  row.names = rownames(dge_v$targets)
)
# Order the annotation so samples from the same group appear together
col_annot_df <- col_annot_df[order(col_annot_df$Treatment), , drop = FALSE]

# B. Select Data
# 1. Identify IDs of significant genes (Using V1/GeneID to match expression matrix)
sig_gene_ids <- deg_tbl$V1[abs(deg_tbl$logFC) > log2(FC_threshold) & deg_tbl$adj.P.Val < P_threshold]

# 2. Extract expression values for these genes
# We also re-order the columns to match our sorted annotation frame
expr_deg <- dge_v$E[sig_gene_ids, rownames(col_annot_df), drop = FALSE]

# C. Generate Heatmap
heatmap_custom <- pheatmap(
  expr_deg,
  scale = "row",            # Z-score scaling (highlights relative differences)
  cluster_rows = TRUE,      # Group similar genes
  cluster_cols = FALSE,     # Don't cluster samples (we want them grouped by Treatment)
  show_rownames = FALSE,    # Hide row names if there are too many genes
  show_colnames = TRUE,
  color = colorRampPalette(c("navy", "white", "red"))(100),
  annotation_col = col_annot_df,
  main = paste("Significant DEGs (n =", length(sig_gene_ids), ")")
)

# D. Save Heatmap
outfile_hm <- file.path(output_dir, "Heatmap_Significant_Genes.pdf")

# Custom save function needed for pheatmap objects
save_pheatmap_pdf <- function(x, filename, width=8, height=8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(heatmap_custom, outfile_hm, width = 6, height = 8)
print(paste("Saved:", outfile_hm))

print("Analysis Complete.")
