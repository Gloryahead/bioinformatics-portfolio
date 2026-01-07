# ==============================================================================
# RNA-SEQ ANALYSIS PIPELINE (EdgeR -> Limma-Voom)
# ==============================================================================
# This script performs the following:
# 1. Imports a raw FeatureCounts matrix.
# 2. Cleans sample names and removes non-count metadata.
# 3. Filters out lowly expressed genes.
# 4. Annotates genes (converts Ensembl IDs to Gene Symbols).
# 5. Normalizes data using TMM and transforms using Voom.
# 6. Generates Quality Control plots (PCA).
# 7. Performs Differential Expression Analysis.
# 8. Saves all results to the input directory.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. INSTALLATION & LIBRARY LOADING
# ------------------------------------------------------------------------------
# Check if BiocManager is installed; if not, install it.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required BioConductor packages (only if missing)
BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db", 
                       "org.Mm.eg.db", "org.Rn.eg.db", "AnnotationDbi"), update=FALSE)

# Install required CRAN packages (only if missing)
install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "stats", "tidyverse"))

# Load the libraries
library(data.table)     # Fast file reading/writing
library(edgeR)          # DGEList object and TMM normalization
library(limma)          # Linear modeling (voom, lmFit, eBayes)
library(ggplot2)        # Plotting
library(ggrepel)        # Non-overlapping text labels
library(ggfortify)      # PCA plotting support
library(org.Mm.eg.db)   # Mouse gene annotation database
library(AnnotationDbi)  # Functions for mapping IDs

# ------------------------------------------------------------------------------
# 2. DATA IMPORT & PATH SETUP
# ------------------------------------------------------------------------------

# DEFINITION: Path to your raw count matrix
# Update this path for new projects.
matrix_file <- "/path/to/your/matrix/file"

# SETUP: Define the Output Directory
# We extract the folder path from your input file. All results will be saved here.
output_dir <- dirname(matrix_file)
print(paste("Output Directory set to:", output_dir))

# READ: Load the data fast
# skip=1: Skips the command line info (lines starting with #) often found in FeatureCounts
# header=TRUE: Ensures we capture column names like Geneid, Chr, etc.
counts_data <- fread(matrix_file, skip = 1, header = TRUE)

# CHECK: View first few rows to confirm import worked
print("Raw data loaded:")
print(head(counts_data))


# ------------------------------------------------------------------------------
# 3. DATA CLEANING & FORMATTING
# ------------------------------------------------------------------------------

# SUBSET: Keep only GeneID and Sample Counts
# We keep column 1 (Geneid) and columns 7 through the end (Samples).
# Columns 2-6 (Chr, Start, End, Strand, Length) are removed.
counts_clean <- counts_data[, c(1, 7:ncol(counts_data)), with = FALSE]

# RENAME: Clean up sample names
# 1. Get current names (e.g., "/path/to/SRR28119110_Aligned.sorted.bam")
old_names <- colnames(counts_clean)
# 2. Remove directory paths (leaves "SRR28119110_Aligned.sorted.bam")
new_names <- basename(old_names)
# 3. Remove the suffix to get clean Sample IDs (leaves "SRR28119110")
new_names <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam", "", new_names)
# 4. Apply new names
colnames(counts_clean) <- new_names

# FORMAT: Convert to standard Data Frame
counts_final <- as.data.frame(counts_clean)
# Set the Row Names to be the GeneIDs
rownames(counts_final) <- counts_final$Geneid
# Remove the Geneid column now that it is the row name
counts_final$Geneid <- NULL 

# EXPORT: Save the clean raw counts
# We use file.path() to robustly combine the output directory and filename.
outfile_raw <- file.path(output_dir, "counts_raw_cleaned.tsv")
fwrite(counts_final, outfile_raw, sep = "\t", row.names = TRUE)
print(paste("Cleaned counts saved to:", outfile_raw))


# ------------------------------------------------------------------------------
# 4. FILTERING LOW COUNTS
# ------------------------------------------------------------------------------

# SETTING: Define filter strictness
# 0.8 means a gene must be non-zero in 80% of samples (e.g., 4 out of 4)
perc_keep <- 0.8 

# LOGIC: Identify genes to keep
# rowSums(counts_final > 0): How many samples have reads for this gene?
# ceiling(perc_keep * ncol): The threshold number of samples required.
gene_keep <- rowSums(counts_final > 0) >= ceiling(perc_keep * ncol(counts_final))

# APPLY: Subset the data
counts_final_low_rm <- counts_final[gene_keep, ]

print(paste("Genes before filtering:", nrow(counts_final)))
print(paste("Genes after filtering:", nrow(counts_final_low_rm)))


# ------------------------------------------------------------------------------
# 5. METADATA & ANNOTATION
# ------------------------------------------------------------------------------

# METADATA: Create table defining experimental groups
# WARNING: Ensure the order of 'Treatment' matches the columns in 'counts_final' EXACTLY.
# Check your columns: print(colnames(counts_final))
meta <- data.frame(SampleID = colnames(counts_final),
                   Treatment = c("KRAS_SPIB", "KRAS_SPIB", "KRAS", "KRAS"))
rownames(meta) <- meta$SampleID

# ANNOTATION: Convert Ensembl IDs to Gene Symbols
# 1. Strip version numbers from IDs (e.g., ENSMUSG...1 -> ENSMUSG...)
clean_ids <- gsub("\\..*", "", rownames(counts_final_low_rm))

# 2. Map IDs using the Mouse database
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = clean_ids,
                       column = "SYMBOL",      # Target: Gene Name
                       keytype = "ENSEMBL",    # Source: Ensembl ID
                       multiVals = "first")    # Handling duplicates

# 3. Handle NAs (Genes not found in database)
# If a name is NA, fill it with the original Ensembl ID so it isn't blank.
na_indices <- is.na(gene_symbols)
gene_symbols[na_indices] <- rownames(counts_final_low_rm)[na_indices]

# 4. Create Gene Info Table
# We separate annotation from counts to keep the count matrix numeric.
gene_info <- data.frame(
  GeneID = rownames(counts_final_low_rm),
  Symbol = gene_symbols
)


# ------------------------------------------------------------------------------
# 6. CREATE DGEList & NORMALIZE
# ------------------------------------------------------------------------------

# OBJECT: Create the edgeR container
# counts: The filtered numeric matrix (using 1:4 assumes 4 samples)
dge <- DGEList(counts = counts_final_low_rm[, 1:ncol(counts_final_low_rm)], 
               samples = meta, 
               genes = gene_info)

# ASSIGN GROUP: Map Treatment to the 'group' slot for edgeR
dge$samples$group <- factor(meta$Treatment, levels = c("KRAS", "KRAS_SPIB"))

# DESIGN MATRIX: Define the statistical model
# ~0 + group: Fit means for each group separately (no intercept)
design <- model.matrix(~0 + group, data = dge$samples)
colnames(design) <- levels(dge$samples$group) # Clean column names (KRAS, KRAS_SPIB)

# NORMALIZE: Calculate TMM scaling factors
dge <- calcNormFactors(dge, method = "TMM")

# TRANSFORM: Run Voom (Mean-Variance modeling)
# This prepares the data for linear modeling
dge_v <- voom(dge, design, plot = FALSE) # Set plot=TRUE to see the trend in RStudio

# SAVE: Save the processed object for later use
saveRDS(dge_v, file.path(output_dir, "dge_v.rds"))


# ------------------------------------------------------------------------------
# 7. QUALITY CONTROL (PCA PLOT)
# ------------------------------------------------------------------------------

# PREPARE: Extract data for plotting
meta_table <- dge_v$targets 
count_table_t <- as.data.frame(t(dge_v$E)) # Transpose for PCA (Samples = Rows)

# CALCULATE: Run PCA
pca_prep <- prcomp(count_table_t, scale. = TRUE)

# PLOT: Generate the PCA graph
# We disable automatic labels (label=FALSE) and add custom ones with geom_text_repel
pca_plot <- autoplot(pca_prep, 
                     data = meta_table, 
                     colour = "Treatment", 
                     shape = "Treatment", 
                     label = FALSE) + 
  geom_text_repel(aes(label = rownames(meta_table)), size = 4) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())

# SAVE: Save PCA to PDF
ggsave(file.path(output_dir, "PCA_Plot.pdf"), plot = pca_plot, 
       device = "pdf", units = "cm", width = 16, height = 14)
print("PCA Plot saved.")


# ------------------------------------------------------------------------------
# 8. DIFFERENTIAL EXPRESSION ANALYSIS
# ------------------------------------------------------------------------------

# DEFINE CONTRAST: What comparison are we making?
# We compare KRAS_SPIB minus KRAS (Positive logFC = Upregulated in KRAS_SPIB)
comparison <- "KRAS_SPIB-KRAS"
contrast_matrix <- makeContrasts(contrasts = comparison, levels = design)

# FIT MODEL: Run the linear model steps
fit <- lmFit(dge_v, design)               # Fit model to every gene
fit_2 <- contrasts.fit(fit, contrast_matrix) # Apply the specific comparison
fit_2 <- eBayes(fit_2)                    # Apply Empirical Bayes smoothing

# EXTRACT RESULTS: Get the table of all genes
# n = Inf: Keep all genes
# sort.by = "p": Sort by significance
deg_tbl <- topTable(fit_2, coef = colnames(contrast_matrix), n = Inf, sort.by = "p")

# SAVE: Write the full results to file
outfile_deg <- file.path(output_dir, "DEG_KRAS_SPIB_vs_KRAS.tsv")
fwrite(deg_tbl, outfile_deg, sep = "\t", row.names = TRUE)

print(paste("Analysis Complete. Results saved to:", outfile_deg))
