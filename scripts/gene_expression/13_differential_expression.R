# ============================================================================
# 13_differential_expression.R
# ============================================================================
#
# What this script does
# ---------------------
# Runs differential expression analysis comparing wildtype (WT) against each
# translocated cell line (T1 and C1) using DESeq2.
#
# For each comparison (WT vs T1, WT vs C1), the script:
#   1. Loads raw gene counts and sample metadata.
#   2. Builds a DESeq2 dataset and filters very low-count genes.
#   3. Estimates size factors and dispersions (DESeq2 standard pipeline).
#   4. Shrinks log2 fold changes using the apeglm method, which reduces
#      noise for genes with low counts or high variability.
#   5. Saves results as TSV files for use in downstream analysis (script 14).
#   6. Runs basic sanity checks: result summaries and MA plots.
#
# Why apeglm shrinkage?
# ---------------------
# Raw log2 fold changes from DESeq2 are noisy for genes with low counts —
# a gene with only 3 counts in one condition and 1 in another will have a
# large raw fold change but very little statistical support. apeglm shrinks
# these unreliable estimates toward zero, producing more conservative and
# interpretable fold changes that are better suited for downstream ranking
# and visualisation.
#
# Output
# ------
#   DE_WT_vs_T1.tsv  — DESeq2 results for WT vs T1 (one row per gene)
#   DE_WT_vs_C1.tsv  — DESeq2 results for WT vs C1
#
# Columns in output files:
#   gene_id        : Ensembl gene ID
#   gene_name      : HGNC gene symbol
#   baseMean       : mean normalised count across all samples
#   log2FoldChange : apeglm-shrunken log2 fold change (condition vs WT)
#   lfcSE          : standard error of the shrunken log2FC
#   svalue         : s-value (analogous to p-value for shrunken estimates)
#   padj           : BH-adjusted p-value
#
# Usage
# -----
#   1. Edit the CONFIG section below to point to your files.
#   2. Run: Rscript 13_differential_expression.R
#      or source the script from RStudio.
#
# Dependencies
# ------------
#   DESeq2  (Bioconductor)
#   Install with:
#     if (!requireNamespace("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#     BiocManager::install("DESeq2")
#
# ============================================================================


# ============================================================================
# CONFIG – edit all paths here before running
# ============================================================================

# Raw gene count matrix (genes x samples, tab-separated).
# Must include a 'gene_id' column and a 'gene_name' column.
COUNTS_FILE <- "/path/to/GSE246689_gene_counts.tsv"

# Sample metadata table (tab-separated).
# Must have sample IDs as row names and at least a 'condition' column
# with values matching the column names of the counts file.
# Expected conditions: WT, T1, C1.
COLDATA_FILE <- "/path/to/coldata.tsv"

# Where to save the output TSV files (created automatically if needed)
OUTPUT_DIR <- "/path/to/results/DE"

# Minimum total count across all samples for a gene to be included.
# Genes with fewer counts than this are filtered before fitting the model.
MIN_COUNT <- 10

# ============================================================================


# ============================================================================
# Load libraries
# ============================================================================

library(DESeq2)


# ============================================================================
# Load and prepare input data
# ============================================================================

cat("Loading count matrix...\n")
counts_raw <- read.table(COUNTS_FILE, header = TRUE, sep = "\t")

# Store gene ID → gene name mapping for adding back after analysis
gene_map <- counts_raw[, c("gene_id", "gene_name")]

# Use gene_id as row names and remove annotation columns
rownames(counts_raw) <- counts_raw$gene_id
counts <- counts_raw[, -(1:2)]   # removes gene_id and gene_name columns

# DESeq2 requires integer counts — round to be safe
counts <- round(counts)
stopifnot(all(counts >= 0))
stopifnot(all(counts == floor(counts)))

cat("Loading sample metadata...\n")
coldata <- read.table(COLDATA_FILE, header = TRUE, row.names = 1, sep = "\t")

# Set WT as the reference level so all comparisons are condition vs WT
coldata$condition <- factor(coldata$condition, levels = c("WT", "T1", "C1"))

# Ensure sample order in the count matrix matches the metadata
counts <- counts[, rownames(coldata)]

cat(sprintf("Count matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("Conditions: %s\n", paste(levels(coldata$condition), collapse = ", ")))


# ============================================================================
# Build DESeq2 dataset and run analysis
# ============================================================================

cat("\nBuilding DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ condition
)

# Filter genes with very low total counts across all samples.
# These genes have insufficient data to estimate reliable fold changes
# and their inclusion would inflate multiple testing correction burden.
dds <- dds[rowSums(counts(dds)) >= MIN_COUNT, ]
cat(sprintf("Genes after low-count filtering (>= %d counts): %d\n",
            MIN_COUNT, nrow(dds)))

cat("\nRunning DESeq2...\n")
dds <- DESeq(dds)

# Print the coefficient names so the correct coef can be verified
cat("\nCoefficient names:\n")
print(resultsNames(dds))


# ============================================================================
# Extract and shrink results
# ============================================================================

cat("\nExtracting results with apeglm shrinkage...\n")

# apeglm shrinkage reduces noise in log2FC estimates for low-count genes,
# producing more conservative and interpretable fold changes
res_T1 <- lfcShrink(dds, coef = "condition_T1_vs_WT", type = "apeglm")
res_C1 <- lfcShrink(dds, coef = "condition_C1_vs_WT", type = "apeglm")


# ============================================================================
# Save results
# ============================================================================

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

save_results <- function(res, cond, gene_map, output_dir) {
  df          <- as.data.frame(res)
  df$gene_id  <- rownames(df)
  df          <- merge(df, gene_map, by = "gene_id")
  out_path    <- file.path(output_dir, sprintf("DE_WT_vs_%s.tsv", cond))
  write.table(df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Saved: %s  (%d genes)\n", out_path, nrow(df)))
  invisible(df)
}

save_results(res_T1, "T1", gene_map, OUTPUT_DIR)
save_results(res_C1, "C1", gene_map, OUTPUT_DIR)


# ============================================================================
# Sanity checks
# ============================================================================

cat("\n--- WT vs T1 summary ---\n")
print(summary(res_T1))
cat(sprintf("Significant genes (padj < 0.05): %d\n",
            sum(res_T1$padj < 0.05, na.rm = TRUE)))

cat("\n--- WT vs C1 summary ---\n")
print(summary(res_C1))
cat(sprintf("Significant genes (padj < 0.05): %d\n",
            sum(res_C1$padj < 0.05, na.rm = TRUE)))

# MA plots: shows log2FC vs mean expression.
# Points in blue are significantly differentially expressed.
# A well-behaved MA plot should be centred around 0 with symmetric scatter.
plotMA(res_T1, main = "MA plot: WT vs T1")
plotMA(res_C1, main = "MA plot: WT vs C1")

cat("\nDone.\n")
