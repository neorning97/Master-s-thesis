# ============================================================================
# 13_differential_expression.R
# ============================================================================
#
# This script runs differential expression analysis using DESeq2.
# The goal: figure out which genes are turned UP or DOWN in the translocated
# cell lines (T1 and C1) compared to the wildtype (WT, "normal") cells.
#
# What's differential expression?
# - For each gene, we compare its expression level (how much RNA is made
#   from it) between two groups: WT vs T1, and WT vs C1.
# - If a gene is much more expressed in T1 than WT, it has a positive
#   log2 fold change (it went UP).
# - If it's much less expressed, it has a negative log2 fold change
#   (it went DOWN).
# - DESeq2 also gives us p-values telling us if the change is statistically
#   significant (probably real, not just random noise).
#
# What this script does:
# 1. Load the gene counts (a big table: genes in rows, samples in columns).
# 2. Load the metadata (which sample is which condition).
# 3. Build a DESeq2 dataset and remove genes with very low counts.
# 4. Run DESeq2's standard analysis pipeline.
# 5. Use "apeglm shrinkage" to get cleaner fold change estimates.
#    (Genes with very low counts can have huge but unreliable fold changes;
#    apeglm pulls those toward zero so they don't mislead us.)
# 6. Save the results as TSV files for use in the next script (14).
# 7. Print summaries and make MA plots.
#
# Output files have one row per gene with:
#   gene_id, gene_name, baseMean (avg expression), log2FoldChange,
#   lfcSE (uncertainty in the fold change), svalue, padj (adjusted p-value)
#
# Edit the file paths in the CONFIG section before running.
#
# Required package: DESeq2 (from Bioconductor, not CRAN)
# Install with:
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install("DESeq2")
# ============================================================================


# ============================================================================
# CONFIG SECTION - Edit these paths before running
# ============================================================================

# -----------------------------------------------------------------------------
# Counts file (genes x samples)
# -----------------------------------------------------------------------------
# Tab-separated file with one row per gene and one column per sample,
# plus a "gene_id" column and "gene_name" column.
COUNTS_FILE <- "/path/to/GSE246689_gene_counts.tsv"

# -----------------------------------------------------------------------------
# Sample metadata file
# -----------------------------------------------------------------------------
# Tab-separated file with sample IDs as row names.
# Must have a "condition" column with values: WT, T1, C1
COLDATA_FILE <- "/path/to/coldata.tsv"

# -----------------------------------------------------------------------------
# Output folder for the results
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "/path/to/results/DE"

# -----------------------------------------------------------------------------
# Minimum total counts for a gene to be kept
# -----------------------------------------------------------------------------
# Genes with very few counts across all samples don't have enough data
# to be reliably analyzed, we filter them out.
MIN_COUNT <- 10


# ============================================================================
# Load the libraries we need
# ============================================================================
# DESeq2 is the main package for differential expression analysis.

library(DESeq2)


# ============================================================================
# Step 1: Load the count data
# ============================================================================

cat("Loading count matrix...\n")

# read.table reads a tab-separated file
# header = TRUE means the first row is column names
# sep = "\t" means tab-separated
counts_raw <- read.table(COUNTS_FILE, header = TRUE, sep = "\t")

# -----------------------------------------------------------------------------
# Save a small mapping of gene_id -> gene_name for later
# -----------------------------------------------------------------------------
# We'll need to add gene names back to the results at the end.
# We pull out just the gene_id and gene_name columns now.
# In R, df[, c("col1", "col2")] selects multiple columns.
gene_map <- counts_raw[, c("gene_id", "gene_name")]

# -----------------------------------------------------------------------------
# Set gene_id as the row names of the table
# -----------------------------------------------------------------------------
# This means each row is identified by its gene_id.
# DESeq2 expects this for the count matrix.
rownames(counts_raw) <- counts_raw$gene_id

# -----------------------------------------------------------------------------
# Remove the gene_id and gene_name columns
# -----------------------------------------------------------------------------
# Now counts_raw has gene_id as row names, so we don't need it as a column.
# We also remove gene_name because DESeq2 only wants count numbers.
# In R, df[, -(1:2)] keeps all columns EXCEPT the first 2.
# (The minus sign means "remove these columns".)
# This works because gene_id and gene_name are the first 2 columns.
counts <- counts_raw[, -(1:2)]

# -----------------------------------------------------------------------------
# Make sure all counts are integers
# -----------------------------------------------------------------------------
# DESeq2 requires integer counts. Some count tables might have small
# decimal numbers from prior processing, round them to be safe.
counts <- round(counts)

# Sanity check: counts should never be negative, and they should be whole
# numbers after rounding.
# stopifnot stops the script if the condition is FALSE.
# all() checks that ALL values in a vector are TRUE.
stopifnot(all(counts >= 0))
stopifnot(all(counts == floor(counts)))


# ============================================================================
# Step 2: Load the sample metadata (which sample is which condition)
# ============================================================================

cat("Loading sample metadata...\n")

# row.names = 1 tells R: the first column is row names (sample IDs).
coldata <- read.table(
  COLDATA_FILE,
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

# -----------------------------------------------------------------------------
# Make condition a factor with WT as the reference level
# -----------------------------------------------------------------------------
# DESeq2 will compare other levels TO the reference level.
# We want comparisons like "T1 vs WT" and "C1 vs WT", not the other way
# around. Setting WT first in the levels makes WT the reference.
coldata$condition <- factor(
  coldata$condition,
  levels = c("WT", "T1", "C1")
)

# -----------------------------------------------------------------------------
# Make sure the columns of "counts" are in the same order as the rows of
# "coldata", DESeq2 will complain if they're not aligned
# -----------------------------------------------------------------------------
# rownames(coldata) gives the sample IDs.
# counts[, those_names] reorders the counts columns to match.
counts <- counts[, rownames(coldata)]

# Print a summary of what we've loaded
cat(sprintf(
  "Count matrix: %d genes x %d samples\n",
  nrow(counts),
  ncol(counts)
))

# paste() with collapse = ", " joins a vector into a single string
condition_list <- paste(levels(coldata$condition), collapse = ", ")
cat(sprintf("Conditions: %s\n", condition_list))


# ============================================================================
# Step 3: Build the DESeq2 dataset
# ============================================================================
# A DESeq2 dataset bundles together:
#   - the count matrix
#   - the sample metadata
#   - the "design formula" telling DESeq2 what to compare

cat("\nBuilding DESeq2 dataset...\n")

# The design formula "~ condition" means:
# "model the counts as a function of condition"
# (the ~ is R's notation for "depends on")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ condition
)


# ============================================================================
# Step 4: Filter out very low-count genes
# ============================================================================
# Genes with very few total counts don't have enough data for a reliable
# analysis, we drop them.

# rowSums() adds up the counts in each row (each gene).
# counts(dds) gets the count matrix back out of the DESeq2 dataset.
# We keep only the rows where the row sum is at least MIN_COUNT.
total_counts_per_gene <- rowSums(counts(dds))
keep_these_genes <- (total_counts_per_gene >= MIN_COUNT)
dds <- dds[keep_these_genes, ]

cat(sprintf(
  "Genes after low-count filtering (>= %d counts): %d\n",
  MIN_COUNT,
  nrow(dds)
))


# ============================================================================
# Step 5: Run the DESeq2 analysis
# ============================================================================
# This one function call does a LOT:
# - Calculates "size factors" (normalizes for differences in sequencing depth)
# - Estimates how variable each gene is across samples
# - Fits a statistical model to each gene
# - Tests if each gene is significantly different between conditions

cat("\nRunning DESeq2...\n")
dds <- DESeq(dds)

# -----------------------------------------------------------------------------
# Print the coefficient names
# -----------------------------------------------------------------------------
# DESeq2 fits coefficients for each comparison. The names will be like
# "condition_T1_vs_WT" and "condition_C1_vs_WT".
# We need these names in the next step to extract the right results.
cat("\nCoefficient names:\n")
print(resultsNames(dds))


# ============================================================================
# Step 6: Extract results with apeglm shrinkage
# ============================================================================
# Why shrinkage?
# - For genes with low counts, the raw log2 fold change can be huge but
#   very unreliable. For example, going from 1 count to 4 counts is a 4x
#   change, but with that few counts it could just be random.
# - apeglm "shrinks" these unreliable fold changes toward 0.
# - The end result: more trustworthy fold changes, especially for low-count
#   genes.
#
# We do this twice, once for T1 vs WT, once for C1 vs WT.

cat("\nExtracting results with apeglm shrinkage...\n")

# Get the T1 vs WT results
res_T1 <- lfcShrink(
  dds,
  coef = "condition_T1_vs_WT",
  type = "apeglm"
)

# Get the C1 vs WT results
res_C1 <- lfcShrink(
  dds,
  coef = "condition_C1_vs_WT",
  type = "apeglm"
)


# ============================================================================
# Step 7: Save the results to TSV files
# ============================================================================

# Make sure the output folder exists
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Helper function: convert one DESeq2 result object to a TSV file
# -----------------------------------------------------------------------------
# Both T1 and C1 need the same processing, so we put it in a function.

save_results <- function(res, cond_name) {

  # Convert the DESeq2 results object to a regular data frame
  df <- as.data.frame(res)

  # The gene IDs are currently the row names. Add them as a proper column
  # so they show up in the output file.
  df$gene_id <- rownames(df)

  # Add the gene_name column by merging with our gene_map.
  # merge() matches rows where gene_id is the same in both tables.
  df <- merge(df, gene_map, by = "gene_id")

  # Build the output file path.
  # sprintf with %s inserts a string, here we put cond_name into the filename.
  filename <- sprintf("DE_WT_vs_%s.tsv", cond_name)
  out_path <- file.path(OUTPUT_DIR, filename)

  # Save as TSV.
  # sep = "\t"      = tab-separated
  # quote = FALSE   = don't put quotes around text values
  # row.names = FALSE = don't write the row numbers as a column
  write.table(
    df,
    out_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(sprintf("Saved: %s  (%d genes)\n", out_path, nrow(df)))
}

# Now call the function for each comparison
save_results(res_T1, "T1")
save_results(res_C1, "C1")


# ============================================================================
# Step 8: Print summaries (sanity checks)
# ============================================================================

# -----------------------------------------------------------------------------
# Summary for WT vs T1
# -----------------------------------------------------------------------------
cat("\n--- WT vs T1 summary ---\n")
print(summary(res_T1))

# Count how many genes are significantly differentially expressed.
# padj is the adjusted p-value (corrected for multiple testing).
# A common cutoff is padj < 0.05.
# We use na.rm = TRUE to ignore NA values when summing.
n_sig_T1 <- sum(res_T1$padj < 0.05, na.rm = TRUE)
cat(sprintf("Significant genes (padj < 0.05): %d\n", n_sig_T1))

# -----------------------------------------------------------------------------
# Summary for WT vs C1
# -----------------------------------------------------------------------------
cat("\n--- WT vs C1 summary ---\n")
print(summary(res_C1))

n_sig_C1 <- sum(res_C1$padj < 0.05, na.rm = TRUE)
cat(sprintf("Significant genes (padj < 0.05): %d\n", n_sig_C1))


# ============================================================================
# Step 9: MA plots
# ============================================================================
# An MA plot shows log2 fold change (y-axis) versus mean expression (x-axis).
# - Each dot is one gene
# - Significant genes are highlighted in colour
# - A "good" MA plot should be roughly symmetric around the y=0 line,
#   with most dots near zero (most genes don't change)

plotMA(res_T1, main = "MA plot: WT vs T1")
plotMA(res_C1, main = "MA plot: WT vs C1")

# Done!
cat("\nDone!\n")
