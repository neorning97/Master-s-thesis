# ============================================================================
# 18_go_enrichment.R
# ============================================================================
#
# This script runs Gene Ontology (GO) enrichment analysis on the
# differentially expressed translocated genes.
#
# What's GO enrichment?
# - "GO" = Gene Ontology, a giant database that labels every gene with
#   the biological processes it's involved in (e.g. "cell division",
#   "immune response").
# - Given a list of "interesting" genes (here, the ones differentially
#   expressed in translocated regions), enrichment analysis asks:
#   "Are any biological processes showing up MORE often in this list
#   than we'd expect by chance?"
# - If yes, that suggests the translocation specifically affects those
#   processes.
#
# What this script does:
# 1. Loads the up- and down-regulated translocated gene lists (from
#    script 16 or 17).
# 2. Combines up + down into one big list per condition.
# 3. Converts gene names (like "ABC1") to Entrez IDs (like "12345"),
#    because the GO analysis package needs Entrez IDs.
# 4. Runs GO enrichment for "Biological Process" (BP), the most
#    commonly used type.
# 5. Saves a dot plot showing the top enriched GO terms.
#
# Why combine up and down into one list?
# - The hypothesis is that the translocation broadly disrupts certain
#   biological processes. We don't really care if a gene went up or down,
#   as long as it changed.
#
# What's BP?
# - GO has 3 ontologies: BP (Biological Process), MF (Molecular Function),
#   CC (Cellular Component).
# - BP = "what biological process is this gene involved in?"
# - MF = "what does the protein do at a molecular level?"
# - CC = "where in the cell is the protein found?"
# - BP is usually the most informative for cancer biology questions.
#
# What's BH correction?
# - When testing thousands of GO terms, some will look "significant" just
#   by chance (false positives).
# - BH (Benjamini-Hochberg) correction adjusts the p-values so we control
#   the false discovery rate (FDR) at 5%.
# - This means: of all the GO terms we call enriched, we expect at most
#   5% to be false positives.
#
# Run scripts 16 or 17 first to make the gene lists.
# Then edit the file paths and run this script.
#
# Required packages: clusterProfiler, org.Hs.eg.db, ggplot2 (Bioconductor)
# Install with:
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
#   install.packages("ggplot2")
# ============================================================================


# ============================================================================
# CONFIG SECTION - Edit these settings before running
# ============================================================================

# -----------------------------------------------------------------------------
# Input gene lists (from scripts 16 or 17)
# -----------------------------------------------------------------------------
# These TSV files must have a "gene_name" column.
# We need 4 files: up + down for both T1 and C1.
T1_UP_FILE   <- "/path/to/all_up_transloc_genes_WT_vs_T1.tsv"
T1_DOWN_FILE <- "/path/to/all_down_transloc_genes_WT_vs_T1.tsv"
C1_UP_FILE   <- "/path/to/all_up_transloc_genes_WT_vs_C1.tsv"
C1_DOWN_FILE <- "/path/to/all_down_transloc_genes_WT_vs_C1.tsv"

# -----------------------------------------------------------------------------
# Output folder for the dot plots
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "/path/to/results/plots"

# -----------------------------------------------------------------------------
# How many top GO terms to show in each plot
# -----------------------------------------------------------------------------
TOP_N_CATEGORIES <- 10

# -----------------------------------------------------------------------------
# Which GO ontology to use
# -----------------------------------------------------------------------------
# Options:
#   "BP" = Biological Process (most common for cancer biology)
#   "MF" = Molecular Function
#   "CC" = Cellular Component
GO_ONTOLOGY <- "BP"

# -----------------------------------------------------------------------------
# Adjusted p-value cutoff for calling a GO term enriched
# -----------------------------------------------------------------------------
# Standard cutoff is 0.05 (= 5% false discovery rate)
PADJ_CUTOFF <- 0.05


# ============================================================================
# Load libraries
# ============================================================================
# clusterProfiler does the GO enrichment analysis.
# org.Hs.eg.db is the human gene annotation database (needed to look up
# what GO terms each gene belongs to).
# ggplot2 is for making plots.

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


# ============================================================================
# Helper function: read gene names from a TSV file
# ============================================================================
# All four input files have the same format, a "gene_name" column.
# This little helper just reads the file and pulls out that column.

read_gene_names <- function(file_path) {

  # Read the tab-separated file
  df <- read.table(file_path, header = TRUE, sep = "\t")

  # Return just the gene_name column (a vector of names)
  return(df$gene_name)
}


# ============================================================================
# Step 1: Make sure the output folder exists
# ============================================================================

# showWarnings = FALSE means: don't warn if the folder already exists
# recursive = TRUE means: also create parent folders if needed
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# Step 2: Load the gene lists
# ============================================================================

cat("Loading gene lists...\n")

# Read each file
T1_up_genes   <- read_gene_names(T1_UP_FILE)
T1_down_genes <- read_gene_names(T1_DOWN_FILE)
C1_up_genes   <- read_gene_names(C1_UP_FILE)
C1_down_genes <- read_gene_names(C1_DOWN_FILE)

# -----------------------------------------------------------------------------
# Combine up + down genes per condition
# -----------------------------------------------------------------------------
# c() concatenates vectors (like Python's list +).
# unique() removes duplicates (in case some genes appear in both up and down).

T1_all_genes <- unique(c(T1_up_genes, T1_down_genes))
C1_all_genes <- unique(c(C1_up_genes, C1_down_genes))

# Group them in a list so we can loop through them later.
# A "list" in R is like a Python dict, it can hold named values.
gene_lists <- list(
  T1 = T1_all_genes,
  C1 = C1_all_genes
)

# Print a summary
# sprintf is like Python's format, %d is an integer placeholder.
cat(sprintf(
  "  T1: %d genes total (%d up, %d down)\n",
  length(T1_all_genes),
  length(T1_up_genes),
  length(T1_down_genes)
))
cat(sprintf(
  "  C1: %d genes total (%d up, %d down)\n",
  length(C1_all_genes),
  length(C1_up_genes),
  length(C1_down_genes)
))


# ============================================================================
# Step 3: Convert gene SYMBOLS to Entrez IDs
# ============================================================================
# Why?
# - Our input gene lists use HGNC symbols (like "TP53", "ABC1").
# - clusterProfiler's GO functions need Entrez IDs (like "7157", "999").
# - bitr() ("Biological ID Translator") does this conversion using the
#   org.Hs.eg.db database.
# - It's normal for some genes to not be mapped, bitr drops them and
#   prints a warning.

cat("\nConverting gene symbols to Entrez IDs...\n")

# Empty list, we'll fill it with one mapping table per condition
entrez_lists <- list()

# names(gene_lists) gives us the keys: c("T1", "C1")
for (cond in names(gene_lists)) {

  # Get the gene names for this condition
  genes <- gene_lists[[cond]]

  # Convert symbols to Entrez IDs
  # bitr returns a data frame with two columns: SYMBOL and ENTREZID
  mapped <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  # Print how many were successfully mapped
  cat(sprintf(
    "  %s: %d/%d genes successfully mapped to Entrez IDs\n",
    cond,
    nrow(mapped),
    length(genes)
  ))

  # Store the mapping in our list
  entrez_lists[[cond]] <- mapped
}


# ============================================================================
# Step 4: Run GO enrichment analysis
# ============================================================================
# enrichGO does the actual statistical test:
# - For each GO term, it counts how many of OUR genes are in that term.
# - It compares that to how many genes from the GENOME are in that term.
# - It uses a hypergeometric test to ask: "Is our count higher than expected
#   by chance?"
# - It then applies BH correction across all GO terms tested.

cat(sprintf("\nRunning GO enrichment (ontology: %s)...\n", GO_ONTOLOGY))

# Empty list, we'll fill it with one enrichment result per condition
go_results <- list()

for (cond in names(entrez_lists)) {

  # Get the Entrez IDs for this condition
  df <- entrez_lists[[cond]]

  # Run the enrichment analysis
  ego <- enrichGO(
    gene = df$ENTREZID,        # our genes of interest
    OrgDb = org.Hs.eg.db,      # the annotation database
    ont = GO_ONTOLOGY,         # which ontology (BP, MF, or CC)
    pAdjustMethod = "BH",      # multiple testing correction method
    pvalueCutoff = PADJ_CUTOFF, # only keep terms below this padj
    readable = TRUE            # convert Entrez IDs back to gene names
                               # in the output (much easier to read)
  )

  # ego could be NULL if no enrichment was found.
  # The if/else here handles both cases.
  if (is.null(ego)) {
    n_terms <- 0
  } else {
    n_terms <- nrow(ego)
  }

  cat(sprintf("  %s: %d enriched GO terms found\n", cond, n_terms))

  # Store the result
  go_results[[cond]] <- ego
}


# ============================================================================
# Step 5: Make and save dot plots
# ============================================================================
# A dot plot shows GO terms on the y-axis.
# - X-axis: "gene ratio" (what fraction of our input genes are in this term)
# - Dot size: how many of our genes are in this term
# - Dot colour: the adjusted p-value (more significant = different colour)

cat("\nSaving GO dot plots...\n")

# Loop through each condition's results
for (cond in names(go_results)) {

  ego <- go_results[[cond]]

  # ---------------------------------------------------------------------------
  # Skip if there are no enriched terms
  # ---------------------------------------------------------------------------
  # message() prints to stderr (the console), similar to cat() but
  # marked as a "message" not regular output.
  if (is.null(ego) || nrow(ego) == 0) {
    message(sprintf("  No enriched GO terms for %s — skipping plot.", cond))
    next   # skip to the next iteration (like Python's continue)
  }

  # ---------------------------------------------------------------------------
  # Make the dot plot
  # ---------------------------------------------------------------------------
  # dotplot() is from clusterProfiler, it knows how to make a nice
  # dot plot from an enrichGO result object.
  # showCategory = how many top terms to show.
  # We add a title and centre it using ggplot2's theme() function.
  # The + here is ggplot2's way of adding layers to a plot.
  p <- dotplot(ego, showCategory = TOP_N_CATEGORIES) +
    ggtitle(sprintf("%s — GO %s enrichment", cond, GO_ONTOLOGY)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  # ---------------------------------------------------------------------------
  # Build the output path and save
  # ---------------------------------------------------------------------------
  filename <- sprintf("%s_GO_%s_dotplot.png", cond, GO_ONTOLOGY)
  out_path <- file.path(OUTPUT_DIR, filename)

  # ggsave saves a ggplot to a file
  ggsave(
    out_path,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )

  cat(sprintf("  Saved: %s\n", out_path))
}


# ============================================================================
# Step 6: Print the top GO terms to the terminal
# ============================================================================
# This gives us a quick text view of the results without having to open
# the dot plots.

cat("\nTop GO terms per condition:\n")

for (cond in names(go_results)) {

  ego <- go_results[[cond]]

  cat(sprintf("\n%s:\n", cond))

  # If there's no result, just say so
  if (is.null(ego) || nrow(ego) == 0) {
    cat("  No significant enrichment.\n")
    next
  }

  # Show only the columns we care about: Description, GeneRatio, p.adjust.
  # ego@result is the data frame inside the ego object.
  # The @ syntax is R's way of accessing slots in S4 objects (a special
  # type of object).
  result_subset <- ego@result[, c("Description", "GeneRatio", "p.adjust")]

  # head() returns the first n rows
  print(head(result_subset, n = TOP_N_CATEGORIES))
}

# Done!
cat("\nDone!\n")
