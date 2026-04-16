# ============================================================================
# 18_go_enrichment.R
# ============================================================================
#
# What this script does
# ---------------------
# Runs Gene Ontology (GO) enrichment analysis on the differentially expressed
# translocated genes identified in scripts 16/17, and produces dot plots
# showing the most significantly enriched biological processes.
#
# GO enrichment analysis asks: given a list of genes of interest (here,
# differentially expressed genes in translocated regions), which biological
# processes are represented more often than expected by chance?
#
# For each condition (T1, C1), the script:
# 1. Loads the ranked gene tables produced by script 16 or 17 (up- and
#    downregulated translocated genes separately).
# 2. Combines up- and downregulated genes into a single list per condition.
# 3. Converts HGNC gene symbols to Entrez IDs, which are required by
#    clusterProfiler.
# 4. Runs GO enrichment for Biological Process (BP) ontology using a
#    hypergeometric test with Benjamini-Hochberg multiple testing correction.
# 5. Saves a dot plot showing the top 10 enriched GO terms per condition.
#
# Why combine up and down?
# ------------------------
# Combining up- and downregulated genes gives a broader picture of which
# biological processes are affected by the translocation overall. 
#
# Why Biological Process (BP)?
# ----------------------------
# The three GO ontologies are Biological Process (BP), Molecular Function
# (MF), and Cellular Component (CC). BP is the most commonly used in
# transcriptomic studies because it describes higher-level biological
# pathways and processes, which are more interpretable in the context of
# cancer biology than molecular function or subcellular localisation.
#
# Why BH correction?
# ------------------
# With thousands of GO terms tested simultaneously, many will appear
# significant by chance. Benjamini-Hochberg (BH) correction controls the
# false discovery rate (FDR), ensuring that at most 5% of called enrichments
# are expected to be false positives.
#
# Usage
# -----
#   1. Run scripts 16 or 17 first to generate the input gene tables.
#   2. Edit the CONFIG section below to point to your files.
#   3. Run: Rscript 18_go_enrichment.R
#      or source the script from RStudio.
#
# Dependencies
# ------------
#   clusterProfiler, org.Hs.eg.db, ggplot2  (Bioconductor)
#   Install with:
#     if (!requireNamespace("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#     BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
#     install.packages("ggplot2")
#
# ============================================================================


# ============================================================================
# CONFIG – edit all paths and settings here before running
# ============================================================================

# Input gene tables — TSV files with a 'gene_name' column.
# These are the ranked gene tables exported by scripts 16.
# Point to the up- and downregulated files for each condition.
INPUT_FILES <- list(
  T1_up   = "/path/to/all_up_transloc_genes_WT_vs_T1.tsv",
  T1_down = "/path/to/all_down_transloc_genes_WT_vs_T1.tsv",
  C1_up   = "/path/to/all_up_transloc_genes_WT_vs_C1.tsv",
  C1_down = "/path/to/all_down_transloc_genes_WT_vs_C1.tsv"
)

# Where to save the output dot plots (created automatically)
OUTPUT_DIR <- "/path/to/results/plots"

# Number of top GO terms to show in each dot plot
TOP_N_CATEGORIES <- 10

# GO ontology to use:
#   "BP" = Biological Process (recommended — describes biological pathways)
#   "MF" = Molecular Function (describes what the gene product does)
#   "CC" = Cellular Component (describes where the gene product is located)
GO_ONTOLOGY <- "BP"

# Adjusted p-value threshold for calling a GO term enriched
PADJ_CUTOFF <- 0.05

# ============================================================================


# ============================================================================
# Load libraries
# ============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


# ============================================================================
# Helper function: read gene names from a TSV file
# ============================================================================

read_gene_names <- function(file_path) {
  # Read the file and extract the gene_name column.
  # The input files are the ranked gene tables from scripts 16/17, which
  # always include a 'gene_name' column with HGNC gene symbols.
  df <- read.table(file_path, header = TRUE, sep = "\t")
  return(df$gene_name)
}


# ============================================================================
# Main analysis
# ============================================================================

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load gene lists per condition ─────────────────────────────────────────────
cat("Loading gene lists...\n")

T1_up_genes   <- read_gene_names(INPUT_FILES$T1_up)
T1_down_genes <- read_gene_names(INPUT_FILES$T1_down)
C1_up_genes   <- read_gene_names(INPUT_FILES$C1_up)
C1_down_genes <- read_gene_names(INPUT_FILES$C1_down)

# Combine up- and downregulated genes per condition.
# unique() removes any genes that appear in both lists.
gene_lists <- list(
  T1 = unique(c(T1_up_genes, T1_down_genes)),
  C1 = unique(c(C1_up_genes, C1_down_genes))
)

cat(sprintf("  T1: %d genes total (%d up, %d down)\n",
            length(gene_lists$T1), length(T1_up_genes), length(T1_down_genes)))
cat(sprintf("  C1: %d genes total (%d up, %d down)\n",
            length(gene_lists$C1), length(C1_up_genes), length(C1_down_genes)))

# ── Convert HGNC symbols to Entrez IDs ───────────────────────────────────────
# clusterProfiler requires Entrez IDs rather than gene symbols.
# bitr() performs the conversion using the human annotation database
# (org.Hs.eg.db). Genes that cannot be mapped are automatically dropped
# with a warning — this is normal for a small fraction of gene symbols.
cat("\nConverting gene symbols to Entrez IDs...\n")

entrez_lists <- lapply(names(gene_lists), function(cond) {
  genes   <- gene_lists[[cond]]
  mapped  <- bitr(genes,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Hs.eg.db)
  cat(sprintf("  %s: %d/%d genes successfully mapped to Entrez IDs\n",
              cond, nrow(mapped), length(genes)))
  return(mapped)
})
names(entrez_lists) <- names(gene_lists)


# ── Run GO enrichment analysis ────────────────────────────────────────────────
# enrichGO tests each GO term for over-representation using a hypergeometric
# test. The genome-wide gene set in org.Hs.eg.db serves as the background.
# Results are filtered to terms with padj < PADJ_CUTOFF after BH correction.
cat(sprintf("\nRunning GO enrichment (ontology: %s)...\n", GO_ONTOLOGY))

go_results <- lapply(names(entrez_lists), function(cond) {
  df  <- entrez_lists[[cond]]
  ego <- enrichGO(
    gene          = df$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = GO_ONTOLOGY,
    pAdjustMethod = "BH",
    pvalueCutoff  = PADJ_CUTOFF,
    readable      = TRUE   # converts Entrez IDs back to gene symbols in output
  )
  n_terms <- if (!is.null(ego)) nrow(ego) else 0
  cat(sprintf("  %s: %d enriched GO terms found\n", cond, n_terms))
  return(ego)
})
names(go_results) <- names(entrez_lists)


# ── Save dot plots ────────────────────────────────────────────────────────────
# A dot plot shows the top enriched GO terms on the y-axis.
# The x-axis shows gene ratio (fraction of input genes in the GO term).
# Dot size represents the number of genes, dot colour represents padj.
cat("\nSaving GO dot plots...\n")

for (cond in names(go_results)) {
  ego <- go_results[[cond]]
  
  # Skip if no enriched terms were found
  if (is.null(ego) || nrow(ego) == 0) {
    message(sprintf("  No enriched GO terms for %s — skipping plot.", cond))
    next
  }
  
  p <- dotplot(ego, showCategory = TOP_N_CATEGORIES) +
    ggtitle(sprintf("%s — GO %s enrichment", cond, GO_ONTOLOGY)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  out_path <- file.path(OUTPUT_DIR, sprintf("%s_GO_%s_dotplot.png", cond, GO_ONTOLOGY))
  ggsave(out_path, plot = p, width = 8, height = 6, dpi = 300)
  cat(sprintf("  Saved: %s\n", out_path))
}


# ── Print top results to terminal ─────────────────────────────────────────────
cat("\nTop GO terms per condition:\n")
for (cond in names(go_results)) {
  ego <- go_results[[cond]]
  cat(sprintf("\n%s:\n", cond))
  if (!is.null(ego) && nrow(ego) > 0) {
    print(head(ego@result[, c("Description", "GeneRatio", "p.adjust")],
               n = TOP_N_CATEGORIES))
  } else {
    cat("  No significant enrichment.\n")
  }
}

cat("\nDone.\n")
