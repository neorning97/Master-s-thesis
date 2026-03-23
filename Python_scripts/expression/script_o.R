# ==============================
# GO enrichment plots for DE genes
# ==============================

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# ------------------------------
# Paths to input files
# ------------------------------
input_files <- list(
  T1_up   = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/DE_volcano_1/all_up_transloc_genes_WT_vs_T1.tsv",
  T1_down = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/DE_volcano_1/all_down_transloc_genes_WT_vs_T1.tsv",
  C1_up   = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/DE_volcano_1/all_up_transloc_genes_WT_vs_C1.tsv",
  C1_down = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/DE_volcano_1/all_down_transloc_genes_WT_vs_C1.tsv"
)

# Output folder for plots
outdir <- "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/GO_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------
# Read gene lists, combine UP + DOWN, and convert to Entrez IDs
# ------------------------------

read_genes <- function(file) {
  read.table(file, header = TRUE, sep = "\t")$gene_name
}

T1_up_genes <- read_genes(input_files$T1_up)
T1_down_genes <- read_genes(input_files$T1_down)
C1_up_genes <- read_genes(input_files$C1_up)
C1_down_genes <- read_genes(input_files$C1_down)

combined_lists <- list(
  T1 = unique(c(T1_up_genes, T1_down_genes)),
  C1 = unique(c(C1_up_genes, C1_down_genes))
)

gene_lists <- lapply(combined_lists, function(genes) {
  bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
})

# ------------------------------
# Run GO enrichment (Biological Process)
# ------------------------------
go_results <- lapply(gene_lists, function(df) {
  enrichGO(
    gene = df$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
})

# ------------------------------
# Plot and save GO dotplots
# ------------------------------
for (name in names(go_results)) {
  ego <- go_results[[name]]
  
  # Skip if no enrichment
  if (is.null(ego) || nrow(ego) == 0) {
    message(paste("No enriched GO terms for", name))
    next
  }
  
  p <- dotplot(ego, showCategory = 10) + ggtitle(paste0(name, " - GO BP"))
  
  outfile <- file.path(outdir, paste0(name, "_GO_dotplot.png"))
  ggsave(outfile, plot = p, width = 8, height = 6, dpi = 300)
  message(paste("Saved GO dotplot:", outfile))
}

lapply(go_results, function(x) if (!is.null(x)) head(x@result) else "No results")

