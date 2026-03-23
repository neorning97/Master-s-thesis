library(DESeq2)

# Read counts
counts_raw <- read.table(
  "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/GSE246689_gene_counts.tsv",
  header = TRUE,
  sep = "\t"
)

# Keep gene_id as rownames
rownames(counts_raw) <- counts_raw$gene_id

# Drop annotation columns
counts <- counts_raw[, -(1:2)]    #removes gene_id and gene_name

# Round to integers
counts <- round(counts)

# Sanity check
stopifnot(all(counts >= 0))
stopifnot(all(counts == floor(counts)))

# Read coldata.tsv file
coldata <- read.table(
  "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/coldata.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

coldata$condition <- factor(coldata$condition, levels = c("WT", "T1", "C1"))

# Ensure column order matches
counts <- counts[, rownames(coldata)]

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Filter very low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

# Check coefficient names
resultsNames(dds)

# Shrink log2FCs
res_T1 <- lfcShrink(
  dds,
  coef = "condition_T1_vs_WT",
  type = "apeglm"
)

res_C1 <- lfcShrink(
  dds,
  coef = "condition_C1_vs_WT",
  type = "apeglm"
)

# Save results with gene names
res_T1_df <- as.data.frame(res_T1)
res_C1_df <- as.data.frame(res_C1)

res_T1_df$gene_id <- rownames(res_T1_df)
res_C1_df$gene_id <- rownames(res_C1_df)

# Bring gene_name back in
gene_map <- counts_raw[, c("gene_id", "gene_name")]

res_T1_df <- merge(res_T1_df, gene_map, by = "gene_id")
res_C1_df <- merge(res_C1_df, gene_map, by = "gene_id")

write.table(
  res_T1_df,
  "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/DE_WT_vs_T1.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  res_C1_df,
  "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/DE_WT_vs_C1.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Sanity check
summary(res_T1)
summary(res_C1)

plotMA(res_T1)
plotMA(res_C1)

sum(res_T1$padj < 0.05, na.rm = TRUE)
sum(res_C1$padj < 0.05, na.rm = TRUE)
