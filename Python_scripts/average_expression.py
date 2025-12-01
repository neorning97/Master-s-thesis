# Script that takes the translocation genes in T1 and C1 (distinguishes the two translocation parts),
# finds the expression data in WT, T1, and C1 for these genes, and
# makes a tsv file that contains these average expressions and the log2 fold change

import pandas as pd
import numpy as np

# -----------------------------------------
# Input files
# -----------------------------------------
t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/T1_genes_expression.tsv"
c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/processed/C1_genes_expression.tsv"
tpm_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/raw/GSE246689_gene_tpm.tsv"

output_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/WT_T1_C1_expression_comparison.tsv"

# -----------------------------------------
# Load T1/C1 translocated genes
# -----------------------------------------
t1 = pd.read_csv(t1_file, sep="\t")
c1 = pd.read_csv(c1_file, sep="\t")

# Keep only genes inside translocation
t1_genes = t1.loc[t1["region_type"] == "inside_transloc", "gene_name"]
c1_genes = c1.loc[c1["region_type"] == "inside_transloc", "gene_name"]

# -----------------------------------------
# Load TPM
# -----------------------------------------
tpm = pd.read_csv(tpm_file, sep="\t")

# -----------------------------------------
# Flag genes in T1/C1 translocations
# -----------------------------------------
tpm["in_T1_transloc"] = tpm["gene_name"].isin(t1_genes)
tpm["in_C1_transloc"] = tpm["gene_name"].isin(c1_genes)

# Keep only genes present in at least one translocation
tpm_filtered = tpm[tpm["in_T1_transloc"] | tpm["in_C1_transloc"]].copy()

# -----------------------------------------
# Compute average expression per condition
# -----------------------------------------
conditions = {
    "WT": "MCF10AWT",
    "T1": "MCF10AT1",
    "C1": "MCF10AC1"
}

for cond_name, prefix in conditions.items():
    rep_cols = [c for c in tpm_filtered.columns if c.startswith(prefix)]
    tpm_filtered[f"{cond_name}_AVG"] = tpm_filtered[rep_cols].mean(axis=1)

# -----------------------------------------
# Compute log2 fold-change vs WT
# -----------------------------------------
tpm_filtered["T1_vs_WT_log2FC"] = np.log2((tpm_filtered["T1_AVG"] + 1) / (tpm_filtered["WT_AVG"] + 1))
tpm_filtered["C1_vs_WT_log2FC"] = np.log2((tpm_filtered["C1_AVG"] + 1) / (tpm_filtered["WT_AVG"] + 1))

# -----------------------------------------
# Keep only relevant columns for output
# -----------------------------------------
output_cols = [
    "gene_name", "in_T1_transloc", "in_C1_transloc",
    "WT_AVG", "T1_AVG", "C1_AVG",
    "T1_vs_WT_log2FC", "C1_vs_WT_log2FC"
]
tpm_filtered[output_cols].to_csv(output_file, sep="\t", index=False)

print(f"Saved translocated gene expression with log2 fold-change and translocation flags to: {output_file}")
