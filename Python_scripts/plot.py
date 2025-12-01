# Script that plots the top 20 genes expressed in the T1 and C1 translocations,
# and compares T1 and C1 to WT

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --------------------------------------------
# File paths
# --------------------------------------------
input_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/WT_T1_C1_expression_comparison.tsv"
output_prefix = "top_genes"

# --------------------------------------------
# 1. Load data
# --------------------------------------------
df = pd.read_csv(input_file, sep="\t")
print(f"Loaded {len(df)} genes.")

top_n = 20  # number of top genes per condition

# --------------------------------------------
# 2. Select top genes per condition, filtered by translocation
# --------------------------------------------
top_wt = (
    df.sort_values('WT_AVG', ascending=False)
      .drop_duplicates(subset='gene_name', keep='first')
      .head(top_n)
)

top_t1 = (
    df[df['in_T1_transloc']]
      .sort_values('T1_AVG', ascending=False)
      .drop_duplicates(subset='gene_name', keep='first')
      .head(top_n)
)

top_c1 = (
    df[df['in_C1_transloc']]
      .sort_values('C1_AVG', ascending=False)
      .drop_duplicates(subset='gene_name', keep='first')
      .head(top_n)
)

# Combine all gene_names
top_genes_names = pd.concat([top_wt, top_t1, top_c1])['gene_name'].unique()
top_genes = df[df['gene_name'].isin(top_genes_names)].copy()

# ---------------------------------------------
# 3. Bar plots for each condition
# ---------------------------------------------
os.makedirs("/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/plots", exist_ok=True)

for condition, top_df in zip(['WT_AVG', 'T1_AVG', 'C1_AVG'], [top_wt, top_t1, top_c1]):
    plt.figure(figsize=(10, 6))
    data = top_df.sort_values(condition, ascending=False)
    plt.barh(data['gene_name'], data[condition])
    plt.xlabel("Expression (TPM)")
    plt.ylabel("Gene")
    plt.title(f"Top {top_n} genes by {condition}")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(
        f"/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/plots/{output_prefix}_{condition}_barplot.png",
        dpi=300
    )
    plt.show()

# --------------------------------------------
# 4. Heatmaps: log2 fold change vs WT
# --------------------------------------------
def plot_log2fc_heatmap(top_genes, fc_column, condition_name):
    heatmap_data = top_genes.set_index('gene_name')[[fc_column]]
    plt.figure(figsize=(6, len(top_genes) * 0.5 + 2))
    sns.heatmap(heatmap_data, cmap="vlag", center=0, annot=True, fmt=".2f", yticklabels=True)
    plt.yticks(rotation=0)
    plt.title(f"Top {top_n} {condition_name} genes (log2 fold change) vs WT")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.savefig(
        f"/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/plots/{output_prefix}_{condition_name}_log2FC_heatmap.png",
        dpi=300
    )
    plt.show()

plot_log2fc_heatmap(top_t1, 'T1_vs_WT_log2FC', 'T1')
plot_log2fc_heatmap(top_c1, 'C1_vs_WT_log2FC', 'C1')

# --------------------------------------------------
# 5. Heatmap: Top 20 WT genes across WT, T1, and C1
# --------------------------------------------------
top20_wt_genes = top_wt
heatmap_data = top20_wt_genes.set_index('gene_name')[['WT_AVG', 'T1_AVG', 'C1_AVG']]
heatmap_data_log2 = np.log2(heatmap_data + 1)  # log2 scale

plt.figure(figsize=(8, 8))
sns.heatmap(
    heatmap_data_log2,
    cmap="vlag",
    center=heatmap_data_log2['WT_AVG'].mean(),
    annot=True,
    fmt=".2f"
)
plt.title("Top 20 WT Genes: Expression across WT, T1, and C1 (log2 scale)")
plt.xlabel("Condition")
plt.ylabel("Gene")
plt.tight_layout()
plt.savefig(
    "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/plots/top20_WT_expression_heatmap.png",
    dpi=300
)
plt.show()
