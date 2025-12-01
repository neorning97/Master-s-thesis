# Script that plots comparisons between WT and T1, WT and C1 using the spatial distances
# between the genes in the translocation regions and the center of the nucleus (0,0,0).

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -------------------------------
# Load aggregated distance files
# -------------------------------
wt_df = pd.read_csv("/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/WT_T1C1_genes_center_distances_aggregated.tsv", sep="\t")
wt_df["condition"] = "WT"

t1_df = pd.read_csv("/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/T1_transloc_gene_center_distances_aggregated.tsv", sep="\t")
t1_df["condition"] = "T1"

c1_df = pd.read_csv("/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/C1_transloc_gene_center_distances_aggregated.tsv", sep="\t")
c1_df["condition"] = "C1"

# Keep only relevant columns
wt_df = wt_df[["gene_name", "mean_distance", "condition"]]
t1_df = t1_df[["gene_name", "mean_distance", "condition"]]
c1_df = c1_df[["gene_name", "mean_distance", "condition"]]

# Combine into one DataFrame
df_plot = pd.concat([wt_df, t1_df, c1_df], axis=0)

# -------------------------------
# Boxplot with stripplot (for many points)
# -------------------------------
plt.figure(figsize=(8,6))
sns.boxplot(x="condition", y="mean_distance", data=df_plot, palette="Set2", fliersize=0)  # fliersize=0 hides outliers
sns.stripplot(x="condition", y="mean_distance", data=df_plot, color=".25", jitter=0.3, alpha=0.5)

plt.ylabel("Distance to Nuclear Center")
plt.title("Comparison of Gene Distances to Nuclear Center")
plt.tight_layout()
plt.show()

# -------------------------------
# Scatter plot for paired comparison
# -------------------------------
# Pivot for paired comparison
df_pivot = df_plot.pivot(index="gene_name", columns="condition", values="mean_distance").reset_index()

plt.figure(figsize=(6,6))
plt.scatter(df_pivot["WT"], df_pivot["T1"], label="WT vs T1", alpha=0.5)
plt.scatter(df_pivot["WT"], df_pivot["C1"], label="WT vs C1", alpha=0.5, color="red")

# Diagonal y=x line
min_val = df_pivot[["WT","T1","C1"]].min().min()
max_val = df_pivot[["WT","T1","C1"]].max().max()
plt.plot([min_val, max_val], [min_val, max_val], "k--")

plt.xlabel("WT distance")
plt.ylabel("T1/C1 distance")
plt.legend()
plt.title("Gene distances: WT vs T1/C1")
plt.tight_layout()
plt.show()

# Count genes below and above diagonal for WT vs T1
below_T1 = (df_pivot["T1"] < df_pivot["WT"]).sum()
above_T1 = (df_pivot["T1"] > df_pivot["WT"]).sum()
total_T1 = len(df_pivot)
fraction_below_T1 = below_T1 / total_T1
fraction_above_T1 = above_T1 / total_T1

print(f"WT vs T1: {fraction_below_T1*100:.2f}% below diagonal, {fraction_above_T1*100:.2f}% above diagonal")

# Count genes below and above diagonal for WT vs C1
below_C1 = (df_pivot["C1"] < df_pivot["WT"]).sum()
above_C1 = (df_pivot["C1"] > df_pivot["WT"]).sum()
total_C1 = len(df_pivot)
fraction_below_C1 = below_C1 / total_C1
fraction_above_C1 = above_C1 / total_C1

print(f"WT vs C1: {fraction_below_C1*100:.2f}% below diagonal, {fraction_above_C1*100:.2f}% above diagonal")
