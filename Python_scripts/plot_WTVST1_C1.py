# Script that takes the distance files and plots the data for analysis

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import scipy.stats as stats
from scipy.stats import wilcoxon



# -------------------------------
# File paths
# -------------------------------
c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/C1_translocation_neighbor_distances_aggregated_filtered.tsv"
t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/T1_translocation_neighbor_distances_aggregated_filtered.tsv"
wt_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/WT_translocation_neighbor_distances_aggregated.tsv"
output_dir = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/plots/distance_plots"
os.makedirs(output_dir, exist_ok=True)

# -------------------------------
# Load data
# -------------------------------
c1_df = pd.read_csv(c1_file, sep="\t")
t1_df = pd.read_csv(t1_file, sep="\t")
wt_df = pd.read_csv(wt_file, sep="\t")

# Keep only relevant columns
cols = ["gene1", "gene2", "chr1", "chr2", "mean_distance", "std_distance"]
c1_df = c1_df[cols].rename(columns={"mean_distance": "mean_C1", "std_distance": "std_C1", "chr1": "gene1_chr", "chr2": "gene2_chr"})
t1_df = t1_df[cols].rename(columns={"mean_distance": "mean_T1", "std_distance": "std_T1", "chr1": "gene1_chr", "chr2": "gene2_chr"})
wt_df = wt_df[cols].rename(columns={"mean_distance": "mean_WT", "std_distance": "std_WT", "chr1": "gene1_chr", "chr2": "gene2_chr"})

# -------------------------------
# Merge all on gene pairs
# -------------------------------
merged_df = wt_df.merge(t1_df, on=["gene1","gene2"], how="inner")
merged_df = merged_df.merge(c1_df, on=["gene1","gene2"], how="inner")

print(f"Number of gene pairs in merged data: {len(merged_df)}")

# -------------------------------
# Scatter plots: WT vs T1 and WT vs C1
# -------------------------------
sns.set(style="whitegrid", font_scale=1.2)

plt.figure(figsize=(8,6))
sns.scatterplot(data=merged_df, x="mean_WT", y="mean_T1")
plt.plot([merged_df["mean_WT"].min(), merged_df["mean_WT"].max()],
         [merged_df["mean_WT"].min(), merged_df["mean_WT"].max()],
         color='red', linestyle='--', label='y=x')
plt.xlabel("WT mean distance")
plt.ylabel("T1 mean distance")
plt.title("WT vs T1 gene pair distances")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "WT_vs_T1_scatter.png"))
plt.close()

plt.figure(figsize=(8,6))
sns.scatterplot(data=merged_df, x="mean_WT", y="mean_C1")
plt.plot([merged_df["mean_WT"].min(), merged_df["mean_WT"].max()],
         [merged_df["mean_WT"].min(), merged_df["mean_WT"].max()],
         color='red', linestyle='--', label='y=x')
plt.xlabel("WT mean distance")
plt.ylabel("C1 mean distance")
plt.title("WT vs C1 gene pair distances")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "WT_vs_C1_scatter.png"))
plt.close()

# -------------------------------
# Optional: Boxplots of distance changes
# -------------------------------
distance_diff_df = pd.DataFrame({
    "WT-T1": merged_df["mean_WT"] - merged_df["mean_T1"],
    "WT-C1": merged_df["mean_WT"] - merged_df["mean_C1"]
})

plt.figure(figsize=(8,6))
sns.boxplot(data=distance_diff_df)
plt.ylabel("Distance difference (WT - T1/C1)")
plt.title("Distribution of gene pair distance changes")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "distance_differences_boxplot.png"))
plt.close()

print("Plots saved to folder:", output_dir)

# --------------------------------
# Gene pairs that move closer in WT
# --------------------------------

# Merge on gene pairs (gene1 and gene2)
merged = pd.merge(
    wt_df,
    c1_df,
    on=["gene1", "gene2"],
    suffixes=("_WT", "_C1")
)

# Filter for gene pairs that are closer in WT
closer_in_WT = merged[merged["mean_WT"] < merged["mean_C1"]]

# Optional: save to a file
closer_in_WT.to_csv("/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/gene_pairs_closer_in_WT_than_C1.tsv", sep="\t", index=False)

print(f"Found {closer_in_WT.shape[0]} gene pairs closer in WT than C1")
print(closer_in_WT[["gene1","gene2","mean_WT","mean_C1"]].head(10))

# -----------------------------------------------------------------------
# Fraction of gene pairs that move closer or further away in T1/C1 vs WT
# -----------------------------------------------------------------------

closer_in_T1 = (merged_df["mean_WT"] > merged_df["mean_T1"]).mean()
further_in_T1 = (merged_df["mean_WT"] < merged_df["mean_T1"]).mean()
closer_in_C1 = (merged_df["mean_WT"] > merged_df["mean_C1"]).mean()
further_in_C1 = (merged_df["mean_WT"] < merged_df["mean_C1"]).mean()
print(f"{closer_in_T1*100:.2f}% of pairs are closer in T1")
print(f"{closer_in_C1*100:.2f}% of pairs are closer in C1")


# ------------------------------------------
# Check if distance changes are significant
# ------------------------------------------
# Define difference arrays for statistical tests
diff_C1 = merged_df["mean_WT"] - merged_df["mean_C1"]
diff_T1 = merged_df["mean_WT"] - merged_df["mean_T1"]

# First check if normally distributed:
plt.hist(diff_C1, bins=100)
plt.xlabel("WT - C1 distance difference")
plt.ylabel("Frequency")
plt.title("Distribution of distance differences")
plt.show()

plt.hist(diff_T1, bins=100)
plt.xlabel("WT - C1 distance difference")
plt.ylabel("Frequency")
plt.title("Distribution of distance differences")
plt.show()

stats.probplot(diff_C1, dist="norm", plot=plt)
plt.title("Q-Q plot of WT - C1 distance differences")
plt.show()
stats.probplot(diff_T1, dist="norm", plot=plt)
plt.title("Q-Q plot of WT - T1 distance differences")
plt.show()

# WT vs C1
stat_C1, p_C1 = wilcoxon(merged_df["mean_WT"], merged_df["mean_C1"])

# WT vs T1
stat_T1, p_T1 = wilcoxon(merged_df["mean_WT"], merged_df["mean_T1"])

print(f"WT vs C1 Wilcoxon signed-rank test: stat={stat_C1:.2f}, p={p_C1:.20f}")
print(f"WT vs T1 Wilcoxon signed-rank test: stat={stat_T1:.2f}, p={p_T1:.20f}")
