# Script that combines the spatial distance files of T1 and C1 into one file.

import pandas as pd

# Load files
T1 = pd.read_csv("/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/T1_translocation_neighbor_distances_aggregated_filtered.tsv", sep="\t")
C1 = pd.read_csv("/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/C1_translocation_neighbor_distances_aggregated_filtered.tsv", sep="\t")

# Keep only relevant columns
cols = ["gene1","gene2","chr1","chr2","mean_distance"]

T1_small = C1[cols].rename(columns={"mean_distance": "mean_distance_T1"})
C1_small = T1[cols].rename(columns={"mean_distance": "mean_distance_C1"})

# Merge on gene pair + chromosomes
merged = pd.merge(
    T1_small, C1_small,
    on=["gene1","gene2","chr1","chr2"],
    how="inner"
)

# Compute distance change
merged["delta_distance(C1-T1)"] = merged["mean_distance_C1"] - merged["mean_distance_T1"]

merged.to_csv("/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/merged_T1_C1_distances.tsv", sep="\t", index=False)
print("Saved merged file: /Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/merged_T1_C1_distances.tsv")
