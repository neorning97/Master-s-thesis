# Script that compares the T1 and C1 distance data in plots

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -------------------------------------------
# Global distribution of distances (T1 vs C1)
# -------------------------------------------

df = pd.read_csv("/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/results/merged_T1_C1_distances.tsv", sep="\t")

plt.figure(figsize=(10,6))
sns.kdeplot(df["mean_distance_T1"], label="T1", linewidth=2)
sns.kdeplot(df["mean_distance_C1"], label="C1", linewidth=2)

plt.xlabel("Mean distance")
plt.ylabel("Density")
plt.title("Distribution of Gene–Gene Distances (T1 vs C1)")
plt.legend()

plt.tight_layout()
plt.show()

# -------------------------------------------
# Scatterplot of T1 vs C1 distances (per gene pair)
# -------------------------------------------


plt.figure(figsize=(7,7))
sns.scatterplot(
    data=df, 
    x="mean_distance_T1", 
    y="mean_distance_C1",
    alpha=0.3,
    s=20
)

plt.plot([df.mean_distance_T1.min(), df.mean_distance_T1.max()],
         [df.mean_distance_T1.min(), df.mean_distance_T1.max()],
         color="red", linestyle="--")  # diagonal

plt.xlabel("C1 mean distance")
plt.ylabel("T1 mean distance")
plt.title("Pairwise Gene Distances: T1 vs C1")

plt.tight_layout()
plt.show()

# -------------------------------------------
# Top N translocated genes by variance
# -------------------------------------------

df["delta"] = df["mean_distance_C1"] - df["mean_distance_T1"]

plt.figure(figsize=(8,6))
sns.histplot(df["delta"], bins=50, kde=True)
plt.axvline(0, color="black", linestyle="--")
plt.xlabel("Δ distance (C1 - T1)")
plt.ylabel("Count")
plt.title("Distribution of Δ-distance")
plt.show()
