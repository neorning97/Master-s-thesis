import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import wilcoxon

# -------------------------------
# Paths
# -------------------------------
t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_distance_to_center_agg.tsv"
c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_distance_to_center_agg.tsv"
wt_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/WT_distance_to_center_agg.tsv"
#t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_distance_to_center_agg_moved.tsv"
#c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_distance_to_center_agg_moved.tsv"
#wt_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_WT_distance_to_center_agg_moved.tsv"

#output_dir = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/verify_center_distance_plots"
output_dir = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/verify_center_distance_plots_moved"
os.makedirs(output_dir, exist_ok=True)

# -------------------------------
# Load data
# -------------------------------
t1 = pd.read_csv(t1_file, sep="\t").rename(columns={"mean_distance":"mean_T1"})
c1 = pd.read_csv(c1_file, sep="\t").rename(columns={"mean_distance":"mean_C1"})
wt = pd.read_csv(wt_file, sep="\t").rename(columns={"mean_distance":"mean_WT"})

# -------------------------------
# Plot function
# -------------------------------
def plot_wt_vs_cond_by_region(wt_df, cond_df, cond_col, label):
    merged = wt_df.merge(
        cond_df,
        on=["gene_name","region_type"],
        how="inner"
    )

    if merged.empty:
        print("No data for", label)
        return

    for region in ["inside_transloc", "control"]:
        sub = merged[merged["region_type"] == region].copy()
        if sub.empty:
            continue

        tag = f"{label}_{region}"

        # ---------------- Fraction closer to condition ----------------
        closer_mask = sub[cond_col] < sub["mean_WT"]
        frac_closer = closer_mask.mean()   # fraction
        n_genes = len(sub)

        print(f"{tag}: {frac_closer:.3f} ({closer_mask.sum()}/{n_genes}) genes closer to {label}")

        # ---------------- Stats ----------------
        stat, p = wilcoxon(sub["mean_WT"], sub[cond_col])
        print(f"{tag}: Wilcoxon p = {p:.3e}")

        # ---------------- Scatter ----------------
        plt.figure(figsize=(6,6))
        sns.scatterplot(
            data=sub,
            x="mean_WT",
            y=cond_col
        )

        lims = [
            sub[["mean_WT", cond_col]].min().min(),
            sub[["mean_WT", cond_col]].max().max()
        ]
        plt.plot(lims, lims, "--r")

        plt.xlabel("WT distance to center")
        plt.ylabel(f"{label} distance to center")
        plt.title(
            f"{label} vs WT ({region})\n"
            f"Fraction closer to {label}: {frac_closer:.2f}"
        )
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{tag}_scatter.png"))
        plt.close()

        # ---------------- Boxplot ----------------
        diff = sub["mean_WT"] - sub[cond_col]
        plt.figure(figsize=(4,5))
        sns.boxplot(y=diff)
        plt.ylabel("WT − condition distance")
        plt.title(
            f"{label} − WT ({region})\n"
            f"Fraction closer to {label}: {frac_closer:.2f}"
        )
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{tag}_boxplot.png"))
        plt.close()

        # ---------------- Density ----------------
        plt.figure(figsize=(5,5))
        sns.kdeplot(sub["mean_WT"], label="WT", fill=True)
        sns.kdeplot(sub[cond_col], label=label, fill=True)
        plt.xlabel("Distance to nuclear center")
        plt.title(
            f"Distance to center ({region})\n"
            f"Fraction closer to {label}: {frac_closer:.2f}"
        )
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{tag}_density.png"))
        plt.close()

# -------------------------------
# Run plots
# -------------------------------
plot_wt_vs_cond_by_region(wt, t1, "mean_T1", "T1")
plot_wt_vs_cond_by_region(wt, c1, "mean_C1", "C1")

print("All center-distance plots saved to:", output_dir)
