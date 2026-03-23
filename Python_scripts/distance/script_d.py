import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from scipy.stats import wilcoxon


# -------------------------------
# File paths
# -------------------------------
c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_translocation_distances_agg.tsv"
t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_translocation_distances_agg.tsv"
wt_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_WT_translocation_distances_agg.tsv"
output_dir = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/distance_plots15Mb_new"
#c1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/C1_translocation_distances_agg_moved.tsv"
#t1_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/T1_translocation_distances_agg_moved.tsv"
#wt_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/WT_translocation_distances_agg_moved.tsv"
#output_dir = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/distance_plots_moved15Mb_test"
os.makedirs(output_dir, exist_ok=True)




# -------------------------------
# Load data
# -------------------------------
c1_df = pd.read_csv(c1_file, sep="\t")
t1_df = pd.read_csv(t1_file, sep="\t")
wt_df = pd.read_csv(wt_file, sep="\t")

# Ensure all chr columns are strings with "chr" prefix
for df in [wt_df, t1_df, c1_df]:
    for col in ["chr1", "chr2"]:
        df[col] = df[col].astype(str)        # convert to string
        df[col] = df[col].apply(lambda x: x if x.startswith("chr") else "chr" + x)

# -------------------------------
# Count unique genes in each category
# -------------------------------
def count_gene_types(df, label):
    trans = set(df.loc[df["gene1_region_type"]=="inside_transloc","gene1"]) | \
            set(df.loc[df["gene2_region_type"]=="inside_transloc","gene2"])

    neighbor = set(df.loc[df["gene1_region_type"]=="neighbor","gene1"]) | \
               set(df.loc[df["gene2_region_type"]=="neighbor","gene2"])

    control = set(df.loc[df["gene1_region_type"]=="control","gene1"]) | \
              set(df.loc[df["gene2_region_type"]=="control","gene2"])

    all_genes = trans | neighbor | control

    print(f"\nGene counts in {label}:")
    print(f"  inside_transloc genes: {len(trans)}")
    print(f"  neighbor genes:       {len(neighbor)}")
    print(f"  control genes:        {len(control)}")
    print(f"  total unique genes:   {len(all_genes)}")


count_gene_types(t1_df, "T1")
count_gene_types(c1_df, "C1")
count_gene_types(wt_df, "WT")

# -------------------------------
# Keep relevant columns
# -------------------------------
cols = ["gene1","gene1_region_type","chr1","gene2","gene2_region_type","chr2","mean_distance","std_distance"]
t1_df = t1_df[cols].rename(columns={"mean_distance":"mean_T1","std_distance":"std_T1"})
c1_df = c1_df[cols].rename(columns={"mean_distance":"mean_C1","std_distance":"std_C1"})
wt_df = wt_df[cols].rename(columns={"mean_distance":"mean_WT","std_distance":"std_WT"})

# -------------------------------
# Split WT into transloc-neighbor and control-neighbor pairs
# -------------------------------
wt_transloc = wt_df[
    ((wt_df["gene1_region_type"]=="inside_transloc") & (wt_df["gene2_region_type"]=="neighbor")) |
    ((wt_df["gene1_region_type"]=="neighbor") & (wt_df["gene2_region_type"]=="inside_transloc"))
]

wt_control = wt_df[
    ((wt_df["gene1_region_type"]=="control") & (wt_df["gene2_region_type"]=="neighbor")) |
    ((wt_df["gene1_region_type"]=="neighbor") & (wt_df["gene2_region_type"]=="control"))
]

# -------------------------------
# Density plot function
# -------------------------------
def plot_density(merged, cond_col, plot_name):
    if merged.empty:
        return

    plt.figure(figsize=(8,6))

    # WT distribution
    sns.kdeplot(
        merged["mean_WT"],
        label="WT",
        fill=True,
        alpha=0.4
    )

    # Condition distribution
    sns.kdeplot(
        merged[cond_col],
        label=plot_name,
        fill=True,
        alpha=0.4
    )

    plt.xlabel("Mean distance")
    plt.ylabel("Density")
    plt.title(f"Distance distribution: WT vs {plot_name}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{plot_name}_density.png"))
    plt.close()

def enforce_equal_pairs(wt_df, cond_df):
        pair_cols = ["gene1","gene2","gene1_region_type","gene2_region_type","chr1","chr2"]

        wt_pairs = wt_df[pair_cols].drop_duplicates()
        cond_pairs = cond_df[pair_cols].drop_duplicates()

        common_pairs = wt_pairs.merge(cond_pairs, on=pair_cols, how="inner")

        wt_eq = wt_df.merge(common_pairs, on=pair_cols, how="inner")
        cond_eq = cond_df.merge(common_pairs, on=pair_cols, how="inner")

        return wt_eq, cond_eq
# -------------------------------
# Function to plot WT vs condition
# -------------------------------
def plot_wt_vs_cond(wt_sub, cond_df, cond_col, plot_name):
    cond_sub = cond_df[
        (cond_df["gene1_region_type"].isin(wt_sub["gene1_region_type"].unique())) &
        (cond_df["gene2_region_type"].isin(wt_sub["gene2_region_type"].unique()))
    ]
    wt_eq, cond_eq = enforce_equal_pairs(wt_sub, cond_sub)

    merged = wt_eq.merge(
        cond_eq,
        on=["gene1","gene2","gene1_region_type","gene2_region_type","chr1","chr2"],
        how="inner"
    )
    
    if merged.empty:
        print(f"No pairs to plot for {plot_name}")
        return merged
    
    # Scatter plot
    plt.figure(figsize=(8,6))
    sns.scatterplot(
        data=merged,
        x="mean_WT",
        y=cond_col,
        hue="gene1_region_type",
        style="gene2_region_type",
        palette={"inside_transloc":"blue","control":"green","neighbor":"orange"}
    )
    plt.plot([merged["mean_WT"].min(), merged["mean_WT"].max()],
             [merged["mean_WT"].min(), merged["mean_WT"].max()],
             color='red', linestyle='--', label='y=x')
    plt.xlabel("WT mean distance")
    plt.ylabel(f"{plot_name} mean distance")
    plt.title(f"WT vs {plot_name}")
    plt.legend(title="Gene type")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{plot_name}_scatter.png"))
    plt.close()
    
    # Boxplot
    plt.figure(figsize=(8,6))
    sns.boxplot(
        x=merged["mean_WT"] - merged[cond_col],
        hue=merged["gene1_region_type"],
        palette={"inside_transloc":"blue","control":"green","neighbor":"orange"}
    )
    plt.ylabel("WT - condition distance")
    plt.title(f"Distance differences: WT - {plot_name}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{plot_name}_boxplot.png"))
    plt.close()

    # Density plot
    plot_density(merged, cond_col, plot_name)



    # Wilcoxon test
    stat, p = wilcoxon(merged["mean_WT"], merged[cond_col])
    print(f"{plot_name} Wilcoxon: stat={stat:.2f}, p={p:.4g}")
    diff = merged["mean_WT"] - merged[cond_col]
    print(diff.describe())
    print((diff==0).sum())

    # Fraction below y=x line
    num_below = (merged[cond_col] < merged["mean_WT"]).sum()
    total_pairs = len(merged)
    fraction_below = num_below / total_pairs
    print(f"Fraction of gene pairs below y=x ({plot_name} closer than WT): {fraction_below:.3f} ({num_below}/{total_pairs})")
    
    return merged



# -------------------------------
# Plot T1 vs WT
# -------------------------------
print("T1 vs WT:")
merged_T1_transloc = plot_wt_vs_cond(wt_transloc, t1_df, "mean_T1", "T1_transloc")
merged_T1_control = plot_wt_vs_cond(wt_control, t1_df, "mean_T1", "T1_control")

# -------------------------------
# Plot C1 vs WT
# -------------------------------
print("C1 vs WT:")
merged_C1_transloc = plot_wt_vs_cond(wt_transloc, c1_df, "mean_C1", "C1_transloc")
merged_C1_control = plot_wt_vs_cond(wt_control, c1_df, "mean_C1", "C1_control")

print("All plots saved to:", output_dir)
