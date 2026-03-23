import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.stats import fisher_exact

# -----------------------------
# Paths
# -----------------------------
DIST_FILE = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/WT_translocation_distances_agg.tsv"
DE_T1_FILE = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/DE_WT_vs_T1.tsv"
DE_C1_FILE = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/DE_WT_vs_C1.tsv"

OUTDIR = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/new_plots/DE_enrichment_new_verify"
os.makedirs(OUTDIR, exist_ok=True)

# -----------------------------
# Load data
# -----------------------------
dist = pd.read_csv(DIST_FILE, sep="\t")
de_T1 = pd.read_csv(DE_T1_FILE, sep="\t")
de_C1 = pd.read_csv(DE_C1_FILE, sep="\t")

de_T1 = de_T1.sort_values("padj").drop_duplicates("gene_name", keep="first")
de_C1 = de_C1.sort_values("padj").drop_duplicates("gene_name", keep="first")

# -----------------------------
# Classify DE genes
# -----------------------------
def classify_DE(df, lfc_thresh=1.0, padj_thresh=0.05):
    conditions = [
        (df["log2FoldChange"] > lfc_thresh) & (df["padj"] < padj_thresh),
        (df["log2FoldChange"] < -lfc_thresh) & (df["padj"] < padj_thresh)
    ]
    choices = ["Up", "Down"]
    df["DE_status"] = np.select(conditions, choices, default="No change")
    return df

de_T1 = classify_DE(de_T1)
de_C1 = classify_DE(de_C1)

# -----------------------------
# Identify translocated genes
# -----------------------------
transloc_genes = set(
    dist.loc[
        dist["gene1_region_type"] == "inside_transloc",
        "gene1"
    ].unique()
)

print("Number of translocation genes:", len(transloc_genes))

# -----------------------------
# Count DE categories
# -----------------------------
def get_counts(gene_set, de_df):
    subset = de_df[de_df["gene_name"].isin(gene_set)]
    counts = subset["DE_status"].value_counts().reindex(
        ["Up", "Down", "No change"], fill_value=0
    )
    return counts

# -----------------------------
# Fisher enrichment test
# -----------------------------
def fisher_test(counts_transloc, counts_control, test_type="DE_vs_NoChange"):
    if test_type == "DE_vs_NoChange":
        # DE vs No-change
        de_transloc = counts_transloc["Up"] + counts_transloc["Down"]
        nde_transloc = counts_transloc["No change"]
        de_control = counts_control["Up"] + counts_control["Down"]
        nde_control = counts_control["No change"]
    elif test_type == "Up_vs_Down":
        # Up vs Down (consider only DE genes)
        de_transloc = counts_transloc["Up"]
        nde_transloc = counts_transloc["Down"]
        de_control = counts_control["Up"]
        nde_control = counts_control["Down"]
    else:
        raise ValueError("Invalid test_type. Must be 'DE_vs_NoChange' or 'Up_vs_Down'")

    table = [
        [de_transloc, nde_transloc],
        [de_control, nde_control]
    ]
    oddsratio, p = fisher_exact(table)
    return table, oddsratio, p

# -----------------------------
# Plot grouped barplot
# -----------------------------
def plot_grouped(counts_transloc, counts_control, cond):
    labels = ["Up", "Down", "No change"]
    trans_vals = counts_transloc[labels].values
    control_vals = counts_control[labels].values
    trans_prop = trans_vals / trans_vals.sum()
    control_prop = control_vals / control_vals.sum()

    x = np.arange(len(labels))
    width = 0.35

    plt.figure(figsize=(6,4))
    plt.bar(x - width/2, trans_prop, width, label="Translocation genes", color="steelblue")
    plt.bar(x + width/2, control_prop, width, label="Genome background", color="orange")

    plt.xticks(x, labels)
    plt.ylabel("Fraction of genes")
    plt.title(f"{cond}: DE categories")

    plt.legend()
    plt.tight_layout()

    outfile = os.path.join(OUTDIR, f"{cond}_grouped_DE_proportions.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print("Saved:", outfile)

# -----------------------------
# Run analysis
# -----------------------------
def run_analysis(de_df, cond):
    all_genes = set(de_df["gene_name"])
    trans_genes = all_genes.intersection(transloc_genes)
    control_genes = all_genes - trans_genes

    print(f"\n{cond}")
    print("Total genes:", len(all_genes))
    print("Translocation genes:", len(trans_genes))
    print("Control genes:", len(control_genes))

    counts_transloc = get_counts(trans_genes, de_df)
    counts_control = get_counts(control_genes, de_df)

    print("\nCounts:")
    print("Translocation:", counts_transloc.to_dict())
    print("Control:", counts_control.to_dict())

    # -----------------------------
    # Fisher test: DE vs No-change
    # -----------------------------
    table, OR, p = fisher_test(counts_transloc, counts_control, test_type="DE_vs_NoChange")
    print("\nDE enrichment (Translocation vs Genome):")
    print(table)
    print("Odds ratio:", OR)
    print("P-value:", p)

    # -----------------------------
    # Fisher test: Up vs Down bias
    # -----------------------------
    table_ud, OR_ud, p_ud = fisher_test(counts_transloc, counts_control, test_type="Up_vs_Down")
    print("\nUp vs Down bias (Translocation vs Genome):")
    print(table_ud)
    print("Odds ratio:", OR_ud)
    print("P-value:", p_ud)

    # -----------------------------
    # Plot barplot
    # -----------------------------
    plot_grouped(counts_transloc, counts_control, cond)

# -----------------------------
# Run analyses
# -----------------------------
run_analysis(de_T1, "T1")
run_analysis(de_C1, "C1")

print(de_T1["gene_name"].duplicated().sum())
print(de_C1["gene_name"].duplicated().sum())
