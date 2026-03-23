import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# -----------------------------
# Paths
# -----------------------------
DE_T1_FILE = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/DE_WT_vs_T1.tsv"
DE_C1_FILE = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/DE_WT_vs_C1.tsv"
DIST_T1_FILE = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_T1_translocation_distances_agg.tsv"
DIST_C1_FILE = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/verify_C1_translocation_distances_agg.tsv"
OUTDIR = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/DE_volcano_verify"
os.makedirs(OUTDIR, exist_ok=True)

# -----------------------------
# Volcano plot function
# -----------------------------
def plot_volcano(df, dist_file, cond_name, 
                 lfc_thresh=1.0, padj_thresh=0.05, 
                 top_n_table=20, top_n_labels=5):
    """
    Plot genome-wide volcano plot with translocated genes highlighted
    and annotate the top N translocated genes.
    """
    df = df.copy()

    # -----------------------------
    # Load translocated genes
    # -----------------------------
    transloc_df = pd.read_csv(dist_file, sep="\t")
    transloc_genes = transloc_df[transloc_df["gene1_region_type"] == "inside_transloc"]["gene1"].unique()

    # Mark translocated genes
    df["is_translocated"] = df["gene_name"].isin(transloc_genes)

    # -----------------------------
    # Clean padj values
    # -----------------------------
    df["padj"] = df["padj"].fillna(1).replace(0, 1e-300)
    df["-log10(padj)"] = -np.log10(df["padj"])

    # -----------------------------
    # Classify differential expression
    # -----------------------------
    conditions = [
        (df["log2FoldChange"] > lfc_thresh) & (df["padj"] < padj_thresh),
        (df["log2FoldChange"] < -lfc_thresh) & (df["padj"] < padj_thresh)
    ]
    choices = ["Up", "Down"]
    df["DE_status"] = np.select(conditions, choices, default="No change")

    # -----------------------------
    # Volcano plot
    # -----------------------------
    plt.figure(figsize=(8,6))

    # Plot all genes (background)
    plt.scatter(df["log2FoldChange"], df["-log10(padj)"],
                c="lightgray", alpha=0.4, s=8)

    # Overlay translocated genes
    trans_df = df[df["is_translocated"]]
    colors = {"Up":"green", "Down":"red", "No change":"blue"}

    for status, color in colors.items():
        subset = trans_df[trans_df["DE_status"] == status]
        plt.scatter(subset["log2FoldChange"], subset["-log10(padj)"],
                    c=color, alpha=0.9, s=30, edgecolor="black", linewidth=0.3,
                    label=f"Translocated {status}")

    # Threshold lines
    plt.axvline(x=lfc_thresh, color='black', linestyle='--', lw=1)
    plt.axvline(x=-lfc_thresh, color='black', linestyle='--', lw=1)
    plt.axhline(y=-np.log10(padj_thresh), color='black', linestyle='--', lw=1)

    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(padj)")
    plt.title(f"WT vs {cond_name}\nTranslocated genes highlighted")
    plt.legend(frameon=False)
    plt.tight_layout()

    # -----------------------------
    # Annotate top translocated genes
    # -----------------------------
    trans_sig = trans_df[(abs(trans_df["log2FoldChange"]) > lfc_thresh) &
                         (trans_df["padj"] < padj_thresh)]

    # Top UP and DOWN by significance
    top_up = trans_sig[trans_sig["DE_status"] == "Up"].nlargest(top_n_labels, "-log10(padj)")
    top_down = trans_sig[trans_sig["DE_status"] == "Down"].nlargest(top_n_labels, "-log10(padj)")
    top_genes = pd.concat([top_up, top_down])

    # Print the top genes in terminal
    print(f"\nTop {top_n_labels} UP translocated genes for {cond_name}:")
    for gene in top_up["gene_name"]:
        print(gene)
    print(f"\nTop {top_n_labels} DOWN translocated genes for {cond_name}:")
    for gene in top_down["gene_name"]:
        print(gene)

    # Annotate on the plot (diagonal / offset)
    for i, (_, row) in enumerate(top_genes.iterrows()):
        x_offset = 0.05 * np.sign(row["log2FoldChange"])
        y_offset = 0.1 + 0.05*i  # stagger labels slightly upward
        plt.text(
            row["log2FoldChange"] + x_offset,
            row["-log10(padj)"] + y_offset,
            row["gene_name"],
            fontsize=6,
            rotation=30,
            ha='left' if row["log2FoldChange"] > 0 else 'right',
            va='bottom'
        )

    # -----------------------------
    # Save plot
    # -----------------------------
    outfile = os.path.join(OUTDIR, f"volcano_WT_vs_{cond_name}_highlighted.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved volcano plot: {outfile}")

    # -----------------------------
    # Export translocated gene tables
    # -----------------------------
    top_up_table = trans_sig[trans_sig["DE_status"] == "Up"].nlargest(top_n_table, "log2FoldChange")
    top_down_table = trans_sig[trans_sig["DE_status"] == "Down"].nsmallest(top_n_table, "log2FoldChange")

    top_up_table.to_csv(os.path.join(OUTDIR, f"top_up_transloc_genes_WT_vs_{cond_name}.tsv"), sep="\t", index=False)
    top_down_table.to_csv(os.path.join(OUTDIR, f"top_down_transloc_genes_WT_vs_{cond_name}.tsv"), sep="\t", index=False)

    print(f"Top {top_n_table} UP translocated genes exported")
    print(f"Top {top_n_table} DOWN translocated genes exported")

    return df


# -----------------------------
# Load DESeq2 results
# -----------------------------
de_T1 = pd.read_csv(DE_T1_FILE, sep="\t")
de_C1 = pd.read_csv(DE_C1_FILE, sep="\t")

# -----------------------------
# Generate plots
# -----------------------------
plot_volcano(de_T1, DIST_T1_FILE, "T1", top_n_table=20, top_n_labels=5)
plot_volcano(de_C1, DIST_C1_FILE, "C1", top_n_table=20, top_n_labels=5)
