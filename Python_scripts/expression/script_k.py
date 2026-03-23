import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import kruskal, chi2_contingency

# ------------------------
# Paths / parameters
# ------------------------
bins_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/new_plots/3d_properties/bins_annotation.csv"
gtf_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/Homo_sapiens.GRCh38.p13.protein_coding_genes.gtf"
tpm_file = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/GSE246689_gene_tpm.tsv"

de_files = {
    "T1": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/DE_WT_vs_T1.tsv",
    "C1": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/results/DE_WT_vs_C1.tsv"
}

output_folder = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/gene_bin_analysis_new"
os.makedirs(output_folder, exist_ok=True)

lfc_threshold = 1.0
padj_threshold = 0.05

# ------------------------
# Load bins
# ------------------------
bins = pd.read_csv(bins_file)
bins[["start","end"]] = bins[["start","end"]].astype(int)

# ------------------------
# Load GTF
# ------------------------
gtf = pd.read_csv(
    gtf_file,
    sep="\t",
    comment="#",
    header=None,
    names=["chrom","source","feature","start","end","score","strand","frame","attributes"]
)

gtf = gtf[gtf["feature"] == "gene"]
gtf["gene_id"] = gtf["attributes"].str.extract('gene_id "([^"]+)"')
gtf["gene_name"] = gtf["attributes"].str.extract('gene_name "([^"]+)"')

genes = gtf[["chrom","start","end","gene_id","gene_name"]].copy()
genes[["start","end"]] = genes[["start","end"]].astype(int)

# ------------------------
# Map genes to bins
# ------------------------
def assign_bin(gene_row, bins_df):
    overlap = bins_df[
        (bins_df.chrom == gene_row.chrom) &
        (bins_df.start < gene_row.end) &
        (bins_df.end > gene_row.start)
    ]

    if not overlap.empty:
        row = overlap.iloc[0]  # NOTE: first overlap (can refine later)
        return (
            row["transloc_id"],
            row["condition"],
            row["behavior"],
            row["behavior_collapsed"],
            row["change_direction"]
        )

    return [np.nan]*5


genes[[
    "transloc_id",
    "condition",
    "behavior",
    "behavior_collapsed",
    "change_direction"
]] = genes.apply(lambda r: pd.Series(assign_bin(r, bins)), axis=1)

genes_bins = genes.dropna(subset=["condition"]).copy()

# ------------------------
# Load TPM
# ------------------------
expr = pd.read_csv(tpm_file, sep="\t")

expr["WT"] = expr[["MCF10AWT_REP1","MCF10AWT_REP2","MCF10AWT_REP3"]].mean(axis=1)
expr["T1"] = expr[["MCF10AT1_REP1","MCF10AT1_REP2","MCF10AT1_REP3"]].mean(axis=1)
expr["C1"] = expr[["MCF10AC1_REP1","MCF10AC1_REP2","MCF10AC1_REP3"]].mean(axis=1)

gene_expr = genes.merge(expr[["gene_id","WT","T1","C1"]], on="gene_id", how="left")
gene_expr_bins = genes_bins.merge(expr[["gene_id","WT","T1","C1"]], on="gene_id", how="left")

# ------------------------
# Log transform + log2FC
# ------------------------
for cond in ["WT","T1","C1"]:
    gene_expr[f"log2_{cond}"] = np.log2(gene_expr[cond] + 1)
    gene_expr_bins[f"log2_{cond}"] = np.log2(gene_expr_bins[cond] + 1)

for cond in ["T1","C1"]:
    gene_expr[f"log2FC_{cond}_vs_WT"] = gene_expr[f"log2_{cond}"] - gene_expr["log2_WT"]
    gene_expr_bins[f"log2FC_{cond}_vs_WT"] = gene_expr_bins[f"log2_{cond}"] - gene_expr_bins["log2_WT"]

# ------------------------
# Load DE results
# ------------------------
for cond, file in de_files.items():
    df = pd.read_csv(file, sep="\t")

    gene_expr = gene_expr.merge(
        df[["gene_id","log2FoldChange","padj"]],
        on="gene_id",
        how="left",
        suffixes=("","_DE")
    )

    gene_expr_bins = gene_expr_bins.merge(
        df[["gene_id","log2FoldChange","padj"]],
        on="gene_id",
        how="left",
        suffixes=("","_DE")
    )

# ------------------------
# DE classification
# ------------------------
def de_call(log2fc, padj):
    if pd.notna(log2fc) and pd.notna(padj):
        if padj < padj_threshold:
            if log2fc > lfc_threshold:
                return "up"
            elif log2fc < -lfc_threshold:
                return "down"
    return "ns"


for cond in ["T1","C1"]:
    gene_expr[f"DE_{cond}"] = gene_expr.apply(
        lambda r: de_call(r[f"log2FC_{cond}_vs_WT"], r["padj"]), axis=1
    )

    gene_expr_bins[f"DE_{cond}"] = gene_expr_bins.apply(
        lambda r: de_call(r[f"log2FC_{cond}_vs_WT"], r["padj"]), axis=1
    )

# ------------------------
# Save annotated tables
# ------------------------
gene_expr.to_csv(os.path.join(output_folder, "gene_expression_genomewide.csv"), index=False)
gene_expr_bins.to_csv(os.path.join(output_folder, "gene_bin_expression_annotation.csv"), index=False)

# ========================
# ANALYSIS FUNCTION
# ========================
def run_analysis(df, cond, category_col, category_order, label):

    df = df[df["condition"] == cond].copy()

    # ---- Kruskal-Wallis ----
    groups = [
        df[df[category_col] == cat][f"log2FC_{cond}_vs_WT"].dropna()
        for cat in category_order if cat in df[category_col].values
    ]

    p_kw = kruskal(*groups).pvalue if len(groups) > 1 else np.nan

    # ---- Chi-square ----
    ct_counts = pd.crosstab(df[category_col], df[f"DE_{cond}"])
    p_chi = chi2_contingency(ct_counts)[1] if ct_counts.shape[0] > 1 else np.nan

    print(f"\n{cond} - {label}")
    print(f"Kruskal-Wallis p = {p_kw:.3e}")
    print(f"Chi-square p = {p_chi:.3e}")
    # ---- Genome-wide control ----
    df_genome = gene_expr[gene_expr[f"DE_{cond}"].notna()]
    ct_genome = pd.crosstab(
        pd.Series(["genome-wide"] * len(df_genome)),
        df_genome[f"DE_{cond}"]
    )

    # ---- Combine tables ----
    ct_combined = ct_counts.copy()
    ct_combined.loc["genome-wide"] = ct_genome.loc["genome-wide"]

    # ---- Print combined ----
    print(ct_combined)

    # ---- Genome-wide control ----
    df_genome = gene_expr[gene_expr[f"DE_{cond}"].notna()]
    ct_genome = pd.crosstab(
        pd.Series(["genome-wide"] * len(df_genome)),
        df_genome[f"DE_{cond}"]
    )

    # ---- Combine + normalize ----
    ct_plot = ct_counts.copy()
    ct_plot.loc["genome-wide"] = ct_genome.loc["genome-wide"]

    ct_frac = ct_plot.div(ct_plot.sum(axis=1), axis=0)

    # Ensure consistent DE order
    ct_frac = ct_frac.reindex(columns=["ns", "up", "down"], fill_value=0)

    # ---- Plot ----
    ct_frac.plot(
        kind="bar",
        stacked=True,
        figsize=(6, 4),
        color=["#FF0000", "#0000FF", "#cccccc"]  # ns, up, down
    )

    plt.ylabel("Fraction of genes")
    plt.title(f"{cond}: {label}")
    plt.tight_layout()

    plt.savefig(
        os.path.join(output_folder, f"{cond}_{category_col}_DE_fraction.png"),
        dpi=300
    )
    plt.close()

# ========================
# RUN ANALYSES
# ========================
for cond in ["T1","C1"]:

    # Collapsed
    run_analysis(
        gene_expr_bins,
        cond,
        "behavior_collapsed",
        ["retained","adopted","other"],
        "Collapsed compartments"
    )

    # Direction
    run_analysis(
        gene_expr_bins,
        cond,
        "change_direction",
        ["A_to_B","B_to_A","retained"],
        "A/B direction"
    )

    # Subcompartment
    run_analysis(
        gene_expr_bins,
        cond,
        "behavior",
        ["retained","adopted","other"],
        "Subcompartment behavior"
    )

print("Analysis complete.")
