"""
07_plot_distance_to_center.py
==============================
Script 3 of 3 for nuclear distance analysis – Visualisation and statistical
testing of distance-to-centre results.
What this script does
---------------------
Loads the aggregated distance-to-centre tables from scripts 05 and 06, then
for each condition (T1, C1) and each gene type (inside_transloc and control)
produces:

  1. A scatter plot of WT distance vs condition distance.
     Each point is one gene. The red dashed line is y = x (equal distance).
     Points below the line = gene moved CLOSER to the nuclear centre
     after the translocation.
     Points above the line = gene moved FURTHER from the nuclear centre.

  2. A Wilcoxon signed-rank test comparing WT and condition distances.
     This non-parametric test is used because distance distributions are
     not normally distributed. It tests whether the median distance to
     the nuclear centre is significantly shifted between WT and the
     translocated condition.

The key biological question
----------------------------
Do translocated genes (inside_transloc) move toward the nuclear centre after
the translocation? In Chrom3D models, the nuclear centre corresponds to the
transcriptionally active interior. Inward movement after translocation would
support the hypothesis that compartment adoption (moving into an active
environment) drives the observed gene expression changes.

The control group (randomly sampled genes from non-translocated chromosomes)
should show no systematic shift, confirming that any effect in the transloc
group is specific to the translocation.

Usage
-----
    1. Run scripts 05 and 06 first.
    2. Edit the CONFIG section below to point to your files.
    3. Run: python 07_plot_distance_to_center.py

Dependencies
------------
    pandas, matplotlib, seaborn, scipy
    Install with: pip install pandas matplotlib seaborn scipy

Author: Nadia Ørning
Date:   08/01/26
"""

import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon


# =============================================================================
# CONFIG – edit all paths here before running
# =============================================================================
CONFIG = {
    # Aggregated distance files from scripts 05 and 06
    "t1_agg_file": "/path/to/T1_distance_to_center_agg.tsv",
    "c1_agg_file": "/path/to/C1_distance_to_center_agg.tsv",
    "wt_agg_file": "/path/to/WT_distance_to_center_agg.tsv",

    # Where to save the output plots (created automatically)
    "plots_dir": "/path/to/plots/distance_to_center",
}
# =============================================================================


# =============================================================================
# Plotting function
# =============================================================================

def plot_wt_vs_condition(
    wt_df: pd.DataFrame,
    cond_df: pd.DataFrame,
    cond_col: str,
    label: str,
    output_dir: str,
) -> None:
    """
    For one condition (T1 or C1), generate scatter plots and run Wilcoxon
    tests separately for translocated genes and control genes.

    Parameters
    ----------
    wt_df      : WT aggregated distances (from script 06).
    cond_df    : Condition (T1 or C1) aggregated distances (from script 05).
    cond_col   : Column name for the condition distances ('mean_T1' or 'mean_C1').
    label      : Condition label used in plot titles and filenames ('T1' or 'C1').
    output_dir : Directory where plots are saved.
    """
    # Merge WT and condition on gene name and region type so each point
    # represents the same gene measured in both WT and the condition
    merged = wt_df.merge(cond_df, on=["gene_name", "region_type"], how="inner")

    if merged.empty:
        print(f"  No shared genes for {label} — skipping.")
        return

    # Analyse inside_transloc and control genes separately
    for region in ["inside_transloc", "control"]:
        sub = merged[merged["region_type"] == region].copy()
        if sub.empty:
            print(f"  No {region} genes for {label} — skipping.")
            continue

        tag = f"{label}_{region}"

        # ── Summary statistics ────────────────────────────────────────────────
        n_genes      = len(sub)
        n_closer     = (sub[cond_col] < sub["mean_WT"]).sum()
        frac_closer  = n_closer / n_genes
        stat, p      = wilcoxon(sub["mean_WT"], sub[cond_col])

        print(f"\n  {tag}")
        print(f"    Wilcoxon signed-rank test: W={stat:.2f}, p={p:.3e}")
        print(f"    Genes closer to nuclear centre in {label}: "
              f"{n_closer}/{n_genes} ({frac_closer:.1%})")

        dist_min = min(sub["mean_WT"].min(), sub[cond_col].min())
        dist_max = max(sub["mean_WT"].max(), sub[cond_col].max())

        # ── Scatter plot ──────────────────────────────────────────────────────
        # Each point is one gene. Points below y=x moved closer to the centre.
        fig, ax = plt.subplots(figsize=(6, 6))
        sns.scatterplot(
            data=sub,
            x="mean_WT", y=cond_col,
            alpha=0.6, s=20, ax=ax,
        )
        ax.plot([dist_min, dist_max], [dist_min, dist_max],
                color="red", linestyle="--", label="y = x")
        ax.set_xlabel("WT distance to nuclear centre")
        ax.set_ylabel(f"{label} distance to nuclear centre")
        ax.set_title(
            f"{label} vs WT — {region}\n"
            f"Fraction closer to centre in {label}: {frac_closer:.2f}  "
            f"(p = {p:.2e})"
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, f"{tag}_scatter.png"))
        plt.close(fig)

        print(f"    Saved: {tag}_scatter.png")


# =============================================================================
# Run
# =============================================================================

if __name__ == "__main__":
    os.makedirs(CONFIG["plots_dir"], exist_ok=True)

    # Load and prepare the three distance tables
    t1_df = pd.read_csv(CONFIG["t1_agg_file"], sep="\t").rename(
        columns={"mean_distance": "mean_T1"}
    )
    c1_df = pd.read_csv(CONFIG["c1_agg_file"], sep="\t").rename(
        columns={"mean_distance": "mean_C1"}
    )
    wt_df = pd.read_csv(CONFIG["wt_agg_file"], sep="\t")
    # Script 06 already writes the WT column as 'mean_WT'

    print("Generating distance-to-centre plots and running statistical tests...")

    # ── T1 vs WT ──────────────────────────────────────────────────────────────
    print("\n── T1 vs WT ──")
    plot_wt_vs_condition(wt_df, t1_df, "mean_T1", "T1", CONFIG["plots_dir"])

    # ── C1 vs WT ──────────────────────────────────────────────────────────────
    print("\n── C1 vs WT ──")
    plot_wt_vs_condition(wt_df, c1_df, "mean_C1", "C1", CONFIG["plots_dir"])

    print(f"\nDone. All plots saved to: {CONFIG['plots_dir']}")
    print("\nOutput plots:")
    for cond in ["T1", "C1"]:
        for region in ["inside_transloc", "control"]:
            print(f"  {CONFIG['plots_dir']}/{cond}_{region}_scatter.png")
