"""
10_subcompartment_classification.py
=====================================
Script for classifying genomic bins in translocated regions as retained,
adopted, or other relative to the wildtype subcompartment state.

What this script does
---------------------
Chromosomal translocations move segments of DNA from one chromosomal
environment to another. One consequence is that the translocated bins may
change which chromatin subcompartment they reside in. They may retain their
original subcompartment identity, or they may adopt the subcompartment
identity of their new chromosomal neighborhood.

This script quantifies that behaviour for each translocation by:

1. Loading subcompartment annotations (A1, A2, B1, B2, B3 etc.) for WT,
   T1, and C1 from a bedGraph file.
2. For each translocated region, finding all 100 kb bins that overlap it
   and classifying each bin as:
     - **retained** : same subcompartment in the condition as in WT
     - **adopted**  : changed to match the dominant subcompartment of the
                      neighboring region on the receiving chromosome
     - **other**    : changed to something else

3. Running the same classification at the collapsed A/B level (A1/A2 -> A,
   B1/B2/B3 -> B) for a simpler view.
4. Also recording the direction of A↔B switches (A_to_B, B_to_A, retained).
5. Building a genome-wide control group (all bins NOT in any translocated
   region) to provide a background rate of subcompartment switching.
6. Running binomial tests comparing the observed fraction of each category
   in each translocation against the background rate from the control.
7. Saving stacked bar plots and a full statistics table.

Why a binomial test?
---------------------
For each translocation and each category (e.g. "adopted"), we ask: is the
fraction of bins in that category significantly higher or lower than we would
expect by chance, given the background switching rate across the rest of the
genome? The binomial test treats the number of adopted bins as a count of
successes out of n trials, where the expected probability of success is the
genome-wide background rate.

Output plots
------------
Per condition (T1 and C1), three stacked bar plots are produced:
  - *_subcompartment.png  : retained / adopted / other at full subcompartment
                            resolution (A1, A2, B1, B2, B3)
  - *_collapsed.png       : same but collapsed to A/B
  - *_direction.png       : A->B, B->A, and retained fractions

A "Genome control" bar is included in each plot for comparison.

Usage
-----
    1. Edit the CONFIG section below to point to your files.
    2. Run: python 10_subcompartment_classification.py

Dependencies
------------
    pandas, numpy, matplotlib, scipy
    Install with: pip install pandas numpy matplotlib scipy

"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binomtest


# =============================================================================
# CONFIG – edit all paths and settings here before running
# =============================================================================
CONFIG = {
    # Subcompartment bedGraph file with one column per condition
    # (WT, T1, C1). Expected columns include e.g. MCF10A_WT.state,
    # MCF10A_T1.state, MCF10A_C1.state
    "subcompartment_file": "/path/to/GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph",

    # BED files defining translocated regions for each condition
    "transloc_beds": {
        "T1": "/path/to/verify_T1_translocations.bed",
        "C1": "/path/to/verify_C1_translocations.bed",
    },

    # BED files defining neighboring regions on the receiving chromosome
    # (new neighbor) and the origin chromosome (old neighbor)
    "neighbor_beds": {
        "T1": {
            "new": "/path/to/verify_T1_neighbor.bed",
            "old": "/path/to/verify_T1_neighbor_originchr.bed",
        },
        "C1": {
            "new": "/path/to/verify_C1_neighbor.bed",
            "old": "/path/to/verify_C1_neighbor_originchr.bed",
        },
    },

    # Maps transloc_id values in the BED files to human-readable labels
    "transloc_name_maps": {
        "T1": {
            "T1": "Der(17)t(3;17)",
            "T2": "Der(3)t(3;17)",
            "T3": "t(6;19)",
        },
        "C1": {
            "T1": "t(2;10)",
            "T2": "Der(17)t(3;17)",
            "T3": "Der(3)t(3;17)",
            "T4": "t(6;19)",
        },
    },

    # Where to save output plots and tables (created automatically)
    "output_dir": "/path/to/results/subcompartment_classification",

    # Small value added to counts when computing odds ratios to avoid
    # division by zero when a category has 0 observations
    "epsilon": 1e-6,
}

CONDITIONS = ["T1", "C1"]
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def get_bins(sub: pd.DataFrame, bed_row: pd.Series) -> pd.DataFrame:
    """
    Return all subcompartment bins that overlap a given genomic region.

    A bin overlaps if any part of it falls within the region, using the
    same half-open interval logic as BED format.
    """
    return sub[
        (sub["chrom"] == bed_row["chrom"]) &
        (sub["start"]  < bed_row["end"]) &
        (sub["end"]    > bed_row["start"])
    ]


def dominant_state(df: pd.DataFrame, condition: str) -> str | None:
    """
    Return the most common subcompartment state in a set of bins for a
    given condition. Used to determine what the 'adopted' state would be
    for a translocated region — i.e. what state dominates its new
    chromosomal neighborhood.

    Returns None if the DataFrame is empty (no bins found).
    """
    col = f"MCF10A_{condition}.state"
    return df[col].value_counts().idxmax() if not df.empty else None


def dominant_state_collapsed(df: pd.DataFrame, condition: str) -> str | None:
    """Same as dominant_state but at the collapsed A/B level."""
    col = f"MCF10A_{condition}.state_collapsed"
    return df[col].value_counts().idxmax() if not df.empty else None


def classify_bin(row: pd.Series, condition: str, dominant: str | None) -> str:
    """
    Classify one bin at full subcompartment resolution as:
      - 'retained' : same state in condition as in WT
      - 'adopted'  : changed to match the dominant state of the new neighborhood
      - 'other'    : changed to something else

    When dominant is None (used for the genome control, which has no
    defined neighborhood), bins are only ever 'retained' or 'other'.
    """
    if row[f"MCF10A_{condition}.state"] == row["MCF10A_WT.state"]:
        return "retained"
    elif dominant is not None and row[f"MCF10A_{condition}.state"] == dominant:
        return "adopted"
    else:
        return "other"


def classify_bin_collapsed(
    row: pd.Series, condition: str, dominant: str | None
) -> str:
    """Same as classify_bin but at the collapsed A/B level."""
    if row[f"MCF10A_{condition}.state_collapsed"] == row["MCF10A_WT.state_collapsed"]:
        return "retained"
    elif dominant is not None and row[f"MCF10A_{condition}.state_collapsed"] == dominant:
        return "adopted"
    else:
        return "other"


def get_change_direction(row: pd.Series, condition: str) -> str:
    """
    Classify the direction of compartment switching at the A/B level:
      - 'A_to_B' : was A in WT, now B
      - 'B_to_A' : was B in WT, now A
      - 'retained' : no A/B switch
    """
    wt   = row["MCF10A_WT.state_collapsed"]
    cond = row[f"MCF10A_{condition}.state_collapsed"]
    if   wt == "A" and cond == "B": return "A_to_B"
    elif wt == "B" and cond == "A": return "B_to_A"
    else:                            return "retained"


def get_background_rates(df: pd.DataFrame, col: str) -> dict[str, float]:
    """
    Compute the fraction of bins in each category across a DataFrame.
    Used to get the genome-wide background rate for the binomial test.
    """
    vc = df[col].value_counts()
    n  = len(df)
    return {k: vc.get(k, 0) / n for k in vc.index}


def plot_stacked_bar(
    df: pd.DataFrame,
    y_cols: list[str],
    filename: str,
    title: str,
    output_dir: str,
    colors: list[str] | None = None,
    use_colormap: bool = False,
) -> None:
    """
    Save a stacked bar plot where each bar is one translocation and the
    height of each segment shows the fraction of bins in that category.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    kwargs = dict(x="transloc_id", kind="bar", stacked=True, y=y_cols, ax=ax)
    if use_colormap:
        df.plot(**kwargs, colormap="tab10")
    else:
        df.plot(**kwargs, color=colors)
    ax.set_ylabel("Fraction of bins")
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close()


# =============================================================================
# Main loop over conditions
# =============================================================================

os.makedirs(CONFIG["output_dir"], exist_ok=True)

# ── Load subcompartment annotations ──────────────────────────────────────────
sub = pd.read_csv(CONFIG["subcompartment_file"], sep="\t")
sub[["start", "end"]] = sub[["start", "end"]].astype(int)

# Collapse full subcompartment labels (A1, A2, B1, B2, B3) to A/B for the
# simpler binary analysis
for cond in ["WT"] + CONDITIONS:
    sub[f"MCF10A_{cond}.state_collapsed"] = sub[f"MCF10A_{cond}.state"].str[0]

all_bins      = []
stats_results = []

for cond in CONDITIONS:
    print(f"\n===== Processing {cond} =====")

    trans  = pd.read_csv(CONFIG["transloc_beds"][cond], sep="\t")
    new_nb = pd.read_csv(CONFIG["neighbor_beds"][cond]["new"], sep="\t")
    old_nb = pd.read_csv(CONFIG["neighbor_beds"][cond]["old"], sep="\t")

    for df_ in (trans, new_nb, old_nb):
        df_[["start", "end"]] = df_[["start", "end"]].astype(int)

    trans_name_map = CONFIG["transloc_name_maps"][cond]

    results_sub       = []
    results_collapsed = []
    results_direction = []
    all_trans_bins    = []

    # Store raw counts for the binomial test (computed after fractions)
    trans_counts_sub       = {}
    trans_counts_collapsed = {}
    trans_counts_dir       = {}

    # ── Per-translocation classification ─────────────────────────────────────
    for tid in trans["transloc_id"].unique():
        trans_rows = trans[trans["transloc_id"] == tid]
        new_rows   = new_nb[new_nb["transloc_id"] == tid]

        # Collect all subcompartment bins overlapping this translocation
        trans_bins = pd.concat(
            [get_bins(sub, r) for _, r in trans_rows.iterrows()],
            ignore_index=True,
        )
        # Collect bins in the new neighborhood to determine the adopted state
        new_bins = pd.concat(
            [get_bins(sub, r) for _, r in new_rows.iterrows()],
            ignore_index=True,
        )

        for df_ in (trans_bins, new_bins):
            df_.drop_duplicates(subset=["chrom", "start", "end"], inplace=True)

        if trans_bins.empty:
            print(f"  {tid}: no overlapping bins, skipping.")
            continue

        # The 'adopted' state is the dominant state in the new neighborhood
        new_dom           = dominant_state(new_bins, cond)
        new_dom_collapsed = dominant_state_collapsed(new_bins, cond)

        # Classify each bin
        trans_bins["behavior"]          = trans_bins.apply(
            lambda r: classify_bin(r, cond, new_dom), axis=1
        )
        trans_bins["behavior_collapsed"] = trans_bins.apply(
            lambda r: classify_bin_collapsed(r, cond, new_dom_collapsed), axis=1
        )
        trans_bins["change_direction"]  = trans_bins.apply(
            lambda r: get_change_direction(r, cond), axis=1
        )

        trans_bins["transloc_id"] = tid
        trans_bins["condition"]   = cond

        all_bins.append(trans_bins)
        all_trans_bins.append(trans_bins[["chrom", "start", "end"]])

        n = len(trans_bins)

        # Store raw counts for statistical testing
        trans_counts_sub[tid]       = trans_bins["behavior"].value_counts()
        trans_counts_collapsed[tid] = trans_bins["behavior_collapsed"].value_counts()
        trans_counts_dir[tid]       = trans_bins["change_direction"].value_counts()

        # Store normalised fractions for plotting
        results_sub.append({
            "transloc_id": tid,
            **(trans_counts_sub[tid] / n),
        })
        results_collapsed.append({
            "transloc_id": tid,
            **(trans_counts_collapsed[tid] / n),
        })
        results_direction.append({
            "transloc_id": tid,
            **(trans_counts_dir[tid] / n),
        })

    # ── Genome-wide control ───────────────────────────────────────────────────
    # All bins NOT in any translocated region serve as the background.
    # Their switching rate gives the expected probability for the binomial test.
    all_trans_bins_df = pd.concat(all_trans_bins).drop_duplicates()
    control_bins = sub.merge(
        all_trans_bins_df, on=["chrom", "start", "end"],
        how="left", indicator=True,
    )
    control_bins = control_bins[control_bins["_merge"] == "left_only"].drop(
        columns="_merge"
    )

    # Classify control bins (dominant=None because they have no defined
    # neighborhood — they can only be retained or other)
    control_bins["behavior"]          = control_bins.apply(
        lambda r: classify_bin(r, cond, None), axis=1
    )
    control_bins["behavior_collapsed"] = control_bins.apply(
        lambda r: classify_bin_collapsed(r, cond, None), axis=1
    )
    control_bins["change_direction"]  = control_bins.apply(
        lambda r: get_change_direction(r, cond), axis=1
    )

    # Background rates used as the null hypothesis in the binomial test
    control_p_sub       = get_background_rates(control_bins, "behavior")
    control_p_collapsed = get_background_rates(control_bins, "behavior_collapsed")
    control_p_dir       = get_background_rates(control_bins, "change_direction")

    # ── Binomial tests ────────────────────────────────────────────────────────
    # For each translocation and each category, test whether the observed
    # fraction differs significantly from the genome-wide background rate.
    eps = CONFIG["epsilon"]

    for tid in trans_counts_sub:
        n = int(trans_counts_sub[tid].sum())

        for level, counts_dict, control_p in [
            ("subcompartment", trans_counts_sub[tid],       control_p_sub),
            ("collapsed",      trans_counts_collapsed[tid], control_p_collapsed),
            ("direction",      trans_counts_dir[tid],       control_p_dir),
        ]:
            for cat in counts_dict.index:
                k  = int(counts_dict.get(cat, 0))
                p0 = control_p.get(cat, 0)

                test = binomtest(k, n, p=p0)

                # Log2 odds ratio: positive = enriched vs background,
                #                  negative = depleted vs background
                odds_trans   = (k + eps) / (n - k + eps)
                odds_control = (p0 + eps) / (1 - p0 + eps)

                stats_results.append({
                    "condition":         cond,
                    "translocation":     trans_name_map.get(tid, tid),
                    "level":             level,
                    "category":          cat,
                    "n_bins":            n,
                    "k_bins":            k,
                    "fraction":          k / n,
                    "expected_fraction": p0,
                    "p_value":           test.pvalue,
                    "log2_odds_ratio":   np.log2(odds_trans / odds_control),
                })

    # ── Add control bar to plots ──────────────────────────────────────────────
    def add_control(results: list, col: str) -> list:
        """Append a 'Genome control' entry showing the background rates."""
        vc = control_bins[col].value_counts(normalize=True)
        return results + [{"transloc_id": "Genome control", **vc.to_dict()}]

    results_sub       = add_control(results_sub,       "behavior")
    results_collapsed = add_control(results_collapsed, "behavior_collapsed")
    results_direction = add_control(results_direction, "change_direction")

    # ── Build plot DataFrames ─────────────────────────────────────────────────
    def make_plot_df(results: list) -> pd.DataFrame:
        """Convert results list to DataFrame and apply human-readable labels."""
        df = pd.DataFrame(results).fillna(0)
        df["transloc_id"] = df["transloc_id"].map(
            lambda x: trans_name_map.get(x, x)
        )
        return df

    final_sub       = make_plot_df(results_sub)
    final_collapsed = make_plot_df(results_collapsed)
    final_direction = make_plot_df(results_direction)

    # ── Save plots ────────────────────────────────────────────────────────────
    plot_stacked_bar(
        final_sub,
        ["retained", "adopted", "other"],
        f"{cond}_subcompartment.png",
        f"{cond} vs WT: Subcompartments",
        CONFIG["output_dir"],
        use_colormap=True,
    )

    plot_stacked_bar(
        final_collapsed,
        ["retained", "adopted", "other"],
        f"{cond}_collapsed.png",
        f"{cond} vs WT: Collapsed A/B",
        CONFIG["output_dir"],
        use_colormap=True,
    )

    plot_stacked_bar(
        final_direction,
        ["A_to_B", "B_to_A", "retained"],
        f"{cond}_direction.png",
        f"{cond} vs WT: A↔B direction",
        CONFIG["output_dir"],
        colors=["#FF0800", "#1A43BF", "#cccccc"],
    )

    print(f"  Plots saved for {cond}")

# =============================================================================
# Save combined outputs
# =============================================================================

bins_out  = os.path.join(CONFIG["output_dir"], "bins_annotation.csv")
stats_out = os.path.join(CONFIG["output_dir"], "all_binomial_stats.csv")

pd.concat(all_bins).to_csv(bins_out, index=False)
pd.DataFrame(stats_results).to_csv(stats_out, index=False)

print(f"\nDone.")
print(f"  Bin annotations: {bins_out}")
print(f"  Statistics:      {stats_out}")
