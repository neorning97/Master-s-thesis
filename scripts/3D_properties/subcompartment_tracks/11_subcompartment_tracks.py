"""
11_subcompartment_tracks.py
============================
Script for visualising chromatin subcompartment states along each derivative
chromosome formed by chromosomal translocations in the T1 and C1 cell lines.

What this script does
---------------------
Each translocation fuses portions of two chromosomes into a new derivative
chromosome. To understand how the translocated segment's chromatin environment
changes, this script draws horizontal subcompartment tracks showing:

  1. The full original chromosome A (the one contributing the fragment)
  2. The full original chromosome B (the one contributing the body)
  3. The assembled derivative chromosome: fragment | body, joined at the
     breakpoint with a dashed vertical line

One figure is produced per translocation per condition. Translocations are
read directly from the BED files.

Translocation structure
-----------------------
Each translocation in the BED file consists of two paired rows sharing a
transloc_id. One row is the "fragment" (the smaller segment that moves onto a
new chromosome) and one is the "body" (the larger receiving segment). The
script identifies which is which by comparing row lengths: the shorter row is
treated as the fragment.

The translocations covered are:

  T1 only:
    - t(3;17)  — chr3:0–58.6Mb fragment -> chr17:22.7–83.3Mb body (reciprocal:
                  chr3:58.6–198.3Mb body <- chr17:0–22.7Mb fragment)
    - t(6;19)  — chr19:0–31.9Mb fragment -> chr6:0–170.8Mb body

  T1 and C1:
    - t(3;17) and t(6;19) as above

  C1 only:
    - t(2;10)  — chr2:0–131.7Mb fragment -> chr10:39.8–133.8Mb body

Subcompartment colour scheme
-----------------------------
A-compartment subcompartments are shown in shades of red (A0 = light,
A3 = dark). B-compartment subcompartments are shown in shades of blue
(B0 = light, B3 = dark). Grey is used for unrecognised states.

Usage
-----
    1. Edit the CONFIG section below to point to your files.
    2. Run: python 11_subcompartment_tracks.py

Dependencies
------------
    pandas, matplotlib
    Install with: pip install pandas matplotlib

"""

import os

import pandas as pd
import matplotlib.pyplot as plt


# =============================================================================
# CONFIG – edit all paths and settings here before running
# =============================================================================
CONFIG = {
    # Subcompartment bedGraph file (100 kb resolution, WT + T1 + C1)
    "subcompartment_file": "/path/to/GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph",

    # Translocation BED files — each row is one segment of a translocation.
    # Paired rows share a transloc_id. The script reads coordinates directly
    # from these files so no breakpoints need to be hardcoded.
    "transloc_beds": {
        "T1": "/path/to/verify_T1_translocations.bed",
        "C1": "/path/to/verify_C1_translocations.bed",
    },

    # Which conditions to plot for each translocation.
    # t(2;10) only occurs in C1, so WT and C1 are plotted for that one.
    # t(3;17) and t(6;19) occur in both T1 and C1, so WT, T1, and C1 are plotted.
    # Conditions listed here must match columns in the subcompartment file.
    "conditions": {
        "WT": "MCF10A_WT.state",
        "T1": "MCF10A_T1.state",
        "C1": "MCF10A_C1.state",
    },

    # Human-readable names for each transloc_id, per condition BED file.
    # Used in plot titles and output filenames.
    "transloc_name_maps": {
        "T1": {
            "T1": "Der(17)t(3;17)",   # chr3 fragment -> chr17 body
            "T2": "Der(3)t(3;17)",    # chr17 fragment -> chr3 body (reciprocal)
            "T3": "t(6;19)",         # chr19 fragment -> chr6 body
        },
        "C1": {
            "T1": "t(2;10)",         # chr2 fragment -> chr10 body
            "T2": "Der(17)t(3;17)",
            "T3": "Der(3)t(3;17)",
            "T4": "t(6;19)",
        },
    },

    # Which conditions to include for each translocation name.
    # Translocations that only appear in one cell line should not include
    # the condition where they don't exist.
    "transloc_conditions": {
        "t(2;10)":       ["WT", "C1"],
        "Der(17)t(3;17)": ["WT", "T1", "C1"],
        "Der(3)t(3;17)":  ["WT", "T1", "C1"],
        "t(6;19)":       ["WT", "T1", "C1"],
    },

    # Where to save output plots (created automatically)
    "output_dir": "/path/to/results/subcompartment_tracks",

    # Colour scheme: A-compartment = shades of red, B = shades of blue
    "subcomp_colors": {
        "A0": "#fcaeae", "A1": "#fb6a6a", "A2": "#ef3b2c", "A3": "#99000d",
        "B0": "#deebf7", "B1": "#9ecae1", "B2": "#4292c6", "B3": "#084594",
    },
}
# =============================================================================


# =============================================================================
# Helper functions
# =============================================================================

def load_transloc_pairs(bed_path: str, name_map: dict) -> list[dict]:
    """
    Load a translocation BED file and pair rows that share a transloc_id.

    Each pair consists of a "fragment" row (the shorter segment that moves)
    and a "body" row (the longer receiving segment). The fragment is identified
    as the row with the smaller genomic span within each pair.

    Returns a list of dicts, each with keys:
        name       : human-readable translocation label
        transloc_id: original transloc_id from the BED file
        fragment   : Series with chrom, start, end for the fragment
        body       : Series with chrom, start, end for the body
    """
    bed = pd.read_csv(bed_path, sep="\t")
    bed[["start", "end"]] = bed[["start", "end"]].astype(int)
    bed["span"] = bed["end"] - bed["start"]

    pairs = []
    for tid in bed["transloc_id"].unique():
        rows = bed[bed["transloc_id"] == tid]
        if len(rows) != 2:
            print(f"  Warning: transloc_id {tid} has {len(rows)} rows, expected 2 — skipping.")
            continue

        # The shorter row is the fragment; the longer is the body
        rows_sorted = rows.sort_values("span")
        fragment    = rows_sorted.iloc[0]
        body        = rows_sorted.iloc[1]

        name = name_map.get(tid, tid)
        pairs.append({
            "name":        name,
            "transloc_id": tid,
            "fragment":    fragment,
            "body":        body,
        })

    return pairs


def plot_track(
    ax,
    segments: list[tuple[pd.DataFrame, int]],
    condition_col: str,
    title: str,
    subcomp_colors: dict,
    x_label_suffix: str = "",
) -> None:
    """
    Draw a horizontal subcompartment track on the given axes.

    Each 100 kb bin is drawn as a filled rectangle coloured by its
    subcompartment state. Multiple genomic segments can be concatenated
    into a single track (used for the derivative chromosome, which joins
    a fragment from one chromosome with a body from another).

    Parameters
    ----------
    ax              : Matplotlib axes to draw on.
    segments        : List of (DataFrame, original_start) tuples. Each tuple
                      is one contiguous genomic segment. original_start is the
                      genomic coordinate where that segment begins, used to
                      correctly position bins along the plot axis.
    condition_col   : Column name for subcompartment states in this condition.
    title           : Plot title.
    subcomp_colors  : Dict mapping subcompartment labels to hex colours.
    x_label_suffix  : Optional suffix for the x-axis label, e.g. "[chr3]".
    """
    cursor = 0  # current x-position along the plot axis

    for df_seg, original_start in segments:
        df_seg = df_seg.copy()

        # Shift genomic coordinates so this segment starts at 'cursor'
        shift = cursor - original_start
        df_seg["plot_start"] = df_seg["start"] + shift
        df_seg["plot_end"]   = df_seg["end"]   + shift

        for _, row in df_seg.iterrows():
            color = subcomp_colors.get(row[condition_col], "#cccccc")
            ax.fill_between(
                [row["plot_start"], row["plot_end"]],
                0, 1,
                color=color,
                edgecolor="none",
            )

        seg_length = df_seg["end"].max() - original_start
        cursor += seg_length

    ax.set_xlim(0, cursor)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(f"Position (bp){x_label_suffix}")

    # Draw breakpoint line when two segments are joined
    if len(segments) == 2:
        seg1_len = segments[0][0]["end"].max() - segments[0][1]
        ax.axvline(
            x=seg1_len, color="green", linestyle="-",
            linewidth=5, label="Breakpoint",
        )
        ax.legend(fontsize=8)


def add_legend(ax, subcomp_colors: dict) -> None:
    """Draw a subcompartment colour legend and hide the axes frame."""
    for state, color in subcomp_colors.items():
        ax.plot([], [], color=color, label=state, linewidth=6)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    ax.axis("off")


def get_chrom_bins(
    sub: pd.DataFrame,
    chrom: str,
    start: int | None = None,
    end: int | None = None,
) -> pd.DataFrame:
    """
    Return all subcompartment bins on a given chromosome, optionally
    filtered to a genomic range [start, end).
    """
    mask = sub["chrom"] == chrom
    if start is not None:
        mask &= sub["start"] >= start
    if end is not None:
        mask &= sub["end"] <= end
    return sub[mask]


# =============================================================================
# Main
# =============================================================================

os.makedirs(CONFIG["output_dir"], exist_ok=True)

# ── Load subcompartment data ──────────────────────────────────────────────────
sub = pd.read_csv(CONFIG["subcompartment_file"], sep="\t")
print(f"Loaded subcompartment data: {len(sub)} bins")

# ── Load and pair translocation rows from BED files ──────────────────────────
# Collect all unique translocation pairs across T1 and C1 BED files.
# We index by name so each translocation is only plotted once even if it
# appears in both T1 and C1 BED files (e.g. t(3;17) and t(6;19)).
all_pairs: dict[str, dict] = {}

for cond_key, bed_path in CONFIG["transloc_beds"].items():
    name_map = CONFIG["transloc_name_maps"][cond_key]
    pairs    = load_transloc_pairs(bed_path, name_map)
    for p in pairs:
        if p["name"] not in all_pairs:
            all_pairs[p["name"]] = p
            print(f"  Loaded: {p['name']}  "
                  f"fragment={p['fragment']['chrom']}:{p['fragment']['start']:,}–{p['fragment']['end']:,}  "
                  f"body={p['body']['chrom']}:{p['body']['start']:,}–{p['body']['end']:,}")

# ── Generate one figure per translocation per condition ───────────────────────
for trans_name, pair in all_pairs.items():

    frag = pair["fragment"]
    body = pair["body"]

    frag_chrom = frag["chrom"]
    body_chrom = body["chrom"]

    # Full original chromosomes (for context tracks)
    frag_chrom_all = get_chrom_bins(sub, frag_chrom)
    body_chrom_all = get_chrom_bins(sub, body_chrom)

    # Sliced segments that form the derivative chromosome
    frag_bins = get_chrom_bins(sub, frag_chrom, frag["start"], frag["end"])
    body_bins = get_chrom_bins(sub, body_chrom, body["start"], body["end"])

    # Determine which conditions to plot for this translocation
    conds_to_plot = CONFIG["transloc_conditions"].get(
        trans_name, list(CONFIG["conditions"].keys())
    )

    for cond_name in conds_to_plot:
        col = CONFIG["conditions"].get(cond_name)
        if col is None:
            print(f"  Skipping {cond_name}: not found in conditions config.")
            continue

        print(f"\nPlotting {trans_name} — {cond_name}")

        fig, axes = plt.subplots(4, 1, figsize=(14, 9))
        fig.suptitle(
            f"Subcompartment tracks {trans_name} — {cond_name}",
            fontsize=13, fontweight="bold",
        )

        # Track 1: full original fragment chromosome
        plot_track(
            axes[0],
            [(frag_chrom_all, 0)],
            col,
            title=f"{frag_chrom} (original) — {cond_name}",
            subcomp_colors=CONFIG["subcomp_colors"],
            x_label_suffix=f" [{frag_chrom}]",
        )

        # Track 2: full original body chromosome
        plot_track(
            axes[1],
            [(body_chrom_all, 0)],
            col,
            title=f"{body_chrom} (original) — {cond_name}",
            subcomp_colors=CONFIG["subcomp_colors"],
            x_label_suffix=f" [{body_chrom}]",
        )

        # Track 3: assembled derivative chromosome
        # Fragment segment joined to body segment, breakpoint marked
        frag_mb_end  = frag["end"] / 1_000_000
        body_mb_start = body["start"] / 1_000_000
        body_mb_end   = body["end"] / 1_000_000

        plot_track(
            axes[2],
            [(frag_bins, frag["start"]), (body_bins, body["start"])],
            col,
            title=(
                f"Derivative {trans_name} — {cond_name} | "
                f"[{frag_chrom}:{frag['start']//1_000_000}–{frag_mb_end:.1f}Mb | "
                f"{body_chrom}:{body_mb_start:.1f}–{body_mb_end:.1f}Mb]"
            ),
            subcomp_colors=CONFIG["subcomp_colors"],
        )

        # Track 4: legend only
        add_legend(axes[3], CONFIG["subcomp_colors"])

        plt.tight_layout()

        safe_name = trans_name.replace(";", "_").replace("(", "").replace(")", "")
        out_path  = os.path.join(
            CONFIG["output_dir"],
            f"subcompartment_{safe_name}_{cond_name}.png",
        )
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out_path}")

print("\nDone.")
