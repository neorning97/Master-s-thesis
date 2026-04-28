"""
07_plot_distance_to_center.py
==============================
This script makes plots and runs statistical tests to compare
the distance-to-centre values from scripts 05 and 06.

The big question:
- Did the translocated genes MOVE in the nucleus after the translocation?
- Specifically, did they move CLOSER to the centre (active region) or
  FURTHER AWAY (inactive region)?

How we answer it:
- We have WT distances (from script 06), these are the "before" values.
- We have T1 and C1 distances (from script 05), these are the "after" values.
- We compare them with:
    1. A scatter plot (visual comparison), one dot per gene.
       The red dashed line is y = x (no change).
       Dots BELOW the line = gene moved CLOSER to the centre.
       Dots ABOVE the line = gene moved FURTHER from the centre.
    2. A Wilcoxon signed-rank test (statistical comparison).
       This tells us if the difference is statistically significant.
       We use this test (not a t-test) because distances aren't normally
       distributed.

We do this twice for each condition:
- For "inside_transloc" genes (the genes that actually got translocated)
- For "control" genes (random genes that should NOT have changed)
The control genes are a sanity check, they should show no big change.

Run scripts 05 and 06 first.
Then edit the file paths below and run this script.

Required libraries: pandas, matplotlib, seaborn, scipy
Install them with: pip install pandas matplotlib seaborn scipy
"""

import os                          
import pandas as pd                
import matplotlib.pyplot as plt    
import seaborn as sns              
from scipy.stats import wilcoxon   


# =============================================================================
# CONFIG SECTION - Edit these file paths before running the script
# =============================================================================
# I'm putting all the paths as variables at the top so they're easy to change.

# The aggregated distance files made by scripts 05 and 06
T1_FILE = "/path/to/T1_distance_to_center_agg.tsv"
C1_FILE = "/path/to/C1_distance_to_center_agg.tsv"
WT_FILE = "/path/to/WT_distance_to_center_agg.tsv"

# Where to save the plots (the folder will be created if it doesn't exist)
PLOTS_FOLDER = "/path/to/plots/distance_to_center"


# =============================================================================
# Helper function: Make plots and run stats for one condition
# =============================================================================
# This function does the work for ONE condition (either T1 or C1).
# We'll call it twice in the main code below, once for T1, once for C1.

def plot_and_test_one_condition(wt_df, cond_df, cond_col, condition_label,
                                output_folder):
    """
    For one condition, make scatter plots and run Wilcoxon tests
    separately for translocated genes and control genes.

    Parameters explained:
    - wt_df: the WT distances table (from script 06)
    - cond_df: the condition's distance table (T1 or C1, from script 05)
    - cond_col: the name of the column with the condition distances
                (either "mean_T1" or "mean_C1")
    - condition_label: a short text label like "T1" or "C1", used in
                       plot titles and filenames
    - output_folder: where to save the plots
    """
    # -------------------------------------------------------------------------
    # Step A: Merge the WT and condition tables on gene_name and region_type
    # -------------------------------------------------------------------------
    # We want each row to be one gene with BOTH its WT distance AND its
    # condition distance side by side.
    # That's what merge does: it joins two tables on shared columns.
    #
    # how="inner" means: only keep rows where the gene appears in BOTH tables.
    # If a gene is only in one of the tables, drop it.
    merged = wt_df.merge(
        cond_df,
        on=["gene_name", "region_type"],
        how="inner"
    )

    # If we got no shared genes, there's nothing to do for this condition
    if len(merged) == 0:
        print("  No shared genes for " + condition_label + " — skipping.")
        return

    # -------------------------------------------------------------------------
    # Step B: Process inside_transloc and control genes separately
    # -------------------------------------------------------------------------
    # We make a list of the gene types we want to look at,
    # then loop through each one.
    region_types_to_check = ["inside_transloc", "control"]

    for region in region_types_to_check:

        # Get only the rows where region_type matches the one we're looking at
        sub = merged[merged["region_type"] == region].copy()

        # If there are no genes of this type, skip and move on
        if len(sub) == 0:
            print("  No " + region + " genes for " + condition_label + " — skipping.")
            continue

        # Make a tag we'll use for printing and naming the output file
        # e.g. "T1_inside_transloc"
        tag = condition_label + "_" + region

        # ---------------------------------------------------------------------
        # Step C: Calculate some summary stats
        # ---------------------------------------------------------------------
        # Total number of genes in this group
        n_genes = len(sub)

        # How many of them moved CLOSER to the centre?
        # (closer = condition distance is SMALLER than WT distance)
        # The comparison gives a True/False for each row, then .sum() counts
        # the Trues (because True counts as 1 and False as 0)
        is_closer = sub[cond_col] < sub["mean_WT"]
        n_closer = is_closer.sum()

        # What fraction is that? (0.0 = none of them, 1.0 = all of them)
        fraction_closer = n_closer / n_genes

        # ---------------------------------------------------------------------
        # Step D: Run the Wilcoxon signed-rank test
        # ---------------------------------------------------------------------
        # This tests whether the WT distances and condition distances are
        # significantly different.
        # It returns two things:
        #   - stat: the test statistic (W value)
        #   - p: the p-value (smaller = more significant difference)
        # A p-value below 0.05 is usually considered "statistically significant".
        stat, p = wilcoxon(sub["mean_WT"], sub[cond_col])

        # Print the results so the user can see them in the terminal
        print("")
        print("  " + tag)
        # We use formatted strings here to control how numbers look:
        #   {:.2f}: show 2 digits after the decimal point
        #   {:.3e}: scientific notation with 3 digits (e.g. 1.234e-05)
        #   {:.1%}: show as a percentage with 1 decimal (e.g. 45.6%)
        print("    Wilcoxon signed-rank test: W={:.2f}, p={:.3e}".format(stat, p))
        print("    Genes closer to nuclear centre in " + condition_label + ": "
              + str(n_closer) + "/" + str(n_genes)
              + " ({:.1%})".format(fraction_closer))

        # ---------------------------------------------------------------------
        # Step E: Figure out the plot range
        # ---------------------------------------------------------------------
        # We need to know the smallest and biggest values across BOTH columns
        # so the y=x reference line spans the whole plot.
        wt_min = sub["mean_WT"].min()
        wt_max = sub["mean_WT"].max()
        cond_min = sub[cond_col].min()
        cond_max = sub[cond_col].max()

        # The overall min is the smaller of the two mins
        # The overall max is the bigger of the two maxes
        dist_min = min(wt_min, cond_min)
        dist_max = max(wt_max, cond_max)

        # ---------------------------------------------------------------------
        # Step F: Make the scatter plot
        # ---------------------------------------------------------------------
        # Create a new figure (the canvas) and an axes (the actual plot area)
        # figsize=(6, 6) makes it a 6x6 inch square
        fig, ax = plt.subplots(figsize=(6, 6))

        # Draw the scatter plot using seaborn
        # x = WT distance, y = condition distance, one dot per gene
        # alpha=0.6 makes dots slightly transparent so we can see overlaps
        # s=20 sets the dot size
        sns.scatterplot(
            data=sub,
            x="mean_WT",
            y=cond_col,
            alpha=0.6,
            s=20,
            ax=ax
        )

        # Draw the y = x reference line in red dashed
        # We give it just two points (the min and the max) and matplotlib
        # connects them with a straight line.
        ax.plot(
            [dist_min, dist_max],   # x values
            [dist_min, dist_max],   # y values (same as x for y=x line)
            color="red",
            linestyle="--",
            label="y = x"
        )

        # Add labels and a title
        ax.set_xlabel("WT distance to nuclear centre")
        ax.set_ylabel(condition_label + " distance to nuclear centre")

        # Build the title with stats info using a multi-line string.
        # The "\n" inside the string creates a line break in the title.
        title_text = (
            condition_label + " vs WT — " + region + "\n"
            + "Fraction closer to centre in " + condition_label + ": "
            + "{:.2f}".format(fraction_closer)
            + "  (p = {:.2e})".format(p)
        )
        ax.set_title(title_text)

        # Show the legend (so the y=x line gets labelled)
        ax.legend()

        # tight_layout adjusts spacing so labels don't get cut off
        fig.tight_layout()

        # ---------------------------------------------------------------------
        # Step G: Save the plot to a file
        # ---------------------------------------------------------------------
        # Build the full path for the output file
        # e.g. /path/to/plots/distance_to_center/T1_inside_transloc_scatter.png
        plot_filename = tag + "_scatter.png"
        plot_path = os.path.join(output_folder, plot_filename)

        # Save the figure
        fig.savefig(plot_path)

        # Close the figure to free up memory
        # (otherwise matplotlib keeps every figure open in memory)
        plt.close(fig)

        print("    Saved: " + plot_filename)


# =============================================================================
# Main code - this is where the actual work happens
# =============================================================================
# The "if __name__ == '__main__':" line is a Python convention.
# It means: only run the code below if this script is being run directly
# (not if it's being imported as a module by another script).
# For our purposes, we can just think of it as "the main code starts here".

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    # Step 1: Make sure the output folder exists
    # -------------------------------------------------------------------------
    # exist_ok=True means "don't crash if the folder already exists"
    os.makedirs(PLOTS_FOLDER, exist_ok=True)

    # -------------------------------------------------------------------------
    # Step 2: Load the three distance tables
    # -------------------------------------------------------------------------
    # The T1 and C1 files have a column called "mean_distance".
    # We rename it to "mean_T1" and "mean_C1" so we can tell them apart
    # after we merge them with the WT table.
    # The .rename() function takes a dict: {old_name: new_name}.

    print("Loading distance files...")

    t1_df = pd.read_csv(T1_FILE, sep="\t")
    t1_df = t1_df.rename(columns={"mean_distance": "mean_T1"})

    c1_df = pd.read_csv(C1_FILE, sep="\t")
    c1_df = c1_df.rename(columns={"mean_distance": "mean_C1"})

    # The WT file already calls its column "mean_WT" (script 06 named it that),
    # so we don't need to rename anything.
    wt_df = pd.read_csv(WT_FILE, sep="\t")

    print("Generating distance-to-centre plots and running statistical tests...")

    # -------------------------------------------------------------------------
    # Step 3: Process T1 vs WT
    # -------------------------------------------------------------------------
    print("")
    print("── T1 vs WT ──")
    plot_and_test_one_condition(
        wt_df=wt_df,
        cond_df=t1_df,
        cond_col="mean_T1",
        condition_label="T1",
        output_folder=PLOTS_FOLDER
    )

    # -------------------------------------------------------------------------
    # Step 4: Process C1 vs WT
    # -------------------------------------------------------------------------
    print("")
    print("── C1 vs WT ──")
    plot_and_test_one_condition(
        wt_df=wt_df,
        cond_df=c1_df,
        cond_col="mean_C1",
        condition_label="C1",
        output_folder=PLOTS_FOLDER
    )

    # -------------------------------------------------------------------------
    # Step 5: Print a summary of where the output files are
    # -------------------------------------------------------------------------
    print("")
    print("Done! All plots saved to: " + PLOTS_FOLDER)
    print("")
    print("Output plots:")

    # Loop through the conditions and gene types to print all 4 expected files
    for cond in ["T1", "C1"]:
        for region in ["inside_transloc", "control"]:
            expected_file = PLOTS_FOLDER + "/" + cond + "_" + region + "_scatter.png"
            print("  " + expected_file)
