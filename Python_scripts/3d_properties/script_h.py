import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import binomtest
import numpy as np

# ------------------------
# Parameters
# ------------------------
conditions = ["T1", "C1"]
output_folder = "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/new_plots/3d_properties_fixed"
os.makedirs(output_folder, exist_ok=True)

epsilon = 1e-6  # for odds ratio stability

# ------------------------
# Load data
# ------------------------
sub = pd.read_csv(
    "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/subcompartments/GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph",
    sep="\t"
)

trans_files = {
    "T1": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/verify_T1_translocations.bed",
    "C1": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/verify_C1_translocations.bed"
}

neighbor_files = {
    "T1": {
        "new": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/verify_T1_neighbor1.bed",
        "old": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/verify_T1_neighbor_originchr1.bed"
    },
    "C1": {
        "new": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/verify_C1_neighbor1.bed",
        "old": "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/verify_C1_neighbor_originchr1.bed"
    }
}

sub[["start","end"]] = sub[["start","end"]].astype(int)

# Collapse to A/B
for cond in ["WT"] + conditions:
    sub[f"MCF10A_{cond}.state_collapsed"] = sub[f"MCF10A_{cond}.state"].str[0]

# ------------------------
# Mapping
# ------------------------
t1_mapping = {"T1":"Der(3)t(3;17)","T2":"Der(17)t(3;17)","T3":"t(6;19)"}
c1_mapping = {"T1":"t(2;10)","T2":"Der(3)t(3;17)","T3":"Der(17)t(3;17)","T4":"t(6;19)"}

# ------------------------
# Functions
# ------------------------
def get_bins(sub, bed_row):
    return sub[(sub.chrom == bed_row.chrom) &
               (sub.start < bed_row.end) &
               (sub.end > bed_row.start)]

def dominant_state(df, condition):
    return df[f"MCF10A_{condition}.state"].value_counts().idxmax() if not df.empty else None

def dominant_state_collapsed(df, condition):
    return df[f"MCF10A_{condition}.state_collapsed"].value_counts().idxmax() if not df.empty else None

def classify_bin(row, condition, dominant):
    if row[f"MCF10A_{condition}.state"] == row["MCF10A_WT.state"]:
        return "retained"
    elif dominant is not None and row[f"MCF10A_{condition}.state"] == dominant:
        return "adopted"
    else:
        return "other"

def classify_bin_collapsed(row, condition, dominant):
    if row[f"MCF10A_{condition}.state_collapsed"] == row["MCF10A_WT.state_collapsed"]:
        return "retained"
    elif dominant is not None and row[f"MCF10A_{condition}.state_collapsed"] == dominant:
        return "adopted"
    else:
        return "other"

def get_change_direction(row, condition):
    wt = row["MCF10A_WT.state_collapsed"]
    cond_val = row[f"MCF10A_{condition}.state_collapsed"]
    if wt == "A" and cond_val == "B": return "A_to_B"
    elif wt == "B" and cond_val == "A": return "B_to_A"
    else: return "retained"

# ------------------------
# Main
# ------------------------
all_bins = []
stats_results = []

for cond in conditions:

    print(f"\nProcessing {cond}")

    trans = pd.read_csv(trans_files[cond], sep="\t")
    new_nb = pd.read_csv(neighbor_files[cond]["new"], sep="\t")
    old_nb = pd.read_csv(neighbor_files[cond]["old"], sep="\t")

    for df_ in [trans,new_nb,old_nb]:
        df_[["start","end"]] = df_[["start","end"]].astype(int)

    trans_name_map = t1_mapping if cond=="T1" else c1_mapping

    results_sub, results_collapsed, results_direction = [], [], []
    all_trans_bins = []

    trans_counts_sub = {}
    trans_counts_collapsed = {}
    trans_counts_dir = {}

    # ------------------------
    # Per translocation
    # ------------------------
    for tid in trans.transloc_id.unique():

        trans_row = trans[trans.transloc_id == tid]
        new_rows = new_nb[new_nb.transloc_id == tid]

        trans_bins = pd.concat([get_bins(sub,r) for _,r in trans_row.iterrows()], ignore_index=True)
        new_bins = pd.concat([get_bins(sub,r) for _,r in new_rows.iterrows()], ignore_index=True)

        for df_ in [trans_bins,new_bins]:
            df_.drop_duplicates(subset=["chrom","start","end"], inplace=True)

        if trans_bins.empty:
            continue

        new_dom = dominant_state(new_bins, cond)
        new_dom_collapsed = dominant_state_collapsed(new_bins, cond)

        # classify
        trans_bins["behavior"] = trans_bins.apply(lambda r: classify_bin(r,cond,new_dom), axis=1)
        trans_bins["behavior_collapsed"] = trans_bins.apply(lambda r: classify_bin_collapsed(r,cond,new_dom_collapsed), axis=1)
        trans_bins["change_direction"] = trans_bins.apply(lambda r: get_change_direction(r,cond), axis=1)

        trans_bins["transloc_id"] = tid
        trans_bins["condition"] = cond

        all_bins.append(trans_bins)
        all_trans_bins.append(trans_bins[["chrom","start","end"]])

        n = len(trans_bins)

        # counts
        trans_counts_sub[tid] = trans_bins["behavior"].value_counts()
        trans_counts_collapsed[tid] = trans_bins["behavior_collapsed"].value_counts()
        trans_counts_dir[tid] = trans_bins["change_direction"].value_counts()

        # normalized for plotting
        results_sub.append({
            "transloc_id": tid,
            **(trans_counts_sub[tid] / n)
        })

        results_collapsed.append({
            "transloc_id": tid,
            **(trans_counts_collapsed[tid] / n)
        })

        results_direction.append({
            "transloc_id": tid,
            **(trans_counts_dir[tid] / n)
        })

    # ------------------------
    # CONTROL
    # ------------------------
    all_trans_bins_df = pd.concat(all_trans_bins).drop_duplicates()
    control_bins = sub.merge(all_trans_bins_df, on=["chrom","start","end"], how="left", indicator=True)
    control_bins = control_bins[control_bins["_merge"]=="left_only"].drop(columns="_merge")

    control_bins["behavior"] = control_bins.apply(lambda r: classify_bin(r,cond,None), axis=1)
    control_bins["behavior_collapsed"] = control_bins.apply(lambda r: classify_bin_collapsed(r,cond,None), axis=1)
    control_bins["change_direction"] = control_bins.apply(lambda r: get_change_direction(r,cond), axis=1)

    # probabilities
    def get_p(df, col):
        vc = df[col].value_counts()
        n = len(df)
        return {k: vc.get(k,0)/n for k in vc.index}

    control_p_sub = get_p(control_bins,"behavior")
    control_p_collapsed = get_p(control_bins,"behavior_collapsed")
    control_p_dir = get_p(control_bins,"change_direction")

    # ------------------------
    # STATISTICS
    # ------------------------
    for tid in trans_counts_sub.keys():
        n = sum(trans_counts_sub[tid])

        for level, counts_dict, control_p in [
            ("subcompartment", trans_counts_sub[tid], control_p_sub),
            ("collapsed", trans_counts_collapsed[tid], control_p_collapsed),
            ("direction", trans_counts_dir[tid], control_p_dir)
        ]:

            for cat in counts_dict.index:
                k = counts_dict.get(cat,0)
                p0 = control_p.get(cat,0)

                test = binomtest(k, n, p=p0)

                odds_trans = (k+epsilon)/(n-k+epsilon)
                odds_control = (p0+epsilon)/(1-p0+epsilon)

                stats_results.append({
                    "condition": cond,
                    "translocation": trans_name_map.get(tid,tid),
                    "level": level,
                    "category": cat,
                    "fraction": k/n,
                    "expected_fraction": p0,
                    "p_value": test.pvalue,
                    "log2_odds_ratio": np.log2(odds_trans/odds_control)
                })

    # ------------------------
    # ADD CONTROL TO PLOTS
    # ------------------------
    def add_control(results, control_bins, col):
        vc = control_bins[col].value_counts(normalize=True)
        return results + [{
            "transloc_id":"Genome control",
            **vc.to_dict()
        }]

    results_sub = add_control(results_sub, control_bins, "behavior")
    results_collapsed = add_control(results_collapsed, control_bins, "behavior_collapsed")
    results_direction = add_control(results_direction, control_bins, "change_direction")

    # ------------------------
    # DATAFRAMES
    # ------------------------
    def make_df(results):
        df = pd.DataFrame(results).fillna(0)
        df["transloc_id"] = df["transloc_id"].map(lambda x: trans_name_map.get(x,x))
        return df

    final_sub = make_df(results_sub)
    final_collapsed = make_df(results_collapsed)
    final_direction = make_df(results_direction)

    # ------------------------
    # PLOTTING
    # ------------------------
    def plot(df, ycols, filename, title, colors=None, use_colormap=False):
        fig, ax = plt.subplots(figsize=(8,5))
        if use_colormap:
            df.plot(
                x="transloc_id",
                kind="bar",
                stacked=True,
                y=ycols,
                colormap="tab10",
                ax=ax
            )
        else:
            df.plot(
                x="transloc_id",
                kind="bar",
                stacked=True,
                y=ycols,
                color=colors,
                ax=ax
            )
        ax.set_ylabel("Fraction of bins")
        ax.set_title(title)
        plt.tight_layout()
        fig.savefig(os.path.join(output_folder, filename), dpi=300)
        plt.close()

    plot(final_sub, ["retained","adopted","other"],
         f"{cond}_subcompartment.png",
         f"{cond} vs WT: Subcompartments", 
         use_colormap=True)

    plot(final_collapsed, ["retained","adopted","other"],
         f"{cond}_collapsed.png",
         f"{cond} vs WT: Collapsed A/B",
         use_colormap=True)

    plot(final_direction, ["A_to_B","B_to_A","retained"],
         f"{cond}_direction.png",
         f"{cond} vs WT: A↔B direction",
         colors=["#FF0800", "#1A43BF", "#cccccc"])

# ------------------------
# Save outputs
# ------------------------
pd.concat(all_bins).to_csv(os.path.join(output_folder,"bins_annotation.csv"), index=False)
pd.DataFrame(stats_results).to_csv(os.path.join(output_folder,"all_binomial_stats.csv"), index=False)

print("DONE: plots + full statistics generated")
