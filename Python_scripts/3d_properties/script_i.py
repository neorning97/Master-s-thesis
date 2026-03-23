import pandas as pd
import matplotlib.pyplot as plt

# Load subcompartments
sub = pd.read_csv("/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/subcompartments/GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph", sep="\t")

# Choose a chromosome to plot (for clarity)
chr_to_plot = "chr10"
df = sub[sub["chrom"] == chr_to_plot]

# Define subcompartment colors
subcomp_colors = {
    "A0": "#fcaeae",
    "A1": "#fb6a6a",
    "A2": "#ef3b2c",
    "A3": "#99000d",
    "B0": "#deebf7",
    "B1": "#9ecae1",
    "B2": "#4292c6",
    "B3": "#084594"
}

fig, ax = plt.subplots(figsize=(12,2))

for _, row in df.iterrows():
    ax.fill_between(
        [row["start"], row["end"]],
        0, 1,
        color=subcomp_colors[row["MCF10A_WT.state"]],
        edgecolor="none"
    )

ax.set_xlim(df["start"].min(), df["end"].max())
ax.set_ylim(0, 1)
ax.set_yticks([])
ax.set_xlabel(f"{chr_to_plot} position (bp)")
ax.set_title("Subcompartment track (WT)")

# Add a legend
for state, color in subcomp_colors.items():
    ax.plot([], [], color=color, label=state, linewidth=6)
ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()
