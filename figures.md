# Figure --> Script Map

This document maps every figure in the thesis to the script that
produced it. Output paths are relative to each script's configured
output directory.

## Chapter 4: Results

### 4.1 Spatial proximity of translocated regions

| Figure | Description | Script | Output file |
|--------|-------------|--------|-------------|
| 4.1.1 | Centroid distance boxplots for t(3;17), t(6;19), t(2;10) | `08_chromosome_centroid_distance.py` | `chromosome_distance_boxplot.png` |
| 4.1.2 | Gene-pair distance, WT vs T1 (translocated) | `04_plot_results.py` | `T1_transloc_scatter.png` |
| 4.1.3 | Gene-pair distance, WT vs C1 (translocated) | `04_plot_results.py` | `C1_transloc_scatter.png` |
| 4.1.4 | Gene-pair distance, WT vs T1 (control) | `04_plot_results.py` | `T1_control_scatter.png` |
| 4.1.5 | Gene-pair distance, WT vs C1 (control) | `04_plot_results.py` | `C1_control_scatter.png` |

### 4.2 Radial nuclear positioning

| Figure | Description | Script | Output file |
|--------|-------------|--------|-------------|
| 4.2.1 | Distance to nuclear centre, translocated genes T1 | `07_plot_distance_to_center.py` | `T1_inside_transloc_scatter.png` |
| 4.2.2 | Distance to nuclear centre, translocated genes C1 | `07_plot_distance_to_center.py` | `C1_inside_transloc_scatter.png` |
| 4.2.3 | Radial trajectory of translocated segments | `09_nuclear_radial_positioning.py` | `trajectory_radial_distance.png` |
| 4.2.4 | Radial trajectory of chr3/chr17 with fragments | `09_nuclear_radial_positioning.py` | `trajectory_chr3_chr17_and_fragments.png` |

### 4.3 Subcompartment reorganisation

| Figure | Description | Script | Output file |
|--------|-------------|--------|-------------|
| 4.3.1 | T1 subcompartment classification | `10_subcompartment_classification.py` | `T1_subcompartment.png` |
| 4.3.2 | C1 subcompartment classification | `10_subcompartment_classification.py` | `C1_subcompartment.png` |
| 4.3.3 | T1 collapsed A/B classification | `10_subcompartment_classification.py` | `T1_collapsed.png` |
| 4.3.4 | C1 collapsed A/B classification | `10_subcompartment_classification.py` | `C1_collapsed.png` |
| 4.3.5 | A<-->B switching direction (T1) | `10_subcompartment_classification.py` | `T1_direction.png` |
| 4.3.6 | A<-->B switching direction (C1) | `10_subcompartment_classification.py` | `C1_direction.png` |
| 4.3.7a | Alluvial: Der(3)t(3;17) whole translocation | `12_subcompartment_alluvial.R` | `alluvial_Der_3_t_3_17_whole_translocation.jpeg` |
| 4.3.7b | Alluvial: Der(17)t(3;17) whole translocation | `12_subcompartment_alluvial.R` | `alluvial_Der_17_t_3_17_whole_translocation.jpeg` |
| 4.3.8 | Subcompartment tracks: Der(17)t(3;17) across WT/T1/C1 | `11_subcompartment_tracks.py` | `subcompartment_Der17t3_17_*.png` |

### 4.4 Gene expression

| Figure | Description | Script | Output file |
|--------|-------------|--------|-------------|
| 4.4.1 | Volcano plot, WT vs T1 (all translocated) | `16_volcano_plots.py` | `volcano_WT_vs_T1_highlighted.png` |
| 4.4.2 | Volcano plot, WT vs C1 (all translocated) | `16_volcano_plots.py` | `volcano_WT_vs_C1_highlighted.png` |
| 4.4.3 | Volcano: Der(17)t(3;17) genes in C1 | `17_translocation_volcano_plots.py` | `volcano_WT_vs_C1_Der17t3_17.png` |
| 4.4.4 | Volcano: Der(3)t(3;17) genes in C1 | `17_translocation_volcano_plots.py` | `volcano_WT_vs_C1_Der3t3_17.png` |
| 4.4.5 | DE proportions: translocated vs genome (T1) | `15_de_enrichment.py` | `T1_grouped_DE_proportions.png` |
| 4.4.6 | DE proportions: translocated vs genome (C1) | `15_de_enrichment.py` | `C1_grouped_DE_proportions.png` |
| 4.4.7 | Subcompartment behaviour vs DE (T1) | `14_gene_bin_expression.py` | `T1_behavior_DE_fraction.png` |
| 4.4.8 | Subcompartment behaviour vs DE (C1) | `14_gene_bin_expression.py` | `C1_behavior_DE_fraction.png` |
| 4.4.9 | GO BP enrichment, T1 | `18_go_enrichment.R` | `T1_GO_BP_dotplot.png` |
| 4.4.10 | GO BP enrichment, C1 | `18_go_enrichment.R` | `C1_GO_BP_dotplot.png` |

## Tables

| Table | Description | Script | Output file |
|-------|-------------|--------|-------------|
| 4.3.3 | Binomial test results for subcompartment classification | `10_subcompartment_classification.py` | `all_binomial_stats.csv` |
| 4.4.x | Top DE translocated genes | `16_volcano_plots.py` | `top_up_transloc_genes_*.tsv`, `top_down_transloc_genes_*.tsv` |
