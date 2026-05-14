# Figure --> Script Map

This document maps every figure in the thesis to the script that
produced it. 

## Chapter 4: Results

### 4.1 Translocated Chromosomes Are Positioned Closer Together in T1 and C1

| Figure | Description | Script | Output file |
|--------|-------------|--------|-------------|
| 4.1.1  | Chrom3D 3D genome models for t(3;17)  | Chrom3D/ChimeraX | *.cmm |
| 4.1.2  | Chrom3D 3D genome models for t(6;19)  | Chrom3D/ChimeraX | *.cmm |
| 4.1.3  | Chrom3D 3D genome models for t(2;10)  | Chrom3D/ChimeraX | *.cmm |
| 4.1.4 | Centroid distance boxplots for t(3;17), t(6;19), t(2;10) | `08_chromosome_centroid_distance.py` | `chromosome_distance_boxplot*.png` |

### 4.2 Translocations Bring Distant Genes into Spatial Proximity

#### 4.2.1 Gene-pair distances (interchromosomal)

| Figure | Description | Script | Output file |
|--------|-------------|--------|-------------|
| 4.2.1 | Interchromosomal gene-pair distance, WT vs T1 (translocated and control) | `04_plot_results.py` | `T1_transloc_scatter.png`, `T1_control_scatter.png` |
| 4.2.2 | Interchromosomal gene-pair distance, WT vs C1 (translocated and control) | `04_plot_results.py` | `C1_transloc_scatter.png`, `C1_control_scatter.png` |
| 4.2.3 | Interchromosomal gene-pair distance, WT vs T1 (karyotype translocated and control) | `04_plot_results.py` | `T1_transloc_scatter.png`, `T1_control_scatter.png` |
| 4.2.4 | Interchromosomal gene-pair distance, WT vs C1 (karyotype translocated and control) | `04_plot_results.py` | `C1_transloc_scatter.png`, `C1_control_scatter.png` |

#### 4.2.2 Gene-pair distances (intrachromosomal)

| Figure | Description | Script | Output file |
|--------|-------------|--------|-------------|
| 4.2.5 | Intrachromosomal gene-pair distance, WT vs T1 (translocated and control) | `04_plot_results.py` | `T1_transloc_scatter.png`, `T1_control_scatter.png` |
| 4.2.6 | Intrachromosomal gene-pair distance, WT vs C1 (translocated and control) | `04_plot_results.py` | `C1_transloc_scatter.png`, `C1_control_scatter.png` |
| 4.2.7 | Intrachromosomal gene-pair distance, WT vs T1 (karyotype translocated and control) | `04_plot_results.py` | `T1_transloc_scatter.png`, `T1_control_scatter.png` |
| 4.2.8 | Intrachromosomal gene-pair distance, WT vs C1 (karyotype translocated and control) | `04_plot_results.py` | `C1_transloc_scatter.png`, `C1_control_scatter.png` |

#### 4.2.3 Radial distances to nuclear centre

| Figure | Description | Script | Output file |
|--------|-------------|--------|-------------|
| 4.2.9 | Distance to nuclear centre, translocated genes T1 and control | `07_plot_distance_to_center.py` | `T1_inside_transloc_scatter.png`, `T1_control_scatter.png` |
| 4.2.10 | Distance to nuclear centre, translocated genes C1 and control | `07_plot_distance_to_center.py` | `C1_inside_transloc_scatter.png`, `C1_control_scatter.png` |
| 4.2.11 | Distance to nuclear centre, karyotype translocated genes T1 and control | `07_plot_distance_to_center.py` | `T1_inside_transloc_scatter.png`, `T1_control_scatter.png` |
| 4.2.12 | Distance to nuclear centre, karyotype translocated genes C1 and control | `07_plot_distance_to_center.py` | `C1_inside_transloc_scatter.png`, `C1_control_scatter.png` |

### 4.3 Translocations Reorganize Subcompartment Identity

| Figure | Description | Script | Output file |
|--------|-------------|--------|-------------|
| 4.3.1 | Compartment classification (T1 and C1) | `10_subcompartment_classification.py` | `T1_collapsed.png`, `C1_collapsed.png` |
| 4.3.2 | A<-->B switching direction (T1 and C1) | `10_subcompartment_classification.py` | `T1_direction.png`, `C1_direction.png` |
| 4.3.3 | Subcompartment classification (T1 and C1) | `10_subcompartment_classification.py` | `T1_subcompartment.png`, `C1_subcompartment.png` |
| 4.3.4 | Subcompartment tracks: t(2;10) across WT/C1 | `11_subcompartment_tracks.py` | `subcompartment_t2_10_*.png` |
| 4.3.5 | Subcompartment tracks: t(6;19) across WT/T1/C1 | `11_subcompartment_tracks.py` | `subcompartment_t6_19_*.png` |
| 4.3.6 | Subcompartment tracks: t(3;17) across WT/T1/C1 | `11_subcompartment_tracks.py` | `subcompartment_t3_17_*.png` |
| 4.3.7 | Alluvial: Der(3)t(3;17) fragment and body | `12_subcompartment_alluvial.R` | `alluvial_chr3_body_Der_3_t_3_17.jpeg`, `alluvial_chr17_fragment_Der_3_t_3_17.jpeg` |
| 4.3.8 | Alluvial: Der(17)t(3;17) fragment and body | `12_subcompartment_alluvial.R` | `alluvial_chr17_body_Der_17_t_3_17.jpeg`, `alluvial_chr3_fragment_Der_17_t_3_17.jpeg` |
| 4.3.9 | Alluvial: t(6;19) fragment and body | `12_subcompartment_alluvial.R` | `alluvial_chr6_body_t_6_19.jpeg`, `alluvial_chr19_fragment_t_6_19.jpeg` |
| 4.3.10 | Alluvial: t(2;10) fragment and body | `12_subcompartment_alluvial.R` | `alluvial_chr10_body_t_2_10.jpeg`, `alluvial_chr2_fragment_t_2_10.jpeg` |
| 4.3.11 | Radial trajectory of translocated segments | `09_nuclear_radial_positioning.py` | `trajectory_radial_distance.png`, `trajectory_chr3_chr17_and_fragments.png` |

### 4.4 Gene Expression

| Figure | Description | Script | Output file |
|--------|-------------|--------|-------------|
| 4.4.1  | Gene expression across compartment behavior (T1 and C1) | `14_gene_bin_expression.py`  | `T1_behavior_collapsed_DE_fraction.png`, `C1_behavior_collapsed_DE_fraction.png` |
| 4.4.2  | Gene expression across A<-->B switching direction behavior (T1 and C1) | `14_gene_bin_expression.py`  | `T1_change_direction_DE_fraction.png`, `C1_change_direction_DE_fraction.png` |
| 4.4.3  | Gene expression across subcompartment behavior (T1 and C1) | `14_gene_bin_expression.py`  | `T1_behavior_DE_fraction.png`, `C1_behavior_DE_fraction.png` |
| 4.4.4 | DE proportions: translocated vs genome (T1 and C1) | `15_de_enrichment.py` | `T1_grouped_DE_proportions.png`, `C1_grouped_DE_proportions.png` |
| 4.4.5 | Volcano plot, WT vs T1/C1 (all translocated) | `16_volcano_plots.py` | `volcano_WT_vs_T1_highlighted.png`, `volcano_WT_vs_C1_highlighted.png` |
| 4.4.6 | Volcano: Der(17)t(3;17) genes in T1 and C1 | `17_translocation_volcano_plots.py` | `volcano_WT_vs_T1_Der17t3_17.png`, `volcano_WT_vs_C1_Der17t3_17.png` |
| 4.4.7 | Volcano: Der(3)t(3;17) genes in T1 and C1 | `17_translocation_volcano_plots.py` | `volcano_WT_vs_T1_Der3t3_17.png`, `volcano_WT_vs_C1_Der3t3_17.png` |
| 4.4.8 | GO BP enrichment, Hi-C and karyotype translocations | `18_go_enrichment.R` | `C1_GO_BP_dotplot.png`, `verify_C1_GO_BP_dotplot.png` |
