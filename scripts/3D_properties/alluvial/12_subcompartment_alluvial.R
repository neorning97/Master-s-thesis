# ============================================================================
# 12_subcompartment_alluvial.R
# ============================================================================
#
# What this script does
# ---------------------
# Produces alluvial (Sankey) diagrams showing how genomic bins flow between
# chromatin subcompartments across conditions (WT -> T1 -> C1 or WT -> C1).
#
# Two types of plot are produced for each translocation:
#
#   1. Per-segment plots — one plot per individual BED row (fragment or body)
#      showing transitions within that specific genomic segment.
#
#   2. Whole-translocation plots — one plot per translocation (e.g.
#      Der(17)t(3;17)) combining both the fragment and body rows into a single
#      region. This gives an overview of the full translocation event.
#
# Translocations and conditions covered:
#
#   Der(17)t(3;17) — chr3:0–58.6Mb + chr17:22.7–83.3Mb — WT, T1, C1
#   Der(3)t(3;17)  — chr3:58.6–198.3Mb + chr17:0–22.7Mb — WT, T1, C1
#   t(6;19)        — chr6:0–170.8Mb + chr19:0–31.9Mb    — WT, T1, C1
#   t(2;10)        — chr2:0–131.7Mb + chr10:39.8–133.8Mb — WT and C1 ONLY
#                    (t(2;10) does not occur in T1)
#
# Segments that appear in both the T1 and C1 BED files are only plotted once.
#
# Why an alluvial plot?
# ---------------------
# A bar chart shows the marginal distribution of subcompartments per condition
# but cannot show which bins changed state and where they went. The alluvial
# diagram makes transitions explicit: the width of each ribbon is proportional
# to the number of 100 kb bins that follow that particular transition path.
#
# Colour scheme
# -------------
# A-compartment = shades of red (A0 = lightest, A3 = darkest).
# B-compartment = shades of blue (B0 = lightest, B3 = darkest).
# Matches the colour scheme used in scripts 10 and 11.
#
# Usage
# -----
#   1. Edit the CONFIG section below to point to your files.
#   2. Run: Rscript 12_subcompartment_alluvial.R
#      or source the script from RStudio.
#
# Dependencies
# ------------
#   tidyverse, ggalluvial
#   Install with: install.packages(c("tidyverse", "ggalluvial"))
#
# ============================================================================


# ============================================================================
# CONFIG – edit all paths and settings here before running
# ============================================================================

# Subcompartment bedGraph file (100 kb resolution, WT + T1 + C1)
INPUT_FILE <- "/path/to/GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph"

# Translocation BED files — one row per segment, paired by transloc_id.
TRANSLOC_BEDS <- list(
  T1 = "/path/to/verify_T1_translocations.bed",
  C1 = "/path/to/verify_C1_translocations.bed"
)

# Maps the 'label' column in the BED files to a human-readable segment name,
# used in plot titles and output filenames.
LABEL_NAMES <- c(
  "T1_transloc_1"         = "chr3 fragment – Der(17)t(3;17)",
  "T1_transloc_1_partner" = "chr17 body – Der(17)t(3;17)",
  "T1_transloc_2"         = "chr3 body – Der(3)t(3;17)",
  "T1_transloc_2_partner" = "chr17 fragment – Der(3)t(3;17)",
  "T1_transloc_3"         = "chr6 body – t(6;19)",
  "T1_transloc_3_partner" = "chr19 fragment – t(6;19)",
  "C1_transloc_1"         = "chr2 fragment – t(2;10)",
  "C1_transloc_1_partner" = "chr10 body – t(2;10)",
  "C1_transloc_2"         = "chr3 fragment – Der(17)t(3;17)",
  "C1_transloc_2_partner" = "chr17 body – Der(17)t(3;17)",
  "C1_transloc_3"         = "chr3 body – Der(3)t(3;17)",
  "C1_transloc_3_partner" = "chr17 fragment – Der(3)t(3;17)",
  "C1_transloc_4"         = "chr6 body – t(6;19)",
  "C1_transloc_4_partner" = "chr19 fragment – t(6;19)"
)

# Maps transloc_id values to whole-translocation names.
# Each transloc_id covers two rows (fragment + body); this name is used for
# the combined whole-translocation plot.
WHOLE_TRANSLOC_NAMES <- c(
  # T1 BED transloc_ids
  "T1" = "Der(17)t(3;17)",
  "T2" = "Der(3)t(3;17)",
  "T3" = "t(6;19)",
  # C1 BED transloc_ids
  "T1" = "t(2;10)",      # C1 T1 = t(2;10) — handled per-BED below
  "T2" = "Der(17)t(3;17)",
  "T3" = "Der(3)t(3;17)",
  "T4" = "t(6;19)"
)

# Translocation name maps per BED file (same as in script 11)
TRANSLOC_NAME_MAPS <- list(
  T1 = c("T1" = "Der(17)t(3;17)", "T2" = "Der(3)t(3;17)", "T3" = "t(6;19)"),
  C1 = c("T1" = "t(2;10)", "T2" = "Der(17)t(3;17)", "T3" = "Der(3)t(3;17)", "T4" = "t(6;19)")
)

# Translocations that should only be plotted for WT and C1 (not T1).
# t(2;10) only occurs in C1, so including T1 would show an unaffected genome.
WJCT_ONLY_CONDITIONS <- c("t(2;10)")   # WT + C1 only

# Where to save output plots (created automatically)
OUTPUT_DIR <- "/path/to/results/subcompartment_alluvial"

# Output dimensions and resolution
PLOT_WIDTH  <- 8    # inches
PLOT_HEIGHT <- 6    # inches
PLOT_DPI    <- 300

# ============================================================================


# ============================================================================
# Load libraries
# ============================================================================

library(tidyverse)
library(ggalluvial)


# ============================================================================
# Colour palette and factor levels (consistent across all plots)
# ============================================================================

subcomp_colors <- c(
  "A3" = "#CB181D", "A2" = "#FB6A4A", "A1" = "#FCAE91", "A0" = "#FEE5D9",
  "B0" = "#DEEBF7", "B1" = "#9ECAE1", "B2" = "#4292C6", "B3" = "#084594"
)

# Ordered from most inactive (B3) to most active (A3)
sub_levels <- c("B3", "B2", "B1", "B0", "A0", "A1", "A2", "A3")


# ============================================================================
# Load subcompartment data
# ============================================================================

cat("Loading subcompartment data...\n")
sub <- read_tsv(INPUT_FILE, show_col_types = FALSE) %>%
  mutate(
    MCF10A_WT.state = factor(MCF10A_WT.state, levels = sub_levels),
    MCF10A_T1.state = factor(MCF10A_T1.state, levels = sub_levels),
    MCF10A_C1.state = factor(MCF10A_C1.state, levels = sub_levels)
  ) %>%
  drop_na()

cat(sprintf("Loaded %d bins total.\n", nrow(sub)))


# ============================================================================
# Helper functions
# ============================================================================

make_safe_name <- function(name) {
  # Convert a human-readable name to a safe filename (no special characters)
  name %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_remove("_$")
}


build_alluvial_plot <- function(flow, title, axes = c("WT", "T1", "C1"),
                                subcomp_colors) {
  # Determine which axes are present in the flow data
  # (for WT+C1 only plots, axis2 = C1 directly)
  n_axes <- length(axes)
  
  if (n_axes == 3) {
    p <- ggplot(flow,
                aes(axis1 = MCF10A_WT.state,
                    axis2 = MCF10A_T1.state,
                    axis3 = MCF10A_C1.state,
                    y = n))
    x_limits <- c("WT", "T1", "C1")
  } else {
    # Two-condition plot: WT and C1 only
    p <- ggplot(flow,
                aes(axis1 = MCF10A_WT.state,
                    axis2 = MCF10A_C1.state,
                    y = n))
    x_limits <- c("WT", "C1")
  }
  
  p +
    geom_alluvium(aes(fill = after_stat(stratum)), alpha = 0.7) +
    geom_stratum(aes(fill = after_stat(stratum)), color = "black", size = 0.2) +
    scale_fill_manual(values = subcomp_colors, name = "Subcompartment") +
    scale_x_discrete(limits = x_limits, labels = x_limits) +
    # Reverse y-axis so A-compartment (active) appears at the top
    scale_y_reverse(breaks = NULL) +
    labs(title = title, x = "Cell line", y = "Number of 100 kb bins") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      panel.grid      = element_blank(),
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.text.x     = element_text(color = "black"),
      axis.title      = element_text(face = "bold"),
      plot.title      = element_text(face = "bold", hjust = 0.5)
    )
}


save_plot <- function(p, name, output_dir, width, height, dpi) {
  out_path <- file.path(output_dir, sprintf("alluvial_%s.jpeg", make_safe_name(name)))
  ggsave(out_path, p, width = width, height = height, dpi = dpi)
  cat(sprintf("  Saved: %s\n", out_path))
}


filter_bins <- function(sub, regions) {
  # Filter subcompartment bins to any row matching one of the given regions.
  # regions is a data frame with columns: chrom, start, end
  map_dfr(seq_len(nrow(regions)), function(i) {
    sub %>% filter(
      chrom >= regions$chrom[i] &   # same chromosome
        chrom <= regions$chrom[i] &
        start >= regions$start[i] &
        end   <= regions$end[i]
    )
  }) %>%
    distinct(chrom, start, end, .keep_all = TRUE)
}


# ============================================================================
# Load BED files and build segment tables
# ============================================================================

cat("\nLoading translocation segments from BED files...\n")

# Individual segments (one row per BED row)
all_segments <- map_dfr(names(TRANSLOC_BEDS), function(cond_key) {
  read_tsv(TRANSLOC_BEDS[[cond_key]], show_col_types = FALSE) %>%
    mutate(source_condition = cond_key)
}) %>%
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  mutate(segment_name = coalesce(LABEL_NAMES[label], label)) %>%
  distinct(chrom, start, end, .keep_all = TRUE)

cat(sprintf("Found %d unique individual segments.\n", nrow(all_segments)))

# Whole-translocation regions: combine fragment + body rows per transloc_id
# by taking the union of their coordinates on each chromosome separately.
whole_translocs <- map_dfr(names(TRANSLOC_BEDS), function(cond_key) {
  bed      <- read_tsv(TRANSLOC_BEDS[[cond_key]], show_col_types = FALSE) %>%
    mutate(start = as.integer(start), end = as.integer(end))
  name_map <- TRANSLOC_NAME_MAPS[[cond_key]]
  
  map_dfr(unique(bed$transloc_id), function(tid) {
    rows <- bed %>% filter(transloc_id == tid)
    transloc_name <- name_map[tid]
    
    # Return one row per chromosome involved in this translocation
    rows %>%
      group_by(chrom) %>%
      summarise(start = min(start), end = max(end), .groups = "drop") %>%
      mutate(
        transloc_name    = transloc_name,
        source_condition = cond_key
      )
  })
}) %>%
  # Drop duplicate chromosome regions for the same translocation
  distinct(transloc_name, chrom, start, end, .keep_all = TRUE)

cat(sprintf("Found %d unique whole-translocation chromosome regions.\n",
            nrow(whole_translocs)))


# ============================================================================
# Part 1: Per-segment alluvial plots
# ============================================================================

cat("\n--- Generating per-segment alluvial plots ---\n")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(all_segments))) {
  seg  <- all_segments[i, ]
  name <- seg$segment_name
  
  # Determine whether this is a t(2;10) segment (WT + C1 only)
  is_wt_c1_only <- any(map_lgl(WJCT_ONLY_CONDITIONS, ~ str_detect(name, fixed(.x))))
  axes <- if (is_wt_c1_only) c("WT", "C1") else c("WT", "T1", "C1")
  
  sub_region <- sub %>%
    filter(chrom == seg$chrom & start >= seg$start & end <= seg$end)
  
  if (nrow(sub_region) == 0) {
    cat(sprintf("  No bins for %s, skipping.\n", name)); next
  }
  cat(sprintf("  %s: %d bins (%s)\n", name, nrow(sub_region),
              paste(axes, collapse = " + ")))
  
  flow <- if (is_wt_c1_only) {
    sub_region %>% count(MCF10A_WT.state, MCF10A_C1.state)
  } else {
    sub_region %>% count(MCF10A_WT.state, MCF10A_T1.state, MCF10A_C1.state)
  }
  
  p <- build_alluvial_plot(flow, title = name, axes = axes,
                           subcomp_colors = subcomp_colors)
  save_plot(p, name, OUTPUT_DIR, PLOT_WIDTH, PLOT_HEIGHT, PLOT_DPI)
}


# ============================================================================
# Part 2: Whole-translocation alluvial plots
# ============================================================================
# Each whole-translocation plot combines all chromosomal regions involved in
# one translocation event (e.g. both the chr3 fragment and the chr17 body
# for Der(17)t(3;17)) into a single count of transitions.

cat("\n--- Generating whole-translocation alluvial plots ---\n")

for (trans_name in unique(whole_translocs$transloc_name)) {
  regions <- whole_translocs %>% filter(transloc_name == trans_name)
  
  is_wt_c1_only <- trans_name %in% WJCT_ONLY_CONDITIONS
  axes <- if (is_wt_c1_only) c("WT", "C1") else c("WT", "T1", "C1")
  
  # Collect bins from all chromosomal regions of this translocation
  sub_all <- map_dfr(seq_len(nrow(regions)), function(i) {
    sub %>% filter(
      chrom == regions$chrom[i] &
        start >= regions$start[i] &
        end   <= regions$end[i]
    )
  }) %>%
    distinct(chrom, start, end, .keep_all = TRUE)
  
  if (nrow(sub_all) == 0) {
    cat(sprintf("  No bins for %s — skipping.\n", trans_name)); next
  }
  cat(sprintf("  %s (whole): %d bins across %d chromosome region(s) (%s)\n",
              trans_name, nrow(sub_all), nrow(regions),
              paste(axes, collapse = " + ")))
  
  flow <- if (is_wt_c1_only) {
    sub_all %>% count(MCF10A_WT.state, MCF10A_C1.state)
  } else {
    sub_all %>% count(MCF10A_WT.state, MCF10A_T1.state, MCF10A_C1.state)
  }
  
  plot_title <- paste0(trans_name, " (whole translocation)")
  p <- build_alluvial_plot(flow, title = plot_title, axes = axes,
                           subcomp_colors = subcomp_colors)
  save_plot(p, plot_title, OUTPUT_DIR, PLOT_WIDTH, PLOT_HEIGHT, PLOT_DPI)
}

cat("\nDone. All plots saved to:", OUTPUT_DIR, "\n")
