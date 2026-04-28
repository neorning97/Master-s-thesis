# ============================================================================
# 12_subcompartment_alluvial.R
# ============================================================================
#
# This script makes "alluvial" plots (also called Sankey diagrams) showing
# how 100 kb genomic bins flow between subcompartments across conditions.
#
# What's an alluvial plot?
# - Imagine a bar chart for each condition (WT, T1, C1) showing how many
#   bins are in each subcompartment.
# - Now imagine RIBBONS connecting the bars, showing which bins moved from
#   one subcompartment to another between conditions.
# - Wider ribbon = more bins took that path.
# - It lets you see WHICH bins changed and where they went.
#
# What we plot:
# 1. Per-segment plots: one plot per individual fragment or body piece.
# 2. Whole-translocation plots: one plot combining the fragment + body
#    of each translocation into a single overall view.
#
# The translocations we plot:
#   Der(17)t(3;17) — chr3 fragment + chr17 body — WT, T1, C1
#   Der(3)t(3;17)  — chr17 fragment + chr3 body — WT, T1, C1
#   t(6;19)        — chr19 fragment + chr6 body — WT, T1, C1
#   t(2;10)        — chr2 fragment + chr10 body — WT and C1 ONLY
#                    (because t(2;10) doesn't exist in T1)
#
# Edit the file paths below before running.
#
# Required packages: tidyverse, ggalluvial
# Install them with: install.packages(c("tidyverse", "ggalluvial"))
# ============================================================================


# ============================================================================
# CONFIG SECTION - Edit these settings before running!
# ============================================================================

# -----------------------------------------------------------------------------
# Subcompartment file
# -----------------------------------------------------------------------------
# Tells us what subcompartment each 100 kb bin is in, for WT/T1/C1
INPUT_FILE <- "/path/to/GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph"

# -----------------------------------------------------------------------------
# BED files for translocations - one for each condition
# -----------------------------------------------------------------------------
T1_BED <- "/path/to/verify_T1_translocations.bed"
C1_BED <- "/path/to/verify_C1_translocations.bed"

# Group them in a list so we can loop through them later.
# A "list" in R is like a dictionary, it can hold named values of
# different types.
TRANSLOC_BEDS <- list(
  T1 = T1_BED,
  C1 = C1_BED
)

# -----------------------------------------------------------------------------
# Map from "label" column in the BED files -> readable segment name
# -----------------------------------------------------------------------------
# This is a "named character vector" in R, basically a key-->value lookup.
# We use it for the per-segment plot titles.
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

# -----------------------------------------------------------------------------
# Maps from transloc_id -> readable translocation name (one map per BED file)
# -----------------------------------------------------------------------------
# We need separate maps because the same transloc_id (like "T1") means
# different things in different BED files.
T1_NAME_MAP <- c(
  "T1" = "Der(17)t(3;17)",
  "T2" = "Der(3)t(3;17)",
  "T3" = "t(6;19)"
)

C1_NAME_MAP <- c(
  "T1" = "t(2;10)",
  "T2" = "Der(17)t(3;17)",
  "T3" = "Der(3)t(3;17)",
  "T4" = "t(6;19)"
)

# Group them in a list so we can pick the right map per BED file.
TRANSLOC_NAME_MAPS <- list(
  T1 = T1_NAME_MAP,
  C1 = C1_NAME_MAP
)

# Note: the original script also defined a WHOLE_TRANSLOC_NAMES vector with
# duplicate names ("T1" appeared twice etc.). When R sees duplicate names,
# only the FIRST one is used on lookup, so it would silently give wrong
# results. The script doesn't actually need it (it uses TRANSLOC_NAME_MAPS
# instead) so I haven't included it here.

# -----------------------------------------------------------------------------
# Translocations to plot for WT and C1 only (skipping T1)
# -----------------------------------------------------------------------------
# t(2;10) doesn't exist in T1, so showing the T1 column would be misleading.
WT_C1_ONLY <- c("t(2;10)")

# -----------------------------------------------------------------------------
# Output folder for plots
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "/path/to/results/subcompartment_alluvial"

# -----------------------------------------------------------------------------
# Plot dimensions
# -----------------------------------------------------------------------------
PLOT_WIDTH  <- 8     # inches
PLOT_HEIGHT <- 6     # inches
PLOT_DPI    <- 300   # dots per inch (resolution)


# ============================================================================
# Load libraries
# ============================================================================
# tidyverse gives us things like dplyr (data manipulation), ggplot2 (plots),
# readr (reading files), stringr (string functions).
# ggalluvial adds the alluvial plot layer for ggplot2.

library(tidyverse)
library(ggalluvial)


# ============================================================================
# Colour palette and ordering
# ============================================================================

# Each subcompartment gets a colour:
# Reds for A (active), blues for B (inactive)
subcomp_colors <- c(
  "A3" = "#CB181D",   # darkest red
  "A2" = "#FB6A4A",
  "A1" = "#FCAE91",
  "A0" = "#FEE5D9",   # lightest red
  "B0" = "#DEEBF7",   # lightest blue
  "B1" = "#9ECAE1",
  "B2" = "#4292C6",
  "B3" = "#084594"    # darkest blue
)

# When we make the subcompartment column a "factor", R uses this order
# in the plots (instead of alphabetical order).
# We order from most inactive (B3) at the bottom to most active (A3) at the top.
sub_levels <- c("B3", "B2", "B1", "B0", "A0", "A1", "A2", "A3")


# ============================================================================
# Helper function: Make a "safe" filename
# ============================================================================
# Translocation names contain characters like "(", ")", and ";" that aren't
# great in filenames. This replaces them with underscores.

make_safe_name <- function(name) {

  # str_replace_all replaces every match of a pattern.
  # The pattern "[^A-Za-z0-9]+" matches one or more characters that are
  # NOT letters or numbers. So things like "(", ")", ";", "-", " ", etc.
  safe <- str_replace_all(name, "[^A-Za-z0-9]+", "_")

  # Remove a trailing underscore if there is one
  safe <- str_remove(safe, "_$")

  return(safe)
}


# ============================================================================
# Helper function: Build an alluvial plot
# ============================================================================
# Takes a "flow" data frame (counts of bins for each path through conditions)
# and returns a ggplot object that can be saved or displayed.

build_alluvial_plot <- function(flow, title, axes) {

  # We use different aesthetics depending on whether we have 2 or 3 axes
  n_axes <- length(axes)

  if (n_axes == 3) {

    # 3 axes: WT -> T1 -> C1
    p <- ggplot(flow, aes(
      axis1 = MCF10A_WT.state,
      axis2 = MCF10A_T1.state,
      axis3 = MCF10A_C1.state,
      y = n
    ))
    x_limits <- c("WT", "T1", "C1")

  } else {

    # 2 axes: WT -> C1 (used for t(2;10))
    p <- ggplot(flow, aes(
      axis1 = MCF10A_WT.state,
      axis2 = MCF10A_C1.state,
      y = n
    ))
    x_limits <- c("WT", "C1")
  }

  # ---------------------------------------------------------------------------
  # Build the plot by adding layers with +
  # ---------------------------------------------------------------------------
  # geom_alluvium draws the ribbons between bars.
  # geom_stratum draws the bars themselves.
  # after_stat(stratum) gets the subcompartment label so we can colour by it.

  p <- p + geom_alluvium(aes(fill = after_stat(stratum)), alpha = 0.7)
  p <- p + geom_stratum(aes(fill = after_stat(stratum)),
                        color = "black", size = 0.2)

  # Use our colour palette and add a legend title
  p <- p + scale_fill_manual(values = subcomp_colors,
                             name = "Subcompartment")

  # Set the x-axis labels
  p <- p + scale_x_discrete(limits = x_limits, labels = x_limits)

  # Reverse the y-axis so A-compartment (active) appears on top.
  # breaks = NULL hides the y-axis numbers (they don't mean much here).
  p <- p + scale_y_reverse(breaks = NULL)

  # Add the title and axis labels
  p <- p + labs(
    title = title,
    x = "Cell line",
    y = "Number of 100 kb bins"
  )

  # Use the minimal theme as a starting point
  p <- p + theme_minimal(base_size = 14)

  # Customize various theme elements:
  # - put legend on the right
  # - hide grid lines and y-axis text/ticks
  # - make x-axis text black, axis titles bold, plot title bold and centred
  p <- p + theme(
    legend.position = "right",
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.x     = element_text(color = "black"),
    axis.title      = element_text(face = "bold"),
    plot.title      = element_text(face = "bold", hjust = 0.5)
  )

  return(p)
}


# ============================================================================
# Helper function: Save a plot to a JPEG file
# ============================================================================

save_plot <- function(p, name) {

  # Make a safe filename
  safe <- make_safe_name(name)

  # Build the full file path.
  # sprintf is like Python's format: %s gets replaced by the string.
  filename <- sprintf("alluvial_%s.jpeg", safe)
  out_path <- file.path(OUTPUT_DIR, filename)

  # Save the plot
  ggsave(out_path, p,
         width = PLOT_WIDTH,
         height = PLOT_HEIGHT,
         dpi = PLOT_DPI)

  cat(sprintf("  Saved: %s\n", out_path))
}


# ============================================================================
# Step 1: Load the subcompartment data
# ============================================================================

cat("Loading subcompartment data...\n")

# read_tsv reads a tab-separated file
# show_col_types = FALSE quiets the messages about column types
sub <- read_tsv(INPUT_FILE, show_col_types = FALSE)

# Convert the subcompartment columns to factors with our ordering.
# A "factor" is R's way of storing categorical data with a defined order.
sub$MCF10A_WT.state <- factor(sub$MCF10A_WT.state, levels = sub_levels)
sub$MCF10A_T1.state <- factor(sub$MCF10A_T1.state, levels = sub_levels)
sub$MCF10A_C1.state <- factor(sub$MCF10A_C1.state, levels = sub_levels)

# Drop rows with any NA (missing) values
sub <- drop_na(sub)

cat(sprintf("Loaded %d bins total.\n", nrow(sub)))


# ============================================================================
# Step 2: Load and combine all BED files into one big segments table
# ============================================================================
# We loop through each BED file, read it, and stack the rows into one table.

cat("\nLoading translocation segments from BED files...\n")

# Empty list - we'll fill it with one data frame per BED file
loaded_beds <- list()

# names(TRANSLOC_BEDS) gives us the keys: "T1" and "C1"
for (cond_key in names(TRANSLOC_BEDS)) {

  # Get the file path for this condition
  bed_path <- TRANSLOC_BEDS[[cond_key]]

  # Read the BED file
  bed_data <- read_tsv(bed_path, show_col_types = FALSE)

  # Add a column saying which condition this came from
  bed_data$source_condition <- cond_key

  # Add it to our list
  loaded_beds[[cond_key]] <- bed_data
}

# bind_rows stacks data frames on top of each other (like pd.concat in Python)
all_segments <- bind_rows(loaded_beds)

# Make sure start and end are integers
all_segments$start <- as.integer(all_segments$start)
all_segments$end   <- as.integer(all_segments$end)

# Add a "segment_name" column with the readable name.
# LABEL_NAMES[label] looks up each label in our LABEL_NAMES vector.
# coalesce() returns the first non-NA value. If a label isn't in the
# lookup, we use the label itself as a fallback.
all_segments$segment_name <- coalesce(
  LABEL_NAMES[all_segments$label],
  all_segments$label
)

# Remove duplicate rows (same chrom/start/end across both BED files)
all_segments <- distinct(all_segments, chrom, start, end, .keep_all = TRUE)

cat(sprintf("Found %d unique individual segments.\n", nrow(all_segments)))


# ============================================================================
# Step 3: Build the "whole translocation" regions table
# ============================================================================
# For each translocation, we want to combine its fragment + body rows
# into a single region per chromosome.
# This gives us, for example, "Der(17)t(3;17)" with one chr3 region and
# one chr17 region.

cat("\nBuilding whole-translocation regions...\n")

# We'll collect rows for the final data frame in this list
whole_transloc_pieces <- list()

# Loop through each BED file
for (cond_key in names(TRANSLOC_BEDS)) {

  # Read the BED file (again, it's a small file, that's fine)
  bed_path <- TRANSLOC_BEDS[[cond_key]]
  bed <- read_tsv(bed_path, show_col_types = FALSE)
  bed$start <- as.integer(bed$start)
  bed$end   <- as.integer(bed$end)

  # Get the name map for this condition
  name_map <- TRANSLOC_NAME_MAPS[[cond_key]]

  # Loop through each unique transloc_id in this BED file
  unique_tids <- unique(bed$transloc_id)

  for (tid in unique_tids) {

    # Get all rows for this translocation
    rows <- bed[bed$transloc_id == tid, ]

    # Look up the readable name
    transloc_name <- name_map[tid]

    # Group by chromosome and take min(start) and max(end) to get the
    # whole region this translocation covers on each chromosome.
    summary_per_chrom <- rows %>%
      group_by(chrom) %>%
      summarise(
        start = min(start),
        end = max(end),
        .groups = "drop"
      )

    # Add the translocation name and source condition columns
    summary_per_chrom$transloc_name <- transloc_name
    summary_per_chrom$source_condition <- cond_key

    # Add it to our pieces list.
    # We use a unique name so we don't overwrite earlier entries.
    list_key <- paste0(cond_key, "_", tid)
    whole_transloc_pieces[[list_key]] <- summary_per_chrom
  }
}

# Stack all the pieces into one big data frame
whole_translocs <- bind_rows(whole_transloc_pieces)

# Drop duplicate (translocation x chromosome region) entries.
# A translocation might appear in both T1 and C1 BED files; we only need
# its region info once.
whole_translocs <- distinct(
  whole_translocs,
  transloc_name, chrom, start, end,
  .keep_all = TRUE
)

cat(sprintf("Found %d unique whole-translocation chromosome regions.\n",
            nrow(whole_translocs)))


# ============================================================================
# Step 4: Make the output folder
# ============================================================================

# showWarnings = FALSE means: don't warn if the folder already exists
# recursive = TRUE means: also create parent folders if needed
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# Step 5: Make per-segment alluvial plots
# ============================================================================
# One plot per row of all_segments

cat("\n--- Generating per-segment alluvial plots ---\n")

# seq_len(N) gives us 1, 2, 3, ..., N (R is 1-indexed!)
for (i in seq_len(nrow(all_segments))) {

  # Get this row's info
  seg <- all_segments[i, ]
  name <- seg$segment_name
  seg_chrom <- seg$chrom
  seg_start <- seg$start
  seg_end <- seg$end

  # ---------------------------------------------------------------------------
  # Decide whether this is a WT+C1-only plot
  # ---------------------------------------------------------------------------
  # Check if any of the WT_C1_ONLY names appears in this segment's name.
  # str_detect returns TRUE/FALSE for each value.
  is_wt_c1_only <- FALSE
  for (only_name in WT_C1_ONLY) {
    if (str_detect(name, fixed(only_name))) {
      is_wt_c1_only <- TRUE
    }
  }

  # Set the axes based on whether it's WT+C1 only
  if (is_wt_c1_only) {
    axes <- c("WT", "C1")
  } else {
    axes <- c("WT", "T1", "C1")
  }

  # ---------------------------------------------------------------------------
  # Get the bins inside this segment
  # ---------------------------------------------------------------------------
  sub_region <- sub[
    sub$chrom == seg_chrom &
    sub$start >= seg_start &
    sub$end   <= seg_end,
  ]

  # Skip if there are no bins
  if (nrow(sub_region) == 0) {
    cat(sprintf("  No bins for %s, skipping.\n", name))
    next   # next is like Python's continue, skip to the next iteration
  }

  cat(sprintf("  %s: %d bins (%s)\n",
              name,
              nrow(sub_region),
              paste(axes, collapse = " + ")))

  # ---------------------------------------------------------------------------
  # Count bins for each path through the conditions
  # ---------------------------------------------------------------------------
  # count() groups by the given columns and counts how many rows are in
  # each group. The count is in a column called "n" by default.
  if (is_wt_c1_only) {
    flow <- count(sub_region, MCF10A_WT.state, MCF10A_C1.state)
  } else {
    flow <- count(sub_region,
                  MCF10A_WT.state, MCF10A_T1.state, MCF10A_C1.state)
  }

  # ---------------------------------------------------------------------------
  # Build and save the plot
  # ---------------------------------------------------------------------------
  p <- build_alluvial_plot(flow, title = name, axes = axes)
  save_plot(p, name)
}


# ============================================================================
# Step 6: Make whole-translocation alluvial plots
# ============================================================================
# These combine all chromosomal regions of a translocation into one plot.
# Example: Der(17)t(3;17) has bits on chr3 AND chr17, the whole-translocation
# plot pools both into a single count.

cat("\n--- Generating whole-translocation alluvial plots ---\n")

# Get the list of unique translocation names
unique_names <- unique(whole_translocs$transloc_name)

for (trans_name in unique_names) {

  # Get all the chromosome regions for this translocation
  regions <- whole_translocs[whole_translocs$transloc_name == trans_name, ]

  # Decide if it's WT+C1 only
  is_wt_c1_only <- trans_name %in% WT_C1_ONLY

  if (is_wt_c1_only) {
    axes <- c("WT", "C1")
  } else {
    axes <- c("WT", "T1", "C1")
  }

  # ---------------------------------------------------------------------------
  # Collect bins from ALL chromosome regions for this translocation
  # ---------------------------------------------------------------------------
  # We loop through each region and find the matching bins, then stack
  # them into one big data frame.
  region_bins_pieces <- list()

  for (j in seq_len(nrow(regions))) {

    region_chrom <- regions$chrom[j]
    region_start <- regions$start[j]
    region_end <- regions$end[j]

    matching_bins <- sub[
      sub$chrom == region_chrom &
      sub$start >= region_start &
      sub$end   <= region_end,
    ]

    region_bins_pieces[[j]] <- matching_bins
  }

  # Stack them all together
  sub_all <- bind_rows(region_bins_pieces)

  # Drop duplicate bins (in case there's any overlap)
  sub_all <- distinct(sub_all, chrom, start, end, .keep_all = TRUE)

  # Skip if no bins were found
  if (nrow(sub_all) == 0) {
    cat(sprintf("  No bins for %s — skipping.\n", trans_name))
    next
  }

  cat(sprintf("  %s (whole): %d bins across %d chromosome region(s) (%s)\n",
              trans_name,
              nrow(sub_all),
              nrow(regions),
              paste(axes, collapse = " + ")))

  # Count bins by path
  if (is_wt_c1_only) {
    flow <- count(sub_all, MCF10A_WT.state, MCF10A_C1.state)
  } else {
    flow <- count(sub_all,
                  MCF10A_WT.state, MCF10A_T1.state, MCF10A_C1.state)
  }

  # Build the plot title
  plot_title <- paste0(trans_name, " (whole translocation)")

  # Build and save the plot
  p <- build_alluvial_plot(flow, title = plot_title, axes = axes)
  save_plot(p, plot_title)
}

# Done!
cat("\nDone! All plots saved to:", OUTPUT_DIR, "\n")
