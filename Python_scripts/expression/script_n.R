# ===============================
# Load libraries
# ===============================
library(tidyverse)
library(ggalluvial)

# ===============================
# File paths
# ===============================
input_file <- "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/raw/subcompartments/GSE246947_MCF10A_WT_T1_C1_100000.subcompartments.bedGraph"
output_file <- "/Users/nadiaorning/Desktop/UiO/Høst2025/data/RNA-seq/folder/plots/subcompartment_alluvial_t317_2.pdf"

# ===============================
# Load data
# ===============================
sub <- read_tsv(input_file, show_col_types = FALSE)
sub_levels <- c("B3","B2","B1","B0","A0","A1","A2","A3")

# Keep only the columns we need
sub <- sub %>%
  mutate(MCF10A_WT.state = factor(MCF10A_WT.state, levels = sub_levels),
         MCF10A_T1.state = factor(MCF10A_T1.state, levels = sub_levels),
         MCF10A_C1.state = factor(MCF10A_C1.state, levels = sub_levels)) %>%
  drop_na()

# ===============================
# Filter to t(3;17) regions
# ===============================

sub_t317 <- sub %>%
  filter(
    (chrom == "chr3" & start >= 58600001 & end <= 198295559)|
      (chrom == "chr17" & start >= 22700001 & end <= 83257441)
  )

cat("Number of bins in t(3;17):", nrow(sub_t317), "\n")

# ===============================
# Count transitions
# ===============================
flow <- sub_t317 %>%
  count(MCF10A_WT.state,
        MCF10A_T1.state,
        MCF10A_C1.state)

# ===============================
# Define color palette for subcompartments
# ===============================
subcomp_colors <- c(
  "A3" = "#CB181D",
  "A2" = "#FB6A4A",
  "A1" = "#FCAE91",
  "A0" = "#FEE5D9",
  "B0" = "#DEEBF7",
  "B1" = "#9ECAE1",
  "B2" = "#4292C6",
  "B3" = "#084594"
)

# ===============================
# Plot alluvial (wide-format)
# ===============================
p <- ggplot(flow,
            aes(axis1 = MCF10A_WT.state,
                axis2 = MCF10A_T1.state,
                axis3 = MCF10A_C1.state,
                y = n)) +
  
  geom_alluvium(aes(fill = after_stat(stratum)), alpha = 0.7) +
  geom_stratum(aes(fill = after_stat(stratum)), color = "black", size = 0.2) +
  
  scale_fill_manual(values = subcomp_colors) +
  
  scale_x_discrete(limits = c("WT", "T1", "C1"),
                   labels = c("WT", "T1", "C1")) +
  
  scale_y_reverse(breaks = NULL) +
  
  labs(x = "Cell line",
       y = "Subcompartments") +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks.y = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ===============================
# Save high-resolution PDF
# ===============================
ggsave(output_file, p,
       width = 8,
       height = 6,
       dpi = 300)

cat("Alluvial plot for t(3;17) saved to PDF.\n")
