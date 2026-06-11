# R/analysis/06b_plot_stability_vs_significance_poster.R
#
# Poster-ready combined stability vs significance plot for Pipeline A.
# Uses all available datasets from:
#   results/pipelineA/*/06_bootstrap_consensus/final_shortlist.csv
#
# Outputs:
#   results/analysis/figures/stability_vs_significance_pipelineA_poster.png
#
# Run from project root:
#   Rscript R/analysis/06b_plot_stability_vs_significance_poster.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
})

ROOT <- getwd()
IN_PATTERN <- file.path(ROOT, "results", "pipelineA", "*", "06_bootstrap_consensus", "final_shortlist.csv")
OUT_DIR <- file.path(ROOT, "results", "analysis", "figures")
OUT_PLOT <- file.path(OUT_DIR, "stability_vs_significance_pipelineA_poster.png")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

extract_dataset <- function(run_dir_name) {
  m <- str_match(run_dir_name, "(?i)DATASET?_([0-9]+)")
  if (!is.na(m[1, 2])) paste0("dataset_", m[1, 2]) else run_dir_name
}

files <- Sys.glob(IN_PATTERN)
if (length(files) == 0) {
  stop("No files found for pattern: ", IN_PATTERN)
}

plot_df <- map_dfr(files, function(fp) {
  run_dir <- basename(dirname(dirname(fp)))
  ds <- extract_dataset(run_dir)
  tbl <- read_csv(fp, show_col_types = FALSE)

  tbl %>%
    mutate(
      dataset = ds,
      category = case_when(
        selected == TRUE & consensus == TRUE  ~ "Selected + Consensus",
        selected == TRUE & consensus == FALSE ~ "Selected only",
        selected == FALSE                     ~ "Filtered out",
        TRUE                                  ~ "Filtered out"
      ),
      neg_log10_p_adjust = -log10(p.adjust)
    ) %>%
    select(dataset, ID, Description, collection, selected, consensus, stability, p.adjust, neg_log10_p_adjust, category)
})

plot_df <- plot_df %>%
  mutate(
    category = factor(
      category,
      levels = c("Filtered out", "Selected only", "Selected + Consensus")
    )
  )

# Split layers to reduce clutter and emphasize key categories
df_filtered <- plot_df %>% filter(category == "Filtered out")
df_sel_only <- plot_df %>% filter(category == "Selected only")
df_cons <- plot_df %>% filter(category == "Selected + Consensus")

p <- ggplot() +
  # Subtle background cloud: filtered-out pathways
  geom_point(
    data = df_filtered,
    aes(x = stability, y = neg_log10_p_adjust, color = category),
    size = 1.0, alpha = 0.30,
    position = position_jitter(width = 0.003, height = 0.03)
  ) +
  # Selected-only pathways
  geom_point(
    data = df_sel_only,
    aes(x = stability, y = neg_log10_p_adjust, color = category),
    size = 2.0, alpha = 0.75,
    position = position_jitter(width = 0.002, height = 0.02)
  ) +
  # Consensus pathways (dominant visual emphasis)
  geom_point(
    data = df_cons,
    aes(x = stability, y = neg_log10_p_adjust, color = category),
    size = 3.2, alpha = 0.96
  ) +
  # Stability threshold
  geom_vline(xintercept = 0.70, linetype = "dashed", linewidth = 0.7, color = "#4A4A4A") +
  # Significance threshold: p.adjust = 0.05
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.7, color = "#4A4A4A") +
  annotate(
    "text",
    x = 0.705, y = max(plot_df$neg_log10_p_adjust, na.rm = TRUE) * 0.98,
    label = "Stability threshold",
    hjust = 0, vjust = 1, size = 3.6, color = "#4A4A4A"
  ) +
  annotate(
    "text",
    x = 0.03, y = -log10(0.05) + 0.18,
    label = "Significance threshold",
    hjust = 0, vjust = 0, size = 3.6, color = "#4A4A4A", alpha = 0.95
  ) +
  annotate(
    "text",
    x = 0.05, y = 0.55,
    label = "Low-confidence pathways",
    hjust = 0, vjust = 0, size = 3.8, color = "#8A8A95"
  ) +
  annotate(
    "text",
    x = 1.045, y = max(plot_df$neg_log10_p_adjust, na.rm = TRUE) * 0.995,
    label = "High-confidence pathways",
    hjust = 1, vjust = 1, size = 3.8, color = "#3D2B7A"
  ) +
  scale_color_manual(values = c(
    "Filtered out" = "#C7C7D1",
    "Selected only" = "#3361BA",
    "Selected + Consensus" = "#6B40A6"
  )) +
  labs(
    title = "Stable and Significant Pathways Are Selected",
    x = "Stability",
    y = "Significance (\u2212log10 adjusted p-value)",
    color = "Pathway status"
  ) +
  coord_cartesian(xlim = c(0, 1.02)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 21, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E7E7EF", linewidth = 0.35),
    legend.spacing.y = grid::unit(6, "pt")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = c(1.0, 2.0, 3.2), alpha = c(0.30, 0.75, 0.96)),
      byrow = TRUE
    )
  )

ggsave(OUT_PLOT, p, width = 11.0, height = 7.6, dpi = 400)

cat("Wrote poster plot: ", OUT_PLOT, "\n", sep = "")
cat("Rows plotted: ", nrow(plot_df), "\n", sep = "")
