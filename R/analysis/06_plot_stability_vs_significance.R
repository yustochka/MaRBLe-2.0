# R/analysis/06_plot_stability_vs_significance.R
#
# Build a stability vs significance scatter plot from Pipeline A step06 outputs
# across all available datasets.
#
# Inputs:
#   results/pipelineA/*/06_bootstrap_consensus/final_shortlist.csv
#
# Outputs:
#   results/analysis/stability_significance_scatter_data.csv
#   results/analysis/figures/stability_vs_significance_pipelineA.png
#
# Run from project root:
#   Rscript R/analysis/06_plot_stability_vs_significance.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
})

ROOT <- getwd()
IN_PATTERN <- file.path(ROOT, "results", "pipelineA", "*", "06_bootstrap_consensus", "final_shortlist.csv")

OUT_DIR <- file.path(ROOT, "results", "analysis")
OUT_FIG_DIR <- file.path(OUT_DIR, "figures")
OUT_DATA <- file.path(OUT_DIR, "stability_significance_scatter_data.csv")
OUT_PLOT <- file.path(OUT_FIG_DIR, "stability_vs_significance_pipelineA.png")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

extract_dataset <- function(run_dir_name) {
  # Handles both DATASE_0_* and DATASET_0_* naming.
  m <- str_match(run_dir_name, "(?i)DATASET?_([0-9]+)")
  if (!is.na(m[1, 2])) {
    paste0("dataset_", m[1, 2])
  } else {
    run_dir_name
  }
}

files <- Sys.glob(IN_PATTERN)
if (length(files) == 0) {
  stop("No files found for pattern: ", IN_PATTERN)
}

plot_df <- map_dfr(files, function(fp) {
  run_dir <- basename(dirname(dirname(fp)))
  ds <- extract_dataset(run_dir)

  tbl <- read_csv(fp, show_col_types = FALSE)
  needed <- c("ID", "Description", "p.adjust", "stability", "selected", "consensus", "collection")
  missing <- setdiff(needed, names(tbl))
  if (length(missing) > 0) {
    stop("Missing columns in ", fp, ": ", paste(missing, collapse = ", "))
  }

  tbl %>%
    mutate(
      dataset = ds,
      run_dir = run_dir,
      category = case_when(
        selected == TRUE & consensus == TRUE  ~ "Selected + Consensus",
        selected == TRUE & consensus == FALSE ~ "Selected only",
        selected == FALSE                     ~ "Filtered out",
        TRUE                                  ~ "Filtered out"
      ),
      neg_log10_p_adjust = -log10(p.adjust)
    ) %>%
    select(dataset, run_dir, ID, Description, collection, selected, consensus, stability, p.adjust, neg_log10_p_adjust, category)
})

plot_df <- plot_df %>%
  mutate(
    dataset = factor(dataset, levels = sort(unique(dataset))),
    category = factor(category, levels = c("Selected + Consensus", "Selected only", "Filtered out"))
  )

write_csv(plot_df, OUT_DATA)

p <- ggplot(plot_df, aes(x = stability, y = neg_log10_p_adjust, color = category)) +
  geom_point(alpha = 0.8, size = 1.8) +
  geom_vline(xintercept = 0.70, linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_color_manual(values = c(
    "Selected + Consensus" = "#1b9e77",
    "Selected only" = "#d95f02",
    "Filtered out" = "#7570b3"
  )) +
  labs(
    title = "Pipeline A: Stability vs Significance",
    x = "Stability",
    y = "-log10(p.adjust)",
    color = "Category"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(OUT_PLOT, p, width = 11, height = 7, dpi = 300)

cat("Wrote data: ", OUT_DATA, "\n", sep = "")
cat("Wrote plot: ", OUT_PLOT, "\n", sep = "")
cat("Rows: ", nrow(plot_df), " | Datasets: ", paste(levels(plot_df$dataset), collapse = ", "), "\n", sep = "")
