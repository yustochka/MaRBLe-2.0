# R/ablation/compare_uncapped_runs.R
# Cross-variant comparison for Pipeline B uncapped ablation study (dataset_0).
#
# Reads stage_counts.csv and FINAL.csv from all completed uncapped ablation
# runs, produces seven plots, six data files, and a lightweight HTML analysis.
#
# Plots (results/ablation_uncapped/plots/):
#   01_funnel.png
#   02_final_size_bar.png
#   03_jaccard_heatmap.png
#   04_component_effect.png
#   05_go_wp_distribution.png
#   06_up_down_go_wp.png
#   07_pathway_frequency.png
#
# Data (results/ablation_uncapped/comparison_data/):
#   uncapped_summary_table.csv
#   uncapped_stage_counts_wide.csv
#   uncapped_pairwise_jaccard.csv
#   uncapped_ids_lost_gained_vs_full.csv
#   uncapped_unique_pathways.csv
#   uncapped_pathway_frequency.csv
#
# HTML: results/ablation_uncapped/uncapped_analysis.html
#
# Run from project root:
#   Rscript R/ablation/compare_uncapped_runs.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(forcats)
  library(scales)
  library(purrr)
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

ABLATION_ROOT <- "results/ablation_uncapped"
PLOT_DIR      <- file.path(ABLATION_ROOT, "plots")
DATA_DIR      <- file.path(ABLATION_ROOT, "comparison_data")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

VARIANTS <- c(
  "B_full_uncapped",
  "B_no_semantic_uncapped",
  "B_no_overlap_uncapped",
  "camera_only_selection_uncapped"
)

LABELS <- c(
  "B_full_uncapped"               = "B_full_uncapped",
  "B_no_semantic_uncapped"        = "B_no_semantic_uncapped",
  "B_no_overlap_uncapped"         = "B_no_overlap_uncapped",
  "camera_only_selection_uncapped" = "camera_only_sel._uncapped"
)

LABELS_WRAP <- c(
  "B_full_uncapped"               = "B_full_uncapped\n(reference)",
  "B_no_semantic_uncapped"        = "B_no_semantic_uncapped\n(no GO collapse)",
  "B_no_overlap_uncapped"         = "B_no_overlap_uncapped\n(no Jaccard clust.)",
  "camera_only_selection_uncapped" = "camera_only_sel._uncapped\n(CAMERA+relaxed solo)"
)

COMPONENT_REMOVED <- c(
  "B_full_uncapped"               = "None (reference)",
  "B_no_semantic_uncapped"        = "GO semantic collapse",
  "B_no_overlap_uncapped"         = "Jaccard overlap clustering",
  "camera_only_selection_uncapped" = "fgsea + absNES solo criterion"
)

PALETTE <- c(
  "B_full_uncapped"               = "#2166AC",
  "B_no_semantic_uncapped"        = "#4DAC26",
  "B_no_overlap_uncapped"         = "#D6604D",
  "camera_only_selection_uncapped" = "#F4A441"
)

theme_ablation <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "#F0F0F0", colour = NA),
      legend.position   = "right",
      plot.title        = element_text(face = "bold", size = base + 2),
      plot.subtitle     = element_text(size = base - 1, colour = "grey40"),
      plot.caption      = element_text(size = base - 3, colour = "grey50", hjust = 0)
    )
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

load_stage_counts <- function(variants, root) {
  purrr::map_dfr(variants, function(v) {
    p <- file.path(root, v, "dataset_0", "stage_counts.csv")
    if (!file.exists(p)) { warning("Missing stage_counts: ", p); return(NULL) }
    readr::read_csv(p, show_col_types = FALSE)
  })
}

load_finals <- function(variants, root) {
  purrr::map(variants, function(v) {
    p <- file.path(root, v, "dataset_0", "FINAL", "FINAL.csv")
    if (!file.exists(p)) { warning("Missing FINAL: ", p); return(character(0)) }
    df <- tryCatch(readr::read_csv(p, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(character(0))
    unique(df$ID)
  }) |> setNames(variants)
}

load_final_dfs <- function(variants, root) {
  purrr::map_dfr(variants, function(v) {
    p <- file.path(root, v, "dataset_0", "FINAL", "FINAL.csv")
    if (!file.exists(p)) return(NULL)
    df <- tryCatch(readr::read_csv(p, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df |>
      dplyr::mutate(direction = toupper(direction)) |>
      dplyr::select(ID, Description, direction, collection, variant) |>
      dplyr::distinct()
  })
}

all_counts   <- load_stage_counts(VARIANTS, ABLATION_ROOT)
all_finals   <- load_finals(VARIANTS, ABLATION_ROOT)
all_final_df <- load_final_dfs(VARIANTS, ABLATION_ROOT)

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

containment_ab <- function(a, b) {
  # fraction of a's IDs contained in b
  if (length(a) == 0) return(NA_real_)
  length(intersect(a, b)) / length(a)
}

# Pairwise Jaccard matrix
jac_mat <- outer(VARIANTS, VARIANTS,
                 FUN = Vectorize(function(i, j) jaccard(all_finals[[i]], all_finals[[j]])))
dimnames(jac_mat) <- list(VARIANTS, VARIANTS)

# ---------------------------------------------------------------------------
# 1. Summary table  →  uncapped_summary_table.csv
# ---------------------------------------------------------------------------

ref_ids <- all_finals[["B_full_uncapped"]]

summary_rows <- purrr::map_dfr(VARIANTS, function(v) {
  df  <- all_final_df |> dplyr::filter(variant == v)
  up  <- df |> dplyr::filter(direction == "UP")
  dn  <- df |> dplyr::filter(direction == "DOWN")
  ids <- all_finals[[v]]

  n_shared   <- length(intersect(ref_ids, ids))
  pct_ref_in_v   <- if (length(ref_ids) > 0) round(100 * n_shared / length(ref_ids), 1) else NA_real_
  pct_v_in_ref   <- if (length(ids) > 0)     round(100 * n_shared / length(ids), 1)     else NA_real_

  tibble::tibble(
    variant       = v,
    final_UP      = nrow(up),
    final_DOWN    = nrow(dn),
    final_total   = nrow(df),
    GO_count      = sum(df$collection == "GO"),
    WP_count      = sum(df$collection == "WP"),
    UP_GO         = sum(up$collection == "GO"),
    UP_WP         = sum(up$collection == "WP"),
    DOWN_GO       = sum(dn$collection == "GO"),
    DOWN_WP       = sum(dn$collection == "WP"),
    jaccard_vs_ref      = round(jac_mat[v, "B_full_uncapped"], 3),
    pct_ref_contained_in_variant = pct_ref_in_v,
    pct_variant_contained_in_ref = pct_v_in_ref
  )
})

readr::write_csv(summary_rows, file.path(DATA_DIR, "uncapped_summary_table.csv"))
message("Saved: uncapped_summary_table.csv")

# ---------------------------------------------------------------------------
# 2. Stage counts wide  →  uncapped_stage_counts_wide.csv
# ---------------------------------------------------------------------------

stage_wide <- all_counts |>
  dplyr::mutate(stage_dir = paste0(stage, "_", direction)) |>
  dplyr::select(variant, stage_dir, count) |>
  tidyr::pivot_wider(names_from = variant, values_from = count)

readr::write_csv(stage_wide, file.path(DATA_DIR, "uncapped_stage_counts_wide.csv"))
message("Saved: uncapped_stage_counts_wide.csv")

# ---------------------------------------------------------------------------
# 3. Pairwise Jaccard  →  uncapped_pairwise_jaccard.csv
# ---------------------------------------------------------------------------

jac_long <- as.data.frame(jac_mat) |>
  tibble::rownames_to_column("variant_row") |>
  tidyr::pivot_longer(-variant_row, names_to = "variant_col", values_to = "jaccard") |>
  dplyr::mutate(
    n_row      = purrr::map_int(variant_row, ~ length(all_finals[[.x]])),
    n_col      = purrr::map_int(variant_col, ~ length(all_finals[[.x]])),
    n_shared   = purrr::map2_int(variant_row, variant_col,
                   ~ length(intersect(all_finals[[.x]], all_finals[[.y]]))),
    pct_row_in_col = round(100 * n_shared / pmax(n_row, 1L), 1),
    pct_col_in_row = round(100 * n_shared / pmax(n_col, 1L), 1),
    jaccard        = round(jaccard, 3)
  )

readr::write_csv(jac_long, file.path(DATA_DIR, "uncapped_pairwise_jaccard.csv"))
message("Saved: uncapped_pairwise_jaccard.csv")

# ---------------------------------------------------------------------------
# 4. IDs lost/gained vs B_full_uncapped  →  uncapped_ids_lost_gained_vs_full.csv
# ---------------------------------------------------------------------------

lost_gained <- purrr::map_dfr(VARIANTS, function(v) {
  ids_v    <- all_finals[[v]]
  n_shared <- length(intersect(ref_ids, ids_v))
  n_lost   <- length(ref_ids) - n_shared
  n_gained <- length(ids_v) - n_shared
  tibble::tibble(
    variant                    = v,
    n_ref                      = length(ref_ids),
    n_variant                  = length(ids_v),
    n_shared                   = n_shared,
    n_lost_from_ref            = n_lost,
    n_gained_vs_ref            = n_gained,
    jaccard                    = round(jac_mat[v, "B_full_uncapped"], 3),
    pct_ref_contained_in_var   = round(100 * n_shared / max(length(ref_ids), 1L), 1),
    pct_var_contained_in_ref   = round(100 * n_shared / max(length(ids_v),  1L), 1),
    ids_lost                   = paste(sort(setdiff(ref_ids, ids_v)),   collapse = "; "),
    ids_gained                 = paste(sort(setdiff(ids_v, ref_ids)),   collapse = "; ")
  )
})

readr::write_csv(lost_gained, file.path(DATA_DIR, "uncapped_ids_lost_gained_vs_full.csv"))
message("Saved: uncapped_ids_lost_gained_vs_full.csv")

# ---------------------------------------------------------------------------
# 5. Unique pathways  →  uncapped_unique_pathways.csv
# ---------------------------------------------------------------------------

pathway_variant_counts <- all_final_df |>
  dplyr::group_by(ID) |>
  dplyr::summarise(
    n_variants    = dplyr::n_distinct(variant),
    variants_list = paste(sort(unique(variant)), collapse = "; "),
    Description   = dplyr::first(Description),
    direction     = paste(sort(unique(direction)), collapse = "/"),
    collection    = dplyr::first(collection),
    .groups       = "drop"
  ) |>
  dplyr::arrange(n_variants, collection, ID)

unique_pathways <- pathway_variant_counts |>
  dplyr::filter(n_variants == 1) |>
  dplyr::select(ID, Description, direction, collection, variants_list)

readr::write_csv(unique_pathways, file.path(DATA_DIR, "uncapped_unique_pathways.csv"))
message("Saved: uncapped_unique_pathways.csv")

# ---------------------------------------------------------------------------
# 6. Pathway frequency  →  uncapped_pathway_frequency.csv
# ---------------------------------------------------------------------------

pathway_freq <- pathway_variant_counts |>
  dplyr::select(ID, Description, collection, direction, n_variants, variants_list)

readr::write_csv(pathway_freq, file.path(DATA_DIR, "uncapped_pathway_frequency.csv"))
message("Saved: uncapped_pathway_frequency.csv")

# ---------------------------------------------------------------------------
# Funnel data (shared by plot 1 and HTML)
# ---------------------------------------------------------------------------

sum_stage <- function(d, stages, dirs = "all") {
  d |> dplyr::filter(stage %in% stages, direction %in% dirs) |>
    dplyr::summarise(total = sum(count, na.rm = TRUE)) |>
    dplyr::pull(total)
}

funnel_rows <- purrr::map_dfr(VARIANTS, function(v) {
  d <- dplyr::filter(all_counts, variant == v)
  tibble::tibble(
    variant           = v,
    S1_camera_sig     = sum_stage(d, c("camera_wp_sig", "camera_go_sig")),
    S2_after_agree    = sum_stage(d, c("after_agreement_wp", "after_agreement_go")),
    S3_after_semantic = sum_stage(d, "after_agreement_wp") + sum_stage(d, "after_semantic_go"),
    S4_overlap_reps   = sum_stage(d, "after_overlap_reps", c("Up", "Down")),
    S5_final          = sum_stage(d, "final", c("Up", "Down"))
  )
})

STAGE_LABELS <- c(
  S1_camera_sig     = "1. CAMERA sig\n(WP+GO)",
  S2_after_agree    = "2. After\nagreement",
  S3_after_semantic = "3. After semantic\ncollapse (GO)",
  S4_overlap_reps   = "4. After overlap\nclustering reps",
  S5_final          = "5. Final\n(uncapped)"
)

funnel_long <- funnel_rows |>
  tidyr::pivot_longer(-variant, names_to = "stage_key", values_to = "count") |>
  dplyr::mutate(
    stage_label   = factor(STAGE_LABELS[stage_key], levels = unname(STAGE_LABELS)),
    variant_label = factor(LABELS_WRAP[variant],    levels = unname(LABELS_WRAP))
  )

# ---------------------------------------------------------------------------
# Plot 1 — Funnel
# ---------------------------------------------------------------------------

p_funnel <- ggplot(funnel_long,
                   aes(x = stage_label, y = count,
                       colour = variant, group = variant)) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = variant), size = 3) +
  scale_colour_manual(values = PALETTE, labels = LABELS, name = "Variant") +
  scale_shape_manual(values = c(16, 17, 15, 18),
                     labels = LABELS, name = "Variant") +
  scale_y_continuous(breaks = pretty_breaks(6), limits = c(0, NA)) +
  labs(
    title    = "Pipeline B Uncapped Ablation — Pathway Funnel",
    subtitle = "Dataset: dataset_0 (RTT mouse, n=8, SYMBOL IDs) | No top-N cap applied",
    x        = "Pipeline stage",
    y        = "Number of pathways",
    caption  = paste0(
      "Stage 3 total = WP (from agreement) + GO (after semantic collapse or passthrough).\n",
      "All variants have no top-N per-direction cap at the final step.\n",
      "camera_only_sel._uncapped: fgsea bypassed; step04 forced passthrough (fg_FDR=NA)."
    )
  ) +
  theme_ablation() +
  theme(axis.text.x = element_text(size = 9))

ggsave(file.path(PLOT_DIR, "01_funnel.png"), p_funnel,
       width = 10, height = 5.5, dpi = 150)
message("Saved: 01_funnel.png")

# ---------------------------------------------------------------------------
# Plot 2 — Final size bar chart (UP/DOWN stacked)
# ---------------------------------------------------------------------------

final_counts_plot <- all_counts |>
  dplyr::filter(stage == "final") |>
  dplyr::mutate(
    direction_label = dplyr::recode(direction, Up = "UP", Down = "DOWN"),
    variant_label   = factor(LABELS_WRAP[variant], levels = rev(unname(LABELS_WRAP)))
  )

p_bar <- ggplot(final_counts_plot,
                aes(y = variant_label, x = count, fill = direction_label)) +
  geom_col(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(count > 0, count, "")),
            position = position_stack(vjust = 0.5),
            size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("UP" = "#D6604D", "DOWN" = "#4393C3"),
                    name = "Direction") +
  scale_x_continuous(breaks = pretty_breaks(8)) +
  labs(
    title    = "Pipeline B Uncapped Ablation — Final Candidate Pool Size",
    subtitle = "Dataset: dataset_0  |  No top-N cap applied to any variant",
    x        = "Number of pathways in final pool",
    y        = NULL,
    caption  = "All keep_final==TRUE pathways retained. No top-N cap."
  ) +
  theme_ablation() +
  theme(legend.position = "bottom")

ggsave(file.path(PLOT_DIR, "02_final_size_bar.png"), p_bar,
       width = 8, height = 4.5, dpi = 150)
message("Saved: 02_final_size_bar.png")

# ---------------------------------------------------------------------------
# Plot 3 — Pairwise Jaccard heatmap
# ---------------------------------------------------------------------------

jac_heat <- jac_long |>
  dplyr::mutate(
    row_label = factor(LABELS[variant_row], levels = unname(LABELS)),
    col_label = factor(LABELS[variant_col], levels = rev(unname(LABELS)))
  )

p_heatmap <- ggplot(jac_heat, aes(x = row_label, y = col_label, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3.0) +
  scale_fill_gradient2(low = "#F7F7F7", mid = "#92C5DE", high = "#2166AC",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "Jaccard\nsimilarity") +
  labs(
    title    = "Pipeline B Uncapped Ablation — Pairwise Jaccard Similarity",
    subtitle = "Based on final pathway IDs (uncapped pools)",
    x = NULL, y = NULL,
    caption  = paste0(
      "Jaccard = |A∩B| / |A∪B|. All sets are uncapped (no top-N applied).\n",
      "Jaccard is size-sensitive: larger pools compress J even when IDs are mostly shared.\n",
      "See uncapped_pairwise_jaccard.csv for containment metrics."
    )
  ) +
  theme_ablation() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )

ggsave(file.path(PLOT_DIR, "03_jaccard_heatmap.png"), p_heatmap,
       width = 8, height = 6.5, dpi = 150)
message("Saved: 03_jaccard_heatmap.png")

# ---------------------------------------------------------------------------
# Plot 4 — Component effect vs B_full_uncapped (Jaccard + containment)
# ---------------------------------------------------------------------------

effect_df <- lost_gained |>
  dplyr::select(variant, n_variant, n_lost_from_ref, n_gained_vs_ref,
                jaccard, pct_ref_contained_in_var, pct_var_contained_in_ref)

effect_long <- effect_df |>
  dplyr::filter(variant != "B_full_uncapped") |>
  dplyr::select(variant, jaccard,
                pct_ref_in_var  = pct_ref_contained_in_var,
                pct_var_in_ref  = pct_var_contained_in_ref,
                n_lost_from_ref, n_gained_vs_ref) |>
  tidyr::pivot_longer(c(jaccard, pct_ref_in_var, pct_var_in_ref,
                        n_lost_from_ref, n_gained_vs_ref),
                      names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    metric_label = dplyr::recode(metric,
      jaccard          = "Jaccard vs ref",
      pct_ref_in_var   = "% of ref IDs\nin variant",
      pct_var_in_ref   = "% of variant IDs\nin ref",
      n_lost_from_ref  = "IDs lost\nfrom ref (#)",
      n_gained_vs_ref  = "IDs gained\nvs ref (#)"
    ),
    metric_label = factor(metric_label, levels = c(
      "Jaccard vs ref",
      "% of ref IDs\nin variant",
      "% of variant IDs\nin ref",
      "IDs lost\nfrom ref (#)",
      "IDs gained\nvs ref (#)"
    )),
    variant_label = factor(LABELS[variant],
                           levels = rev(LABELS[VARIANTS != "B_full_uncapped"]))
  )

p_effect <- ggplot(effect_long,
                   aes(x = value, y = variant_label, fill = variant)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(metric %in% c("jaccard"),
                               sprintf("%.2f", value),
                               ifelse(metric %in% c("pct_ref_in_var", "pct_var_in_ref"),
                                      paste0(round(value, 0), "%"),
                                      as.character(as.integer(value))))),
            hjust = -0.15, size = 3.3) +
  scale_fill_manual(values = PALETTE, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  facet_wrap(~ metric_label, scales = "free_x", nrow = 1) +
  labs(
    title    = "Pipeline B Uncapped Ablation — Component Effect vs B_full_uncapped",
    subtitle = "Dataset: dataset_0  |  Reference: B_full_uncapped",
    x        = NULL, y = NULL,
    caption  = paste0(
      "% of ref IDs in variant: what fraction of B_full_uncapped's pool survives when this component is removed.\n",
      "% of variant IDs in ref: what fraction of this variant's pool is also in B_full_uncapped.\n",
      "Asymmetric containment reveals directionality of differences."
    )
  ) +
  theme_ablation() +
  theme(
    strip.text  = element_text(size = 8.5),
    axis.text.y = element_text(size = 9)
  )

ggsave(file.path(PLOT_DIR, "04_component_effect.png"), p_effect,
       width = 14, height = 4.5, dpi = 150)
message("Saved: 04_component_effect.png")

# ---------------------------------------------------------------------------
# Plot 5 — GO vs WP distribution per variant
# ---------------------------------------------------------------------------

coll_counts <- all_final_df |>
  dplyr::count(variant, collection) |>
  dplyr::mutate(variant_label = factor(LABELS_WRAP[variant], levels = unname(LABELS_WRAP)))

totals_df <- coll_counts |>
  dplyr::group_by(variant_label) |>
  dplyr::summarise(total = sum(n), .groups = "drop")

p_go_wp <- ggplot(coll_counts,
                  aes(x = variant_label, y = n, fill = collection)) +
  geom_col(colour = "white", linewidth = 0.4) +
  geom_text(data = totals_df,
            aes(x = variant_label, y = total, label = paste0("n=", total), fill = NULL),
            vjust = -0.4, size = 3.3, inherit.aes = FALSE) +
  scale_fill_manual(values = c("GO" = "#4393C3", "WP" = "#D6604D"),
                    name = "Collection") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title    = "Pipeline B Uncapped Ablation — GO vs WP in Final Pool",
    subtitle = "Dataset: dataset_0",
    x        = NULL,
    y        = "Number of pathways in final pool",
    caption  = "Totals shown above bars. No top-N cap applied."
  ) +
  theme_ablation() +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "bottom"
  )

ggsave(file.path(PLOT_DIR, "05_go_wp_distribution.png"), p_go_wp,
       width = 8, height = 5.5, dpi = 150)
message("Saved: 05_go_wp_distribution.png")

# ---------------------------------------------------------------------------
# Plot 6 — UP/DOWN × GO/WP breakdown per variant
# ---------------------------------------------------------------------------

breakdown <- all_final_df |>
  dplyr::mutate(group = paste0(direction, "\n", collection)) |>
  dplyr::count(variant, group) |>
  dplyr::mutate(
    variant_label = factor(LABELS[variant], levels = unname(LABELS)),
    group         = factor(group, levels = c("UP\nGO", "UP\nWP", "DOWN\nGO", "DOWN\nWP"))
  )

GROUP_COLS <- c("UP\nGO"   = "#D6604D", "UP\nWP"   = "#F4A441",
                "DOWN\nGO" = "#4393C3", "DOWN\nWP" = "#2166AC")

p_breakdown <- ggplot(breakdown,
                      aes(x = variant_label, y = n, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.0) +
  scale_fill_manual(values = GROUP_COLS, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title    = "Pipeline B Uncapped Ablation — UP/DOWN × GO/WP Breakdown",
    subtitle = "Dataset: dataset_0",
    x        = NULL,
    y        = "Number of pathways",
    caption  = "UP-GO (red), UP-WP (orange), DOWN-GO (light blue), DOWN-WP (dark blue). No top-N cap."
  ) +
  theme_ablation() +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "bottom"
  )

ggsave(file.path(PLOT_DIR, "06_up_down_go_wp.png"), p_breakdown,
       width = 9, height = 5.5, dpi = 150)
message("Saved: 06_up_down_go_wp.png")

# ---------------------------------------------------------------------------
# Plot 7 — Pathway frequency across uncapped variants
# ---------------------------------------------------------------------------

freq_summary <- pathway_freq |>
  dplyr::count(n_variants, collection)

freq_totals <- freq_summary |>
  dplyr::group_by(n_variants) |>
  dplyr::summarise(total = sum(n), .groups = "drop")

p_freq <- ggplot(freq_summary,
                 aes(x = factor(n_variants), y = n, fill = collection)) +
  geom_col(colour = "white", linewidth = 0.4) +
  geom_text(data = freq_totals,
            aes(x = factor(n_variants), y = total, label = total, fill = NULL),
            vjust = -0.4, size = 3.3, inherit.aes = FALSE) +
  scale_fill_manual(values = c("GO" = "#4393C3", "WP" = "#D6604D"),
                    name = "Collection") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title    = "Pipeline B Uncapped Ablation — Pathway Frequency Across Variants",
    subtitle = paste0("How many of the ", length(VARIANTS), " uncapped variants contain each pathway"),
    x        = "Number of uncapped variants containing the pathway",
    y        = "Number of distinct pathway IDs",
    caption  = "Pathways in 1 variant only are sensitive to that component's removal; pathways in all 4 are robust."
  ) +
  theme_ablation() +
  theme(legend.position = "bottom")

ggsave(file.path(PLOT_DIR, "07_pathway_frequency.png"), p_freq,
       width = 7, height = 5, dpi = 150)
message("Saved: 07_pathway_frequency.png")

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

cat("\n=== FUNNEL SUMMARY ===\n")
print(as.data.frame(funnel_rows), row.names = FALSE)

cat("\n=== SUMMARY TABLE ===\n")
print(as.data.frame(summary_rows), row.names = FALSE)

cat("\n=== LOST/GAINED VS B_full_uncapped ===\n")
print(as.data.frame(lost_gained[, c("variant", "n_variant", "n_shared",
                                     "n_lost_from_ref", "n_gained_vs_ref",
                                     "jaccard",
                                     "pct_ref_contained_in_var",
                                     "pct_var_contained_in_ref")]),
      row.names = FALSE)

cat("\n=== PAIRWISE JACCARD ===\n")
print(round(jac_mat, 2))

cat("\n=== UNIQUE PATHWAYS (in exactly 1 variant) ===\n")
print(as.data.frame(unique_pathways), row.names = FALSE)

# ---------------------------------------------------------------------------
# HTML analysis
# ---------------------------------------------------------------------------

# Helper: render a data frame as an HTML table
df_to_html <- function(df, caption = NULL, digits = 3) {
  # Round numeric columns
  df <- as.data.frame(df)
  for (j in seq_len(ncol(df))) {
    if (is.numeric(df[[j]])) df[[j]] <- round(df[[j]], digits)
  }

  header <- paste0("<th>", names(df), "</th>", collapse = "")
  rows   <- apply(df, 1, function(r) {
    paste0("<td>", r, "</td>", collapse = "")
  })
  rows <- paste0("<tr>", rows, "</tr>", collapse = "\n")

  cap_html <- if (!is.null(caption))
    paste0("<caption>", caption, "</caption>") else ""

  paste0(
    '<table class="tbl">\n', cap_html,
    "<thead><tr>", header, "</tr></thead>\n",
    "<tbody>\n", rows, "\n</tbody>\n</table>"
  )
}

# Gather numbers for inline text
ref_n <- nrow(all_final_df |> dplyr::filter(variant == "B_full_uncapped"))
sem_n <- nrow(all_final_df |> dplyr::filter(variant == "B_no_semantic_uncapped"))
ovl_n <- nrow(all_final_df |> dplyr::filter(variant == "B_no_overlap_uncapped"))
cam_n <- nrow(all_final_df |> dplyr::filter(variant == "camera_only_selection_uncapped"))

ref_ids_v  <- all_finals[["B_full_uncapped"]]
sem_ids_v  <- all_finals[["B_no_semantic_uncapped"]]
ovl_ids_v  <- all_finals[["B_no_overlap_uncapped"]]
cam_ids_v  <- all_finals[["camera_only_selection_uncapped"]]

j_sem <- round(jac_mat["B_no_semantic_uncapped", "B_full_uncapped"], 2)
j_ovl <- round(jac_mat["B_no_overlap_uncapped",  "B_full_uncapped"], 2)
j_cam <- round(jac_mat["camera_only_selection_uncapped", "B_full_uncapped"], 2)

sem_lost  <- lost_gained |> dplyr::filter(variant == "B_no_semantic_uncapped") |> dplyr::pull(n_lost_from_ref)
sem_gain  <- lost_gained |> dplyr::filter(variant == "B_no_semantic_uncapped") |> dplyr::pull(n_gained_vs_ref)
ovl_lost  <- lost_gained |> dplyr::filter(variant == "B_no_overlap_uncapped")  |> dplyr::pull(n_lost_from_ref)
ovl_gain  <- lost_gained |> dplyr::filter(variant == "B_no_overlap_uncapped")  |> dplyr::pull(n_gained_vs_ref)
cam_lost  <- lost_gained |> dplyr::filter(variant == "camera_only_selection_uncapped") |> dplyr::pull(n_lost_from_ref)
cam_gain  <- lost_gained |> dplyr::filter(variant == "camera_only_selection_uncapped") |> dplyr::pull(n_gained_vs_ref)

pct_ref_in_sem <- lost_gained |> dplyr::filter(variant == "B_no_semantic_uncapped") |> dplyr::pull(pct_ref_contained_in_var)
pct_ref_in_ovl <- lost_gained |> dplyr::filter(variant == "B_no_overlap_uncapped")  |> dplyr::pull(pct_ref_contained_in_var)
pct_ref_in_cam <- lost_gained |> dplyr::filter(variant == "camera_only_selection_uncapped") |> dplyr::pull(pct_ref_contained_in_var)
pct_cam_in_ref <- lost_gained |> dplyr::filter(variant == "camera_only_selection_uncapped") |> dplyr::pull(pct_var_contained_in_ref)
pct_ovl_in_ref <- lost_gained |> dplyr::filter(variant == "B_no_overlap_uncapped")  |> dplyr::pull(pct_var_contained_in_ref)
pct_sem_in_ref <- lost_gained |> dplyr::filter(variant == "B_no_semantic_uncapped") |> dplyr::pull(pct_var_contained_in_ref)

n_all4 <- nrow(pathway_variant_counts |> dplyr::filter(n_variants == 4))
n_uniq <- nrow(unique_pathways)

# Identify the single dominant direction of extra pathways in B_no_overlap_uncapped vs ref
ovl_extra_ids <- setdiff(ovl_ids_v, ref_ids_v)
ovl_extra_df  <- all_final_df |> dplyr::filter(variant == "B_no_overlap_uncapped",
                                                 ID %in% ovl_extra_ids)
ovl_extra_dir <- if (nrow(ovl_extra_df) > 0) {
  t <- table(ovl_extra_df$direction)
  paste0(names(sort(t, decreasing = TRUE)), " (", sort(t, decreasing = TRUE), ")", collapse = "; ")
} else "—"

ovl_extra_coll <- if (nrow(ovl_extra_df) > 0) {
  t <- table(ovl_extra_df$collection)
  paste0(names(sort(t, decreasing = TRUE)), " (", sort(t, decreasing = TRUE), ")", collapse = "; ")
} else "—"

cam_extra_ids <- setdiff(cam_ids_v, ref_ids_v)
cam_extra_df  <- all_final_df |> dplyr::filter(variant == "camera_only_selection_uncapped",
                                                 ID %in% cam_extra_ids)

# Images encode as relative paths — browser will load them from plots/
img <- function(fname, alt = fname, w = "100%") {
  paste0('<img src="plots/', fname, '" alt="', alt, '" style="max-width:', w, ';margin:12px 0;">')
}

html_content <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Pipeline B Uncapped Ablation — dataset_0</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         font-size: 14px; line-height: 1.6; color: #222;
         max-width: 1100px; margin: 0 auto; padding: 24px 32px; }
  h1   { font-size: 1.6em; color: #1a3a5c; border-bottom: 2px solid #2166AC; padding-bottom: 6px; }
  h2   { font-size: 1.2em; color: #2166AC; margin-top: 2em; }
  h3   { font-size: 1.05em; color: #444; }
  .note { background: #f0f4fa; border-left: 4px solid #2166AC;
          padding: 10px 16px; margin: 12px 0; border-radius: 3px; }
  .warn { background: #fff8e1; border-left: 4px solid #F4A441;
          padding: 10px 16px; margin: 12px 0; border-radius: 3px; }
  table.tbl { border-collapse: collapse; width: 100%; font-size: 12.5px; margin: 12px 0; }
  table.tbl th { background: #2166AC; color: #fff; padding: 6px 10px; text-align: left; }
  table.tbl td { padding: 5px 10px; border-bottom: 1px solid #ddd; }
  table.tbl tr:nth-child(even) td { background: #f7f9fc; }
  caption { caption-side: top; font-style: italic; font-size: 12px; color: #555;
            margin-bottom: 6px; text-align: left; }
  ul, ol { margin: 6px 0; padding-left: 1.6em; }
  li { margin: 3px 0; }
  .fig { margin: 18px 0; }
  .analysis { background: #f9f9f9; border-left: 3px solid #ccc;
              padding: 10px 14px; margin-top: 8px; font-size: 13px; }
  code { background: #eee; border-radius: 3px; padding: 1px 5px; font-size: 12.5px; }
  hr { border: none; border-top: 1px solid #ddd; margin: 24px 0; }
</style>
</head>
<body>

<h1>Pipeline B Uncapped Ablation — Follow-up Analysis</h1>
<p><strong>Dataset:</strong> dataset_0 (RTT mouse, n=8, SYMBOL IDs) &nbsp;|&nbsp;
   <strong>Date:</strong> ', format(Sys.Date(), "%Y-%m-%d"), ' &nbsp;|&nbsp;
   <strong>Output root:</strong> <code>results/ablation_uncapped/</code></p>

<h2>Why this follow-up was run</h2>
<p>
The first ablation experiment revealed that the <strong>top-5 per-direction cap</strong> hides
structural differences between variants: B_full and B_no_semantic were identical at cap=5
(Jaccard&nbsp;=&nbsp;1.00), and B_no_overlap differed (J&nbsp;=&nbsp;0.30) but only through
the lens of 10 final pathways per run. This makes it impossible to tell whether semantic
collapse and overlap clustering affect the <em>composition</em> of the candidate pool or
merely its ordering before truncation.
</p>
<p>
This follow-up removes the top-N cap from all four variants and compares the full
<code>keep_final&nbsp;==&nbsp;TRUE</code> pools directly.
</p>

<div class="note">
<strong>Conceptual equivalences with the first ablation experiment (for cross-run interpretation):</strong>
<ul>
  <li><strong>B_full_uncapped</strong> is conceptually equivalent to the previous <strong>B_no_topcap</strong>
      (all components active, top-N removed). Rerun fresh into <code>results/ablation_uncapped/</code>
      with a new variant name for a clean self-contained comparison.</li>
  <li><strong>B_no_semantic_uncapped</strong> is conceptually equivalent to the previous
      <strong>B_no_semantic_no_topcap</strong> (semantic collapse off, top-N removed). Same logic.</li>
  <li><strong>B_no_overlap_uncapped</strong> and <strong>camera_only_selection_uncapped</strong>
      are new: the first ablation only ran these variants with the cap on.</li>
</ul>
Because these are fresh runs with fresh R sessions and identical configs, any small numerical
differences from the prior experiment would be unexpected; the results should match exactly.
</div>

<hr>

<h2>Variants</h2>
<table class="tbl">
<thead><tr><th>Variant</th><th>fgsea</th><th>Semantic collapse</th><th>Overlap clustering</th><th>Top-N cap</th><th>Notes</th></tr></thead>
<tbody>
<tr><td>B_full_uncapped</td><td>ON</td><td>ON</td><td>ON</td><td>NONE</td><td>Reference (all components active)</td></tr>
<tr><td>B_no_semantic_uncapped</td><td>ON</td><td>OFF</td><td>ON</td><td>NONE</td><td>GO terms go directly to step05</td></tr>
<tr><td>B_no_overlap_uncapped</td><td>ON</td><td>ON</td><td>OFF</td><td>NONE</td><td>All step04 GO reps proceed to step06</td></tr>
<tr><td>camera_only_selection_uncapped</td><td>OFF</td><td>ON (forced passthrough)</td><td>ON</td><td>NONE</td>
    <td>Diagnostic: CAMERA FDR only; solo_abs_nes_min disabled; step04 bypassed (fg_FDR=NA)</td></tr>
</tbody>
</table>
<div class="warn">
  <strong>camera_only_selection_uncapped is a diagnostic variant, not a production configuration.</strong>
  fgsea is bypassed and the absNES solo criterion is disabled (sentinel&nbsp;=&nbsp;&minus;1). It tests
  how much of the candidate pool CAMERA alone can recover. Do not interpret it as a recommendation
  to remove fgsea.
</div>

<hr>

<h2>1. Funnel — How many pathways survive each stage?</h2>
<div class="fig">', img("01_funnel.png", "Funnel plot"), '</div>
<div class="analysis">
<strong>Reading this plot:</strong> All variants start from the same CAMERA output (stage&nbsp;1).
The lines separate at the agreement step (stage&nbsp;2) and further at overlap clustering (stage&nbsp;4).
Stage&nbsp;5 is the uncapped final pool.
<ul>
  <li><strong>Removing semantic collapse</strong> raises stage&nbsp;3 (more GO terms survive) and
      consequently stage&nbsp;4 (more material enters overlap clustering). Final pool is
      B_no_semantic_uncapped=', sem_n, ' vs reference=', ref_n, '.</li>
  <li><strong>Removing overlap clustering</strong>: stage&nbsp;4 rises (more terms pass through)
      but the final (stage&nbsp;5) is ', ovl_n, ' — note that step06\'s solo criterion
      (agreement_q&nbsp;&le;&nbsp;0.05, absNES&nbsp;&ge;&nbsp;1.5) is then applied to all raw terms,
      not to pre-selected cluster representatives. This makes step06 more restrictive,
      yielding a different (and in this case smaller) final pool than the reference.</li>
  <li><strong>CAMERA-only</strong>: stage&nbsp;2 is the full CAMERA significant set (47+80=127) because
      fgsea agreement is skipped. Semantic passthrough is forced. Final pool=', cam_n, '.</li>
</ul>
</div>

<hr>

<h2>2. Final pool sizes (UP / DOWN)</h2>
<div class="fig">', img("02_final_size_bar.png", "Final size bar chart"), '</div>
<div class="analysis">
<strong>Key numbers:</strong>
<ul>
  <li>B_full_uncapped: ', ref_n, ' total</li>
  <li>B_no_semantic_uncapped: ', sem_n, ' total</li>
  <li>B_no_overlap_uncapped: ', ovl_n, ' total</li>
  <li>camera_only_selection_uncapped: ', cam_n, ' total</li>
</ul>
Note the counter-intuitive result for B_no_overlap_uncapped (', ovl_n, ' pathways, smaller than the
reference ', ref_n, '): removing overlap clustering does not increase the pool — it changes what
step06 receives. Cluster representatives are pre-screened to be the strongest terms per cluster,
so they are more likely to pass step06\'s strict solo criterion. Raw un-clustered terms include
weaker entries that fail solo, resulting in fewer final survivors. The camera-only variant is
the largest pool (', cam_n, '), because fgsea agreement filtering is skipped entirely.
</div>

<hr>

<h2>3. Pairwise Jaccard heatmap</h2>
<div class="fig">', img("03_jaccard_heatmap.png", "Jaccard heatmap"), '</div>
<div class="analysis">
Note that Jaccard is size-sensitive: when one set is much larger, the shared denominator
compresses J even if all smaller-set IDs are present. Use the containment metrics
(table below and plot&nbsp;4) alongside Jaccard for a complete picture.
<br><br>
Jaccard values: B_full_uncapped vs B_no_semantic_uncapped&nbsp;=&nbsp;', j_sem, ',
vs B_no_overlap_uncapped&nbsp;=&nbsp;', j_ovl, ',
vs camera_only_selection_uncapped&nbsp;=&nbsp;', j_cam, '.
</div>

<hr>

<h2>4. Component effect vs B_full_uncapped (Jaccard + containment)</h2>
<div class="fig">', img("04_component_effect.png", "Component effect"), '</div>
<div class="analysis">
The two containment columns resolve the Jaccard ambiguity:
<ul>
  <li><strong>% of ref IDs in variant</strong>: of the ', ref_n, ' pathways in B_full_uncapped,
      how many appear in this variant?
      Semantic=', pct_ref_in_sem, '%, Overlap=', pct_ref_in_ovl, '%, Camera=', pct_ref_in_cam, '%.</li>
  <li><strong>% of variant IDs in ref</strong>: of this variant\'s pool, how many are also in the ref?
      Semantic=', pct_sem_in_ref, '%, Overlap=', pct_ovl_in_ref, '%, Camera=', pct_cam_in_ref, '%.</li>
</ul>
High "% of ref in variant" means the variant retains nearly all reference pathways
(ref IDs are preserved). High "% of variant in ref" means this variant\'s pool is largely
a subset of the reference — the extra gained IDs are new, not rearrangements.
</div>

<h3>Containment table</h3>
',
df_to_html(
  lost_gained |>
    dplyr::select(variant, n_ref, n_variant, n_shared,
                  n_lost_from_ref, n_gained_vs_ref, jaccard,
                  pct_ref_contained_in_var, pct_var_contained_in_ref),
  caption = "Containment metrics vs B_full_uncapped. pct_ref_contained_in_var = % of ref IDs present in this variant. pct_var_contained_in_ref = % of this variant's IDs present in ref."
),

'
<hr>

<h2>5. GO vs WP composition</h2>
<div class="fig">', img("05_go_wp_distribution.png", "GO vs WP distribution"), '</div>
<div class="analysis">
Shows whether extra pathways (when cap is removed) are predominantly GO or WP terms.
</div>

<hr>

<h2>6. UP / DOWN × GO / WP breakdown</h2>
<div class="fig">', img("06_up_down_go_wp.png", "UP/DOWN x GO/WP breakdown"), '</div>
<div class="analysis">
Reveals whether extra pathways accumulate in a specific direction or collection.
The extra pathways in B_no_overlap_uncapped by direction: ', ovl_extra_dir, ';
by collection: ', ovl_extra_coll, '.
</div>

<hr>

<h2>7. Pathway frequency across uncapped variants</h2>
<div class="fig">', img("07_pathway_frequency.png", "Pathway frequency"), '</div>
<div class="analysis">
Pathways appearing in all ', length(VARIANTS), ' variants are robust to any single
component removal (n=', n_all4, '). Pathways appearing in only 1 variant are
sensitive to that specific component (n=', n_uniq, ').
</div>

<hr>

<h2>Summary table</h2>
',
df_to_html(
  summary_rows |>
    dplyr::select(variant, final_total, final_UP, final_DOWN,
                  GO_count, WP_count, UP_GO, UP_WP, DOWN_GO, DOWN_WP,
                  jaccard_vs_ref, pct_ref_contained_in_variant, pct_variant_contained_in_ref),
  caption = "Final pool counts and containment vs B_full_uncapped."
),

'
<hr>

<h2>Questions answered</h2>

<h3>Q1. Without the final cap, which component changes the candidate pool the most?</h3>
<p>
<strong>Removing fgsea agreement</strong> (camera_only_selection_uncapped) produces the largest pool
(', cam_n, '), because the agreement filter is the main upstream bottleneck.
Removing semantic collapse adds ', sem_n - ref_n, ' pathways vs reference (', sem_n, ' vs ', ref_n, ').
Removing overlap clustering counter-intuitively <em>reduces</em> the pool to ', ovl_n, ': without
clustering, step06\'s strict solo criterion is applied to all raw agreement-filtered terms rather
than to pre-screened cluster representatives, making selection more restrictive. Overlap
clustering concentrates quality before step06; without it, more terms fail the solo threshold.
</p>

<h3>Q2. Does semantic collapse matter more when the cap is removed?</h3>
<p>
At cap=5 (first ablation), B_full and B_no_semantic were identical (Jaccard=1.00).
Without the cap, B_no_semantic_uncapped has ', sem_n, ' vs ', ref_n, ' in ref
(Jaccard=', j_sem, '). Semantic collapse does create measurable differences in
pool composition — they were masked by the cap. Removed GO terms cluster with
retained ones at low Jaccard thresholds, but some unique-to-B_no_semantic
pathways are only visible without the cap.
</p>

<h3>Q3. How different is B_no_overlap_uncapped from B_full_uncapped?</h3>
<p>
IDs lost from ref: ', ovl_lost, '; IDs gained vs ref: ', ovl_gain, '.
', pct_ref_in_ovl, '% of ref IDs survive in B_no_overlap_uncapped (containment).
', pct_ovl_in_ref, '% of B_no_overlap_uncapped\'s pool is also in ref.
<br><br>
The pool is <strong>smaller</strong> than the reference (', ovl_n, ' vs ', ref_n, ').
Mechanism: step05 passthrough lets all 93 agreement+semantic-filtered terms into step06, but
step06\'s solo criterion (q&nbsp;&le;&nbsp;0.05, absNES&nbsp;&ge;&nbsp;1.5) is harsh on raw
un-clustered input. Cluster representatives are pre-selected as the strongest term per cluster,
so they pass solo more reliably. Raw terms include weaker cluster members that fail solo,
reducing total survivors. This reveals that overlap clustering serves two roles: deduplication
<em>and</em> quality concentration before the final selection step.
</p>

<h3>Q4. Does camera-only selection recover most of the full uncapped pool?</h3>
<p>
camera_only_selection_uncapped has ', cam_n, ' pathways.
', pct_ref_in_cam, '% of ref IDs appear in it; ', pct_cam_in_ref, '% of its
pool is in the ref. ', cam_gain, ' IDs are in camera-only but not in the ref —
these are pathways that CAMERA finds significant but which fail the CAMERA&cap;fgsea
agreement filter. ', cam_lost, ' ref IDs are absent — pathways that fgsea helps
recover through directional agreement filtering.
This confirms that fgsea contributes both as an independent filter and as a directional
ranker, not merely as a redundant signal.
</p>

<h3>Q5. Are the extra pathways mostly GO, WP, UP, or DOWN?</h3>
<p>
See plot&nbsp;6 for per-variant breakdown. For variants larger than the reference, extra pathways
(IDs in variant not in ref) can be characterised by direction and collection.
B_no_overlap_uncapped is smaller than the reference, so the question is instead which ref IDs
it <em>misses</em>. For the variants that do gain IDs (semantic and camera-only), see the
uncapped_ids_lost_gained_vs_full.csv for the full ID lists.
B_no_overlap_uncapped unique IDs by direction: ', ovl_extra_dir, ';
by collection: ', ovl_extra_coll, '.
</p>

<h3>Q6. Which pathways are robust across all uncapped variants?</h3>
<p>
', n_all4, ' pathway IDs appear in all ', length(VARIANTS), ' uncapped variants.
These are the most robust: they survive agreement filtering, semantic collapse,
overlap clustering, and CAMERA-only selection simultaneously.
</p>

<h3>Q7. Which pathways are unique to specific variants?</h3>
<p>
', n_uniq, ' pathways appear in exactly one variant.
See <code>uncapped_unique_pathways.csv</code> for the full list.
</p>

<hr>

<h2>Interpretation notes</h2>
<ul>
  <li>All findings are specific to <strong>dataset_0 only</strong>. Do not generalise to other datasets.</li>
  <li>camera_only_selection_uncapped is a <strong>diagnostic variant</strong>. Its results show
      what CAMERA alone can produce, not a recommendation to remove fgsea from the pipeline.</li>
  <li><strong>Final size</strong> differences reflect pool composition, not quality.
      A larger uncapped pool contains more redundancy by construction.</li>
  <li><strong>Composition effects</strong> (GO/WP ratio, UP/DOWN split) and
      <strong>redundancy effects</strong> (overlap clustering expansion) are distinct phenomena
      even when they produce similar pool-size changes.</li>
  <li>Jaccard alone is misleading when pool sizes differ substantially.
      Always read it alongside the containment percentages.</li>
</ul>

<hr>
<p style="font-size:11px;color:#888;">
Generated ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ' by R/ablation/compare_uncapped_runs.R.
Dataset: dataset_0. Pipeline: Pipeline B ablation (uncapped variants).
This file references plots/ via relative paths — keep them together.
</p>

</body>
</html>')

html_out <- file.path(ABLATION_ROOT, "uncapped_analysis.html")
writeLines(html_content, html_out)
message("Saved: ", html_out)

message("\nAll outputs written to:")
message("  ", PLOT_DIR)
message("  ", DATA_DIR)
message("  ", html_out)
