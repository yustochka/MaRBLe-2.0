# R/ablation/compare_ablation_runs.R
# Cross-variant comparison for Pipeline B ablation study (dataset_0).
#
# Reads stage_counts.csv and FINAL.csv from all completed ablation runs,
# produces seven plots and three data files:
#
# Plots (results/ablation/plots/):
#   01_funnel.png              — pathway count at each pipeline stage
#   02_final_size_bar.png      — final shortlist size per variant (UP/DOWN)
#   03_jaccard_heatmap.png     — pairwise Jaccard similarity heatmap
#   04_component_effect.png    — component effect summary vs B_full
#   05_go_wp_distribution.png  — GO vs WP composition per variant
#   06_up_down_go_wp.png       — UP/DOWN × GO/WP breakdown per variant
#   07_pathway_frequency.png   — how many variants contain each pathway
#
# Data (results/ablation/comparison_data/):
#   funnel_data.csv
#   pairwise_jaccard.csv
#   all_pathways_variant_counts.csv
#   unique_pathways.csv
#   main_findings_table.csv
#
# Run from project root:
#   Rscript R/ablation/compare_ablation_runs.R

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

ABLATION_ROOT <- "results/ablation"
PLOT_DIR      <- file.path(ABLATION_ROOT, "plots")
DATA_DIR      <- file.path(ABLATION_ROOT, "comparison_data")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# Canonical variant order and display labels
VARIANTS <- c(
  "B_full",
  "B_no_semantic",
  "B_no_overlap",
  "B_no_topcap",
  "B_no_fgsea",
  "B_no_fgsea_camera_only_selection",
  "B_no_semantic_no_topcap"
)

LABELS <- c(
  "B_full"                             = "B_full (baseline)",
  "B_no_semantic"                      = "B_no_semantic",
  "B_no_overlap"                       = "B_no_overlap",
  "B_no_topcap"                        = "B_no_topcap",
  "B_no_fgsea"                         = "B_no_fgsea",
  "B_no_fgsea_camera_only_selection"   = "camera_only_sel.",
  "B_no_semantic_no_topcap"            = "B_no_sem_no_topcap"
)

LABELS_WRAP <- c(
  "B_full"                             = "B_full\n(baseline)",
  "B_no_semantic"                      = "B_no_semantic\n(no GO collapse)",
  "B_no_overlap"                       = "B_no_overlap\n(no Jaccard clust.)",
  "B_no_topcap"                        = "B_no_topcap\n(no top-5 cap)",
  "B_no_fgsea"                         = "B_no_fgsea\n(CAMERA only, empty)",
  "B_no_fgsea_camera_only_selection"   = "camera_only_sel.\n(CAMERA+relaxed sel.)",
  "B_no_semantic_no_topcap"            = "B_no_sem_no_topcap\n(no collapse, no cap)"
)

COMPONENT_REMOVED <- c(
  "B_full"                             = "None (baseline)",
  "B_no_semantic"                      = "GO semantic collapse",
  "B_no_overlap"                       = "Jaccard overlap clustering",
  "B_no_topcap"                        = "Top-5 per direction cap",
  "B_no_fgsea"                         = "fgsea (CAMERA only, selection fails)",
  "B_no_fgsea_camera_only_selection"   = "fgsea + absNES solo criterion",
  "B_no_semantic_no_topcap"            = "GO semantic collapse + top-5 cap"
)

KEY_INTERP <- c(
  "B_full"                             = "Reference; all components active",
  "B_no_semantic"                      = "No effect on final shortlist (capped); absorbed by Jaccard step",
  "B_no_overlap"                       = "7/10 IDs change; dominant driver of pathway identity",
  "B_no_topcap"                        = "Exposes 45 hidden valid candidates; identity order unchanged",
  "B_no_fgsea"                         = "Selection fails (absNES=NA); 0 output is mechanical not biological",
  "B_no_fgsea_camera_only_selection"   = "CAMERA alone recovers 9/10 output; fgsea acts as ranking gate",
  "B_no_semantic_no_topcap"            = "8 additional GO terms visible; Jaccard=0.84 vs B_no_topcap"
)

PALETTE <- c(
  "B_full"                             = "#2166AC",
  "B_no_semantic"                      = "#4DAC26",
  "B_no_overlap"                       = "#D6604D",
  "B_no_topcap"                        = "#8E44AD",
  "B_no_fgsea"                         = "#999999",
  "B_no_fgsea_camera_only_selection"   = "#F4A441",
  "B_no_semantic_no_topcap"            = "#1A7A5E"
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

# Pairwise Jaccard matrix
jac_mat <- outer(VARIANTS, VARIANTS,
                 FUN = Vectorize(function(i, j) jaccard(all_finals[[i]], all_finals[[j]])))
dimnames(jac_mat) <- list(VARIANTS, VARIANTS)

# ---------------------------------------------------------------------------
# 1. Funnel data
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
    S1_camera_sig     = sum_stage(d, c("camera_wp_sig","camera_go_sig")),
    S2_after_agree    = sum_stage(d, c("after_agreement_wp","after_agreement_go")),
    S3_after_semantic = sum_stage(d, "after_agreement_wp") + sum_stage(d, "after_semantic_go"),
    S4_overlap_reps   = sum_stage(d, "after_overlap_reps", c("Up","Down")),
    S5_final          = sum_stage(d, "final", c("Up","Down"))
  )
})

STAGE_LABELS <- c(
  S1_camera_sig     = "1. CAMERA sig\n(WP+GO)",
  S2_after_agree    = "2. After\nagreement",
  S3_after_semantic = "3. After semantic\ncollapse (GO)",
  S4_overlap_reps   = "4. After overlap\nclustering reps",
  S5_final          = "5. Final\nshortlist"
)

funnel_long <- funnel_rows |>
  tidyr::pivot_longer(-variant, names_to = "stage_key", values_to = "count") |>
  dplyr::mutate(
    stage_label   = factor(STAGE_LABELS[stage_key], levels = unname(STAGE_LABELS)),
    variant_label = factor(LABELS_WRAP[variant],    levels = unname(LABELS_WRAP))
  )

readr::write_csv(
  funnel_long |> dplyr::select(variant, stage_key, stage_label, count),
  file.path(DATA_DIR, "funnel_data.csv")
)

p_funnel <- ggplot(funnel_long,
                   aes(x = stage_label, y = count,
                       colour = variant, group = variant)) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = variant), size = 3) +
  scale_colour_manual(values = PALETTE, labels = LABELS, name = "Variant") +
  scale_shape_manual(values = c(16,17,15,18,4,8,3),
                     labels = LABELS, name = "Variant") +
  scale_y_continuous(breaks = pretty_breaks(6), limits = c(0, NA)) +
  labs(
    title    = "Pipeline B Ablation — Pathway Funnel",
    subtitle = "Dataset: dataset_0 (RTT mouse, n=8, SYMBOL IDs)",
    x        = "Pipeline stage",
    y        = "Number of pathways",
    caption  = paste0(
      "Stage 3 total = WP (from agreement) + GO (after semantic collapse or passthrough).\n",
      "B_no_fgsea final = 0: solo_high_conf cannot activate without absNES.\n",
      "Uncapped variants (B_no_topcap, B_no_sem_no_topcap) retain all step06 candidates."
    )
  ) +
  theme_ablation() +
  theme(axis.text.x = element_text(size = 9))

ggsave(file.path(PLOT_DIR, "01_funnel.png"), p_funnel,
       width = 10, height = 5.5, dpi = 150)
message("Saved: 01_funnel.png")

# ---------------------------------------------------------------------------
# 2. Final size bar chart (UP/DOWN stacked)
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
    title    = "Pipeline B Ablation — Final Shortlist Size",
    subtitle = "Dataset: dataset_0  |  standard cap = 5 per direction",
    x        = "Number of pathways in final shortlist",
    y        = NULL,
    caption  = paste0(
      "B_no_fgsea = 0 (solo_high_conf requires absNES, which is NA when fgsea is bypassed).\n",
      "B_no_topcap and B_no_sem_no_topcap have no per-direction cap applied."
    )
  ) +
  theme_ablation() +
  theme(legend.position = "bottom")

ggsave(file.path(PLOT_DIR, "02_final_size_bar.png"), p_bar,
       width = 8, height = 5.5, dpi = 150)
message("Saved: 02_final_size_bar.png")

# ---------------------------------------------------------------------------
# 3. Pairwise Jaccard heatmap
# ---------------------------------------------------------------------------

jac_long <- as.data.frame(jac_mat) |>
  tibble::rownames_to_column("variant_row") |>
  tidyr::pivot_longer(-variant_row, names_to = "variant_col", values_to = "jaccard") |>
  dplyr::mutate(
    row_label = factor(LABELS[variant_row], levels = unname(LABELS)),
    col_label = factor(LABELS[variant_col], levels = rev(unname(LABELS)))
  )

readr::write_csv(
  jac_long |> dplyr::select(variant_row, variant_col, jaccard),
  file.path(DATA_DIR, "pairwise_jaccard.csv")
)

p_heatmap <- ggplot(jac_long, aes(x = row_label, y = col_label, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3.0) +
  scale_fill_gradient2(low = "#F7F7F7", mid = "#92C5DE", high = "#2166AC",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "Jaccard\nsimilarity") +
  labs(
    title    = "Pipeline B Ablation — Pairwise Jaccard Similarity",
    subtitle = "Based on final shortlist pathway IDs",
    x = NULL, y = NULL,
    caption  = paste0(
      "Jaccard = |A∩B| / |A∪B|. B_no_fgsea has empty set → J = 0 with all others.\n",
      "Uncapped variants have larger denominators, reducing J even when they contain all smaller-set IDs.\n",
      "B_full (J=1.00 with B_no_semantic) confirms no effect from semantic collapse on the capped shortlist."
    )
  ) +
  theme_ablation() +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 8)
  )

ggsave(file.path(PLOT_DIR, "03_jaccard_heatmap.png"), p_heatmap,
       width = 9, height = 7.5, dpi = 150)
message("Saved: 03_jaccard_heatmap.png")

# ---------------------------------------------------------------------------
# 4. Component effect summary chart (vs B_full)
# ---------------------------------------------------------------------------

b_full_ids <- all_finals[["B_full"]]
n_b_full   <- length(b_full_ids)

effect_df <- purrr::map_dfr(VARIANTS, function(v) {
  ids_v     <- all_finals[[v]]
  n_v       <- length(ids_v)
  n_shared  <- length(intersect(b_full_ids, ids_v))
  n_lost    <- n_b_full - n_shared        # in B_full, not in variant
  n_gained  <- n_v - n_shared             # in variant, not in B_full
  tibble::tibble(
    variant      = v,
    final_size   = n_v,
    size_diff    = n_v - n_b_full,
    jaccard_full = jac_mat[v, "B_full"],
    ids_lost     = as.integer(n_lost),
    ids_gained   = as.integer(n_gained)
  )
})

# Long form for faceted plot
effect_long <- effect_df |>
  dplyr::filter(variant != "B_full") |>
  dplyr::select(variant, jaccard_full, ids_lost, ids_gained) |>
  tidyr::pivot_longer(c(jaccard_full, ids_lost, ids_gained),
                      names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    metric_label = dplyr::recode(metric,
      jaccard_full = "Jaccard similarity\nvs B_full (0–1)",
      ids_lost     = "Pathways from B_full\nnot in variant (#)",
      ids_gained   = "New pathways not\nin B_full (#)"
    ),
    metric_label = factor(metric_label, levels = c(
      "Jaccard similarity\nvs B_full (0–1)",
      "Pathways from B_full\nnot in variant (#)",
      "New pathways not\nin B_full (#)"
    )),
    variant_label = factor(LABELS[variant], levels = rev(LABELS[VARIANTS != "B_full"]))
  )

p_effect <- ggplot(effect_long,
                   aes(x = value, y = variant_label, fill = variant)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(metric == "jaccard_full",
                               sprintf("%.2f", value),
                               as.character(as.integer(value)))),
            hjust = -0.15, size = 3.3) +
  scale_fill_manual(values = PALETTE, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  facet_wrap(~ metric_label, scales = "free_x", nrow = 1) +
  labs(
    title    = "Pipeline B Ablation — Component Effect vs B_full",
    subtitle = "Dataset: dataset_0",
    x        = NULL, y = NULL,
    caption  = paste0(
      "B_no_fgsea: ids_lost=10, ids_gained=0 (empty output).\n",
      "Uncapped variants (B_no_topcap, B_no_sem_no_topcap): ids_lost=0; large ids_gained due to cap removal.\n",
      "Jaccard is asymmetric: uncapped sets have larger denominators, compressing J even when all B_full IDs are present."
    )
  ) +
  theme_ablation() +
  theme(
    strip.text  = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )

ggsave(file.path(PLOT_DIR, "04_component_effect.png"), p_effect,
       width = 12, height = 5, dpi = 150)
message("Saved: 04_component_effect.png")

# ---------------------------------------------------------------------------
# 5. GO vs WP distribution per variant
# ---------------------------------------------------------------------------

coll_counts <- all_final_df |>
  dplyr::count(variant, collection) |>
  dplyr::mutate(
    variant_label = factor(LABELS_WRAP[variant], levels = unname(LABELS_WRAP)),
    total_v = sum(n), .by = variant
  )

# Total labels per variant
totals_df <- coll_counts |>
  dplyr::group_by(variant_label) |>
  dplyr::summarise(total = sum(n), .groups = "drop")

p_go_wp <- ggplot(coll_counts,
                  aes(x = variant_label, y = n, fill = collection)) +
  geom_col(colour = "white", linewidth = 0.4) +
  geom_text(data = totals_df,
            aes(x = variant_label, y = total, label = paste0("n=", total),
                fill = NULL),
            vjust = -0.4, size = 3.3, inherit.aes = FALSE) +
  scale_fill_manual(values = c("GO" = "#4393C3", "WP" = "#D6604D"),
                    name = "Collection") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Pipeline B Ablation — GO vs WP in Final Shortlist",
    subtitle = "Dataset: dataset_0",
    x        = NULL,
    y        = "Number of pathways in FINAL.csv",
    caption  = "B_no_fgsea excluded (empty output). Totals shown above bars."
  ) +
  theme_ablation() +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "bottom"
  )

ggsave(file.path(PLOT_DIR, "05_go_wp_distribution.png"), p_go_wp,
       width = 9, height = 5.5, dpi = 150)
message("Saved: 05_go_wp_distribution.png")

# ---------------------------------------------------------------------------
# 6. UP/DOWN × GO/WP breakdown per variant
# ---------------------------------------------------------------------------

breakdown <- all_final_df |>
  dplyr::mutate(group = paste0(direction, "\n", collection)) |>
  dplyr::count(variant, group) |>
  dplyr::mutate(
    variant_label = factor(LABELS[variant], levels = unname(LABELS)),
    group         = factor(group, levels = c("UP\nGO","UP\nWP","DOWN\nGO","DOWN\nWP"))
  )

GROUP_COLS <- c("UP\nGO" = "#D6604D", "UP\nWP" = "#F4A441",
                "DOWN\nGO" = "#4393C3", "DOWN\nWP" = "#2166AC")

p_breakdown <- ggplot(breakdown,
                      aes(x = variant_label, y = n, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.0) +
  scale_fill_manual(values = GROUP_COLS, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Pipeline B Ablation — UP/DOWN × GO/WP Breakdown",
    subtitle = "Dataset: dataset_0",
    x        = NULL,
    y        = "Number of pathways",
    caption  = "B_no_fgsea excluded (empty). Grouped bars per variant: UP-GO (red), UP-WP (orange), DOWN-GO (light blue), DOWN-WP (dark blue)."
  ) +
  theme_ablation() +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "bottom"
  )

ggsave(file.path(PLOT_DIR, "06_up_down_go_wp.png"), p_breakdown,
       width = 10, height = 5.5, dpi = 150)
message("Saved: 06_up_down_go_wp.png")

# ---------------------------------------------------------------------------
# 7. Shared pathway frequency plot
# ---------------------------------------------------------------------------

pathway_freq <- all_final_df |>
  dplyr::group_by(ID, collection) |>
  dplyr::summarise(n_variants = dplyr::n_distinct(variant), .groups = "drop")

freq_summary <- pathway_freq |>
  dplyr::count(n_variants, collection)

# Helper for total per bar
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Pipeline B Ablation — Pathway Frequency Across Variants",
    subtitle = "How many of the 7 variants include each pathway in their final shortlist",
    x        = "Number of variants containing the pathway",
    y        = "Number of distinct pathway IDs",
    caption  = paste0(
      "Includes all 7 ablation variants (B_no_fgsea contributes 0 pathways).\n",
      "Pathways in 1 variant only are unique/sensitive; pathways in 5–7 variants are robust."
    )
  ) +
  theme_ablation() +
  theme(legend.position = "bottom")

ggsave(file.path(PLOT_DIR, "07_pathway_frequency.png"), p_freq,
       width = 7, height = 5, dpi = 150)
message("Saved: 07_pathway_frequency.png")

# ---------------------------------------------------------------------------
# Pathway variant count tables
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

readr::write_csv(pathway_variant_counts,
                 file.path(DATA_DIR, "all_pathways_variant_counts.csv"))
readr::write_csv(unique_pathways,
                 file.path(DATA_DIR, "unique_pathways.csv"))

# ---------------------------------------------------------------------------
# Main findings table
# ---------------------------------------------------------------------------

main_findings <- effect_df |>
  dplyr::mutate(
    component_removed = COMPONENT_REMOVED[variant],
    key_interpretation = KEY_INTERP[variant]
  ) |>
  dplyr::select(
    variant,
    component_removed,
    final_size,
    jaccard_vs_b_full = jaccard_full,
    ids_lost_from_b_full = ids_lost,
    ids_gained_vs_b_full = ids_gained,
    key_interpretation
  ) |>
  dplyr::mutate(jaccard_vs_b_full = round(jaccard_vs_b_full, 2))

readr::write_csv(main_findings, file.path(DATA_DIR, "main_findings_table.csv"))

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

cat("\n=== FUNNEL SUMMARY ===\n")
print(as.data.frame(funnel_rows), row.names = FALSE)

cat("\n=== PAIRWISE JACCARD (rounded) ===\n")
print(round(jac_mat, 2))

cat("\n=== MAIN FINDINGS TABLE ===\n")
print(as.data.frame(main_findings[, c("variant","final_size","jaccard_vs_b_full",
                                       "ids_lost_from_b_full","ids_gained_vs_b_full")]),
      row.names = FALSE)

cat("\n=== PATHWAY FREQUENCY ===\n")
print(as.data.frame(pathway_freq |> dplyr::count(n_variants)),
      row.names = FALSE)

message("\nUnique pathways (in exactly 1 variant): ", nrow(unique_pathways))
message("Distinct pathway IDs across all variants: ", nrow(pathway_variant_counts))
message("\nAll outputs written to:")
message("  ", PLOT_DIR)
message("  ", DATA_DIR)
