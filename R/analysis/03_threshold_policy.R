# R/analysis/03_threshold_policy.R
# Reads existing analysis outputs and produces a per-dataset threshold
# recommendation table.
#
# Inputs:
#   results/analysis/dataset_summary.csv
#   results/analysis/sweep_pipelineA.csv
#   results/analysis/expressed_fraction.csv
#
# Output:
#   results/analysis/threshold_policy.csv
#
# Run from project root:
#   Rscript R/analysis/03_threshold_policy.R
#
# ---- Classification logic (applied in order) --------------------------------
#
#  LOW_COVERAGE  : n_universe < 2000  OR  med_frac_in_matrix < 0.20
#                  (focused panel; pathway members mostly unmeasured; ORA unreliable)
#
#  NO_OVERLAP    : n_consensus_pairs == 0 across ALL sweep combos
#                  AND final_source never reaches "tier1" in any sweep row
#                  (GO↔WP gene overlap structurally absent at all tested thresholds)
#
#  BORDERLINE    : baseline (TAU=0.70, CONS_J=0.15) gives final_source != "tier1"
#                  BUT at least one sweep row achieves final_source == "tier1"
#                  (Tier1 recoverable by relaxing CONS_JACCARD_MIN)
#
#  HEALTHY       : baseline gives final_source == "tier1"
#
# ---- Threshold recommendation rules ----------------------------------------
#
#  HEALTHY     : TAU=0.70, CONS_J=0.15  (keep defaults; Tier1 robust)
#  BORDERLINE  : TAU=0.70, CONS_J=best value that achieves final_source=="tier1"
#  NO_OVERLAP  : TAU=0.70, CONS_J=0.10  (most permissive tested; won't help consensus
#                                         but does not hurt; Tier2 is the output)
#  LOW_COVERAGE: TAU=0.50, CONS_J=0.10  (least strict; stability near zero)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(readr))

ROOT <- getwd()

# ---- Constants: baseline thresholds -----------------------------------------
BASELINE_TAU    <- 0.70
BASELINE_CONS_J <- 0.15
LOW_UNIVERSE    <- 2000L
LOW_COVERAGE_FRAC <- 0.20

# ---- Read inputs -------------------------------------------------------------
summary_csv <- file.path(ROOT, "results", "analysis", "dataset_summary.csv")
sweep_csv   <- file.path(ROOT, "results", "analysis", "sweep_pipelineA.csv")
expr_csv    <- file.path(ROOT, "results", "analysis", "expressed_fraction.csv")

stopifnot(
  file.exists(summary_csv),
  file.exists(sweep_csv),
  file.exists(expr_csv)
)

summary_df <- read_csv(summary_csv, show_col_types = FALSE)
sweep_df   <- read_csv(sweep_csv,   show_col_types = FALSE)
expr_df    <- read_csv(expr_csv,    show_col_types = FALSE)

# ---- Derive per-dataset signals from each input ----------------------------

# 1) From dataset_summary.csv — pipelineA rows only
sumA <- summary_df[summary_df$pipeline == "pipelineA", ]

# 2) From expressed_fraction.csv — pipelineA rows, median per dataset
expr_medA <- do.call(rbind, lapply(
  unique(expr_df$dataset[expr_df$pipeline == "pipelineA"]),
  function(d) {
    sub <- expr_df[expr_df$dataset == d & expr_df$pipeline == "pipelineA", ]
    data.frame(
      dataset            = d,
      med_frac_in_matrix = median(sub$frac_in_matrix, na.rm = TRUE),
      med_frac_low_expr  = median(sub$frac_low_expr,  na.rm = TRUE),
      stringsAsFactors   = FALSE
    )
  }
))

# 3) From sweep_pipelineA.csv — compute per-dataset sweep signals
sweep_signals <- do.call(rbind, lapply(
  unique(sweep_df$dataset),
  function(d) {
    sub <- sweep_df[sweep_df$dataset == d, ]

    # Baseline row (TAU=0.70, CONS_J=0.15)
    baseline <- sub[
      abs(sub$tau_stability    - BASELINE_TAU)    < 1e-9 &
      abs(sub$cons_jaccard_min - BASELINE_CONS_J) < 1e-9, ]

    tier1_at_baseline   <- if (nrow(baseline) == 1) as.integer(baseline$n_tier1)     else NA_integer_
    baseline_fs         <- if (nrow(baseline) == 1) as.character(baseline$final_source) else NA_character_
    # Stable counts at baseline (from sweep columns n_go_stable, n_wp_stable)
    n_stable_go_baseline <- if (nrow(baseline) == 1) as.integer(baseline$n_go_stable) else NA_integer_
    n_stable_wp_baseline <- if (nrow(baseline) == 1) as.integer(baseline$n_wp_stable) else NA_integer_

    # Tier1 ever achieved (final_source == "tier1" implies n_tier1 >= 5 via step98 logic)
    tier1_rows    <- sub[!is.na(sub$final_source) & sub$final_source == "tier1", ]
    tier1_achievable <- nrow(tier1_rows) > 0

    # Lowest CONS_J for which Tier1 is achieved (at any TAU)
    min_cons_j_achieving_tier1 <- if (tier1_achievable) min(tier1_rows$cons_jaccard_min) else NA_real_

    # Lowest TAU at that CONS_J — informational only, not used for recommendations
    min_tau_achieving_tier1 <- if (tier1_achievable) {
      min(tier1_rows$tau_stability[
        tier1_rows$cons_jaccard_min == min_cons_j_achieving_tier1])
    } else NA_real_

    max_tier1      <- max(sub$n_tier1,          na.rm = TRUE)
    max_cons_pairs <- max(sub$n_consensus_pairs, na.rm = TRUE)

    data.frame(
      dataset                    = d,
      tier1_at_baseline          = tier1_at_baseline,
      baseline_final_source      = baseline_fs,
      n_stable_go_baseline       = n_stable_go_baseline,
      n_stable_wp_baseline       = n_stable_wp_baseline,
      tier1_achievable           = tier1_achievable,
      min_tau_achieving_tier1    = min_tau_achieving_tier1,
      min_cons_j_achieving_tier1 = min_cons_j_achieving_tier1,
      max_tier1                  = max_tier1,
      max_consensus_pairs        = max_cons_pairs,
      stringsAsFactors           = FALSE
    )
  }
))

# ---- Merge signals into one table per dataset --------------------------------
datasets <- unique(sumA$dataset)

policy_rows <- lapply(datasets, function(d) {

  # --- Pull signals ---
  s    <- sumA[sumA$dataset == d, ]
  expr <- expr_medA[expr_medA$dataset == d, ]
  sw   <- sweep_signals[sweep_signals$dataset == d, ]
  sub  <- sweep_df[sweep_df$dataset == d, ]        # raw sweep rows for this dataset

  n_universe          <- if (nrow(s) > 0) s$n_universe[1]        else NA_integer_
  de_fraction         <- if (nrow(s) > 0) s$de_fraction[1]       else NA_real_
  mean_stab_go        <- if (nrow(s) > 0) s$mean_stability_go[1] else NA_real_
  mean_stab_wp        <- if (nrow(s) > 0) s$mean_stability_wp[1] else NA_real_
  baseline_cons_pairs <- if (nrow(s) > 0) s$consensus_pairs[1]   else NA_integer_

  med_frac_in_matrix  <- if (nrow(expr) > 0) expr$med_frac_in_matrix[1] else NA_real_
  med_frac_low_expr   <- if (nrow(expr) > 0) expr$med_frac_low_expr[1]  else NA_real_

  tier1_at_baseline          <- if (nrow(sw) > 0) sw$tier1_at_baseline[1]          else NA_integer_
  baseline_fs                <- if (nrow(sw) > 0) sw$baseline_final_source[1]       else NA_character_
  n_stable_go_baseline       <- if (nrow(sw) > 0) sw$n_stable_go_baseline[1]       else NA_integer_
  n_stable_wp_baseline       <- if (nrow(sw) > 0) sw$n_stable_wp_baseline[1]       else NA_integer_
  tier1_achievable           <- if (nrow(sw) > 0) isTRUE(sw$tier1_achievable[1])   else FALSE
  min_tau_achieving_tier1    <- if (nrow(sw) > 0) sw$min_tau_achieving_tier1[1]    else NA_real_
  min_cons_j_achieving_tier1 <- if (nrow(sw) > 0) sw$min_cons_j_achieving_tier1[1] else NA_real_
  max_consensus_pairs        <- if (nrow(sw) > 0) sw$max_consensus_pairs[1]        else NA_integer_

  # --- Classify (in priority order) ---
  n_universe_ok <- !is.na(n_universe) && n_universe >= LOW_UNIVERSE
  frac_ok       <- is.na(med_frac_in_matrix) || med_frac_in_matrix >= LOW_COVERAGE_FRAC

  if (!n_universe_ok || !frac_ok) {
    classification <- "LOW_COVERAGE"

  } else if (!tier1_achievable && (!is.na(max_consensus_pairs) && max_consensus_pairs == 0)) {
    classification <- "NO_OVERLAP"

  } else if (!isTRUE(baseline_fs == "tier1") && tier1_achievable) {
    classification <- "BORDERLINE"

  } else {
    classification <- "HEALTHY"
  }

  # --- Assign recommended thresholds ---
  rec_tau <- switch(classification,
    HEALTHY      = 0.70,
    BORDERLINE   = 0.70,
    NO_OVERLAP   = 0.70,
    LOW_COVERAGE = 0.50
  )

  rec_cons_j <- switch(classification,
    HEALTHY      = 0.15,
    BORDERLINE   = if (!is.na(min_cons_j_achieving_tier1)) min_cons_j_achieving_tier1 else 0.10,
    NO_OVERLAP   = 0.10,
    LOW_COVERAGE = 0.10
  )

  allow_tier2 <- TRUE   # always allow Tier2 as fallback

  confidence <- switch(classification,
    HEALTHY      = "high",
    BORDERLINE   = "medium",
    NO_OVERLAP   = "low",
    LOW_COVERAGE = "unreliable"
  )

  # --- n_tier1 at the recommended threshold combination ---
  rec_row <- sub[
    abs(sub$tau_stability    - rec_tau)    < 1e-9 &
    abs(sub$cons_jaccard_min - rec_cons_j) < 1e-9, ]
  tier1_at_recommended <- if (nrow(rec_row) == 1) as.integer(rec_row$n_tier1) else NA_integer_

  # --- Notes ---
  notes <- switch(classification,
    HEALTHY = paste0(
      "Tier1 robust across all sweep combinations (max=", sw$max_tier1[1], "); keep defaults"
    ),
    BORDERLINE = paste0(
      "Tier1 achievable at CONS_J=", min_cons_j_achieving_tier1,
      " (Tier1 at recommended thresholds=", tier1_at_recommended,
      "); baseline CONS_J=0.15 yields only ", tier1_at_baseline,
      " (below step98 min of 5)"
    ),
    NO_OVERLAP = paste0(
      "Tier1 structurally impossible: n_consensus_pairs=0 in all ",
      nrow(sub), " sweep combos (3 TAU x 3 CONS_J tested); ",
      "stable WP reps at baseline=", n_stable_wp_baseline,
      " but no GO-WP gene overlap detected; ",
      "Tier2 is the correct final output for this dataset"
    ),
    LOW_COVERAGE = paste0(
      "Focused panel: universe=", n_universe, " genes; ",
      "median frac_in_matrix=", round(med_frac_in_matrix, 3),
      "; mean_stability_go=", sprintf("%.5f", mean_stab_go),
      "; ORA based on <15% of pathway genes — results unreliable"
    )
  )

  data.frame(
    dataset                    = d,
    n_universe                 = n_universe,
    de_fraction                = round(de_fraction, 4),
    mean_stability_go          = round(mean_stab_go, 4),
    mean_stability_wp          = round(mean_stab_wp, 4),
    n_stable_go_baseline       = n_stable_go_baseline,
    n_stable_wp_baseline       = n_stable_wp_baseline,
    baseline_consensus_pairs   = baseline_cons_pairs,
    med_frac_in_matrix         = round(med_frac_in_matrix, 4),
    med_frac_low_expr          = round(med_frac_low_expr, 4),
    tier1_at_baseline          = tier1_at_baseline,
    tier1_achievable           = tier1_achievable,
    min_tau_achieving_tier1    = min_tau_achieving_tier1,
    min_cons_j_achieving_tier1 = min_cons_j_achieving_tier1,
    max_consensus_pairs        = max_consensus_pairs,
    classification             = classification,
    recommended_TAU_STABILITY    = rec_tau,
    recommended_CONS_JACCARD_MIN = rec_cons_j,
    tier1_at_recommended       = tier1_at_recommended,
    allow_tier2                = allow_tier2,
    confidence_flag            = confidence,
    notes                      = notes,
    stringsAsFactors = FALSE
  )
})

policy <- do.call(rbind, policy_rows)
rownames(policy) <- NULL

# ---- Write output ------------------------------------------------------------
out_dir  <- file.path(ROOT, "results", "analysis")
out_file <- file.path(out_dir, "threshold_policy.csv")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write_csv(policy, out_file)

# ---- Print per-dataset summary -----------------------------------------------
cat("=== 03_threshold_policy.R ===\n\n")
cat(sprintf("%-12s  %-14s  %7s  %5s  %5s  %4s  %4s  %6s  %5s  %4s  %6s\n",
            "dataset", "class", "uni", "stbGO", "stbWP", "stGO", "stWP",
            "medFrc", "T1bl", "T1rc", "recTAU"))
cat(strrep("-", 90), "\n", sep = "")

for (i in seq_len(nrow(policy))) {
  r <- policy[i, ]
  cat(sprintf("%-12s  %-14s  %7d  %5.3f  %5s  %4s  %4s  %6.3f  %5s  %4s  %5.2f/%-4.2f\n",
    r$dataset,
    r$classification,
    r$n_universe,
    r$mean_stability_go,
    if (!is.na(r$mean_stability_wp)) sprintf("%.3f", r$mean_stability_wp) else "  NA ",
    if (!is.na(r$n_stable_go_baseline)) as.character(r$n_stable_go_baseline) else " NA",
    if (!is.na(r$n_stable_wp_baseline)) as.character(r$n_stable_wp_baseline) else " NA",
    r$med_frac_in_matrix,
    if (!is.na(r$tier1_at_baseline))    as.character(r$tier1_at_baseline)    else " NA",
    if (!is.na(r$tier1_at_recommended)) as.character(r$tier1_at_recommended) else " NA",
    r$recommended_TAU_STABILITY,
    r$recommended_CONS_JACCARD_MIN
  ))
}

cat(strrep("-", 90), "\n\n", sep = "")

cat("Column legend:\n")
cat("  uni    = n_universe (measured ∩ pathway-mapped genes)\n")
cat("  stbGO  = mean GO bootstrap stability (all reps)\n")
cat("  stbWP  = mean WP bootstrap stability (all reps)\n")
cat("  stGO   = n_go_stable at baseline (TAU=0.70, CONS_J=0.15)\n")
cat("  stWP   = n_wp_stable at baseline\n")
cat("  medFrc = median frac_in_matrix (fraction of pathway genes measured)\n")
cat("  T1bl   = Tier1 count at baseline\n")
cat("  T1rc   = Tier1 count at recommended thresholds\n")
cat("  recTAU = recommended TAU / CONS_J\n")
cat("\n")

cat("Classifications:\n")
for (i in seq_len(nrow(policy))) {
  r <- policy[i, ]
  cat(sprintf("  %-12s  [%s]  confidence=%s\n",
              r$dataset, r$classification, r$confidence_flag))
  cat(sprintf("             -> TAU=%.2f  CONS_J=%.2f  allow_tier2=%s  tier1_at_recommended=%s\n",
              r$recommended_TAU_STABILITY, r$recommended_CONS_JACCARD_MIN,
              r$allow_tier2,
              if (!is.na(r$tier1_at_recommended)) as.character(r$tier1_at_recommended) else "NA"))
  cat(sprintf("             %s\n\n", r$notes))
}

cat("Output written to:", out_file, "\n")
