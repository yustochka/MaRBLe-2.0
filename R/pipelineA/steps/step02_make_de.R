# R/pipelineA/steps/step02_make_de.R
# Pipeline A — Step 02: Build the DE gene list (intersected with universe)
#
# Reads:
#   paths$step01/universe_genes.txt   (written by step01_make_universe in the same run)
#   cfg$dataset$statistics_file
#
# Column overrides (optional, in config under dataset:):
#   gene_col            exact column name for gene symbols
#   p_col               raw p-value column (fallback)
#   p_adjusted_preferred  prefer adj. p-value column when available (default TRUE)
#   p_col_adj           adj. p-value column name (default "adj. p-value")
#   log2fc_col          log2FC column for the adaptive tightening ladder
#
# Adaptive DE policy (under pipelineA in config):
#   target_de_frac_min / target_de_frac_max   target band for DE fraction
#   de_log2fc_ladder    |log2FC| thresholds tried in order (e.g. [0.5, 1.0])
#   de_p_ladder         p-cutoff tightening steps (e.g. [0.05, 0.01, 0.001])
#
# Writes into paths$step02:
#   de_genes_rawp05.txt   (final DE gene list)
#   de_policy.yml         (calibration metadata: thresholds chosen + history)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

step02_make_de <- function(cfg, paths, ctx) {

  # ---- Inputs ----------------------------------------------------------------
  UNIVERSE_FILE <- file.path(paths$step01, "universe_genes.txt")
  STATS_FILE    <- cfg$dataset$statistics_file
  DE_ALPHA      <- 0.05

  if (!file.exists(UNIVERSE_FILE)) {
    stop(
      "Universe file not found: ", UNIVERSE_FILE,
      "\nRun step01_make_universe() before step02_make_de()."
    )
  }

  universe <- readr::read_lines(UNIVERSE_FILE) |> unique()
  stats    <- readr::read_delim(STATS_FILE, delim = NULL, show_col_types = FALSE)

  # ---- Column detection -------------------------------------------------------
  nms <- names(stats)

  # Exact match after lowercasing column names in the file
  first_match <- function(patterns, original_names) {
    idx <- which(tolower(original_names) %in% patterns)
    if (length(idx) == 0) return(NA_character_)
    original_names[idx[1]]
  }

  # Fuzzy match: strip non-alphanumeric before comparing.
  # "Gene Symbol" -> "genesymbol", "gene_symbol" -> "genesymbol" — both match
  # pattern "genesymbol".
  first_match_norm <- function(patterns, original_names) {
    norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
    idx  <- which(norm(original_names) %in% norm(patterns))
    if (length(idx) == 0) return(NA_character_)
    original_names[idx[1]]
  }

  # -- Gene column: config override > exact match > fuzzy match ----------------
  if (!is.null(cfg$dataset$gene_col) && nzchar(cfg$dataset$gene_col)) {
    gcol <- cfg$dataset$gene_col
    if (!gcol %in% nms)
      stop("Configured gene_col '", gcol, "' not found in statistics file.\n",
           "  Available columns: ", paste(nms, collapse = ", "))
    gcol_source <- "config"
  } else {
    g_exact <- c("gene", "genes", "symbol", "geneid", "gene_id",
                 "external_gene_name", "hgnc_symbol")
    gcol <- first_match(g_exact, nms)
    if (is.na(gcol))
      gcol <- first_match_norm(c("genesymbol", "gene_symbol", "hgncsymbol"), nms)
    gcol_source <- "auto-detected"
  }

  # -- P-value column: config override > exact match > regex fallback ----------
  if (!is.null(cfg$dataset$p_col) && nzchar(cfg$dataset$p_col)) {
    pcol <- cfg$dataset$p_col
    if (!pcol %in% nms)
      stop("Configured p_col '", pcol, "' not found in statistics file.\n",
           "  Available columns: ", paste(nms, collapse = ", "))
    pcol_source <- "config"
  } else {
    # Raw p-value patterns first, adjusted p-value patterns as lower-priority fallback
    p_exact <- c("p", "pvalue", "p.value", "p-value", "raw_p", "p_val", "pval",
                 "pr(>|t|)", "padj", "p.adjust", "fdr", "adj.p.value")
    pcol <- first_match(p_exact, nms)
    if (is.na(pcol))
      pcol <- first_match_norm(c("adjpvalue", "pvalueadj"), nms)   # "adj. p-value" etc.
    if (is.na(pcol)) {
      # final fallback: regex on column names, validated as numeric [0,1]
      lnms     <- tolower(nms)
      rx       <- "(^p$)|(^p[._-]?value$)|(^pval$)|^pr\\(>\\|t\\|\\)$"
      cand_idx <- grep(rx, lnms, perl = TRUE)
      good     <- vapply(stats[cand_idx], function(x)
        is.numeric(x) && all(is.na(x) | (x >= 0 & x <= 1)),
        logical(1))
      cand_idx <- cand_idx[good]
      if (length(cand_idx)) pcol <- nms[cand_idx[1]]
    }
    pcol_source <- "auto-detected"
  }

  if (is.na(pcol)) {
    cat("Available columns:\n"); print(nms)
    stop("Couldn't find a p-value column. Set cfg$dataset$p_col to your exact header.")
  }
  if (is.na(gcol)) {
    cat("Available columns:\n"); print(nms)
    stop("Couldn't find a gene symbol column. Set cfg$dataset$gene_col to your exact header.")
  }

  cat(sprintf("Gene column    (%s): '%s'\n", gcol_source, gcol))
  cat(sprintf("P-value column (%s): '%s'\n", pcol_source, pcol))

  # ---- CHANGE 1: Adjusted p-value preference ---------------------------------
  # If p_adjusted_preferred is TRUE (default) and p_col_adj exists in the file,
  # use it instead of the raw p-value column detected above.
  p_adjusted_preferred <- if (!is.null(cfg$dataset$p_adjusted_preferred))
    isTRUE(cfg$dataset$p_adjusted_preferred) else TRUE
  p_col_adj_cfg <- if (!is.null(cfg$dataset$p_col_adj) && nzchar(cfg$dataset$p_col_adj))
    cfg$dataset$p_col_adj else "adj. p-value"
  if (p_adjusted_preferred && p_col_adj_cfg %in% nms) {
    active_pcol <- p_col_adj_cfg
  } else {
    active_pcol <- pcol
  }
  cat(sprintf("Using p column: '%s' (adjusted_preferred=%s)\n",
              active_pcol, p_adjusted_preferred))

  # ---- Adaptive DE policy -------------------------------------------------------
  # Target band and tightening ladders come from pipelineA config.
  # The policy is dataset-agnostic: it depends only on these config values and
  # which columns are present in the statistics file.

  target_min <- if (!is.null(cfg$pipelineA$target_de_frac_min))
    suppressWarnings(as.numeric(cfg$pipelineA$target_de_frac_min)) else 0.05
  target_max <- if (!is.null(cfg$pipelineA$target_de_frac_max))
    suppressWarnings(as.numeric(cfg$pipelineA$target_de_frac_max)) else 0.25

  # Log2FC column (used for L1/L2 ladder steps)
  lfc_col_cfg   <- if (!is.null(cfg$dataset$log2fc_col) && nzchar(cfg$dataset$log2fc_col))
    cfg$dataset$log2fc_col else "log2FC"
  lfc_available <- lfc_col_cfg %in% nms

  # Ladders from config; defaults match the specification
  lfc_ladder_cfg <- if (!is.null(cfg$pipelineA$de_log2fc_ladder))
    sort(unique(suppressWarnings(as.numeric(
      unlist(cfg$pipelineA$de_log2fc_ladder)))))
  else c(0.5, 1.0)
  lfc_ladder_cfg <- lfc_ladder_cfg[!is.na(lfc_ladder_cfg)]   # drop NAs (YAML nulls)

  p_ladder_cfg <- if (!is.null(cfg$pipelineA$de_p_ladder))
    sort(unique(suppressWarnings(as.numeric(
      unlist(cfg$pipelineA$de_p_ladder)))), decreasing = FALSE)
  else c(0.05, 0.01, 0.001)
  p_ladder_cfg <- p_ladder_cfg[!is.na(p_ladder_cfg)]

  # Helper: compute DE genes for given thresholds, universe-intersected
  genes <- as.character(stats[[gcol]])
  p_vec <- suppressWarnings(as.numeric(stats[[active_pcol]]))
  lfc_vec <- if (lfc_available)
    suppressWarnings(as.numeric(stats[[lfc_col_cfg]])) else NULL

  de_from_thresholds <- function(p_cut, lfc_min) {
    kp <- !is.na(p_vec) & p_vec < p_cut
    if (!is.na(lfc_min) && lfc_available) {
      kl <- !is.na(lfc_vec) & abs(lfc_vec) >= lfc_min
      kp <- kp & kl
    }
    sort(intersect(unique(genes[kp]), universe))
  }

  # ---- Base state ------------------------------------------------------------
  N           <- length(universe)
  current_p   <- DE_ALPHA
  current_lfc <- NA_real_
  de_genes    <- de_from_thresholds(current_p, current_lfc)
  n           <- length(de_genes)
  de_frac     <- if (N > 0) n / N else 0

  cat(sprintf("Base DE: p < %.3f, lfc = none -> n=%d, fraction=%.3f\n",
              current_p, n, de_frac))

  tightening_applied <- character(0)
  target_reached     <- de_frac <= target_max
  ladder_exhausted   <- FALSE
  ladder_history     <- list(
    list(step = "base", p_cutoff = current_p, lfc_min = NA_real_,
         n = n, de_fraction = de_frac)
  )

  # ---- Tightening ladder (only if fraction > target_max) --------------------
  if (!target_reached) {
    cat(sprintf("DE fraction %.3f > target_max %.2f — applying tightening ladder\n",
                de_frac, target_max))

    if (!lfc_available && length(lfc_ladder_cfg) > 0)
      cat(sprintf("  log2FC column '%s' not found — skipping L1/L2 (log2FC) steps\n",
                  lfc_col_cfg))

    # L1, L2 — log2FC tightening
    if (lfc_available) {
      for (lfc_thr in lfc_ladder_cfg) {
        if (de_frac <= target_max) break
        current_lfc <- lfc_thr
        de_genes    <- de_from_thresholds(current_p, current_lfc)
        n           <- length(de_genes)
        de_frac     <- if (N > 0) n / N else 0
        step_label  <- sprintf("L%d: |%s| >= %.1f",
                               length(tightening_applied) + 1, lfc_col_cfg, lfc_thr)
        tightening_applied <- c(tightening_applied, step_label)
        cat(sprintf("  %s -> n=%d, fraction=%.3f\n", step_label, n, de_frac))
        ladder_history[[length(ladder_history) + 1]] <-
          list(step = step_label, p_cutoff = current_p, lfc_min = current_lfc,
               n = n, de_fraction = de_frac)
      }
    }

    # L3, L4 — p-value tightening (skip entries >= current_p)
    for (p_thr in p_ladder_cfg[p_ladder_cfg < current_p]) {
      if (de_frac <= target_max) break
      current_p  <- p_thr
      de_genes   <- de_from_thresholds(current_p, current_lfc)
      n          <- length(de_genes)
      de_frac    <- if (N > 0) n / N else 0
      step_label <- sprintf("L%d: p < %.3f",
                            length(tightening_applied) + 1, p_thr)
      tightening_applied <- c(tightening_applied, step_label)
      cat(sprintf("  %s -> n=%d, fraction=%.3f\n", step_label, n, de_frac))
      ladder_history[[length(ladder_history) + 1]] <-
        list(step = step_label, p_cutoff = current_p, lfc_min = current_lfc,
             n = n, de_fraction = de_frac)
    }

    target_reached   <- de_frac <= target_max
    ladder_exhausted <- !target_reached
    if (ladder_exhausted)
      cat(sprintf("  Ladder exhausted; DE fraction %.3f still above target %.2f (proceeding)\n",
                  de_frac, target_max))
    else
      cat(sprintf("  Target reached: fraction=%.3f\n", de_frac))

  } else if (de_frac < target_min) {
    cat(sprintf("DE fraction %.3f < target_min %.2f (keeping conservative; no loosening)\n",
                de_frac, target_min))
  } else {
    cat(sprintf("DE fraction %.3f within target band [%.2f, %.2f]\n",
                de_frac, target_min, target_max))
  }

  de_base <- de_genes   # already sorted and universe-intersected

  # ---- Write de_policy.yml (calibration metadata) ----------------------------
  policy_reason <- if (length(tightening_applied) > 0) {
    if (ladder_exhausted)
      "DE fraction above target_de_frac_max; ladder exhausted without reaching target"
    else
      "DE fraction above target_de_frac_max; tightening applied"
  } else if (de_frac < target_min) {
    "DE fraction below target_de_frac_min (no loosening applied)"
  } else {
    "DE fraction within target band; no tightening needed"
  }

  policy <- list(
    dataset              = cfg$dataset$name,
    universe_N           = N,
    de_count_n           = length(de_base),
    de_fraction          = round(de_frac, 4),
    p_col_used           = active_pcol,
    p_adjusted_preferred = p_adjusted_preferred,
    final_p_cutoff       = current_p,
    final_lfc_min        = if (is.na(current_lfc)) "none" else current_lfc,
    lfc_col              = lfc_col_cfg,
    lfc_col_available    = lfc_available,
    target_de_frac_min   = target_min,
    target_de_frac_max   = target_max,
    target_reached       = target_reached,
    tightening_steps     = as.list(tightening_applied),
    reason               = policy_reason,
    ladder_history       = ladder_history
  )

  policy_file <- file.path(paths$step02, "de_policy.yml")
  yaml::write_yaml(policy, policy_file)

  # ---- Outputs ---------------------------------------------------------------
  out_file <- file.path(paths$step02, "de_genes_rawp05.txt")
  writeLines(de_base, out_file)

  cat("\n=== Pipeline A Step 02: DE genes ===\n")
  cat("P col:    '", active_pcol, "' (adj preferred: ", p_adjusted_preferred, ")\n", sep = "")
  cat("P cutoff: ", current_p, "\n", sep = "")
  cat("LFC min:  ", if (is.na(current_lfc)) "none" else sprintf("%.1f", current_lfc), "\n", sep = "")
  cat("Universe: ", N, " genes\n", sep = "")
  cat("DE genes: ", length(de_base), " (fraction: ", round(de_frac, 3), ")\n", sep = "")
  if (length(tightening_applied) > 0)
    cat("Tightening:", paste(tightening_applied, collapse = " | "), "\n")
  cat("Wrote:    ", out_file, "\n")
  cat("Policy:   ", policy_file, "\n\n", sep = "")

  invisible(list(de_genes_file  = out_file,
                 de_policy_file = policy_file))
}
