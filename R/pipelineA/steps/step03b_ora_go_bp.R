# R/pipelineA/steps/step03b_ora_go_bp.R
# Pipeline A — Step 03b: ORA on GO Biological Process with data-driven filters
#
# Reads:
#   paths$step01/universe_genes.txt        (written by step01_make_universe)
#   paths$step02/de_genes_rawp05.txt       (written by step02_make_de)
#
# Writes into paths$step03_go_bp:
#   ora_go_bp_all.csv
#   ora_go_bp_filtered.csv
#   run_note.txt

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

step03b_ora_go_bp <- function(cfg, paths, ctx) {

  # ---- Inputs ----------------------------------------------------------------
  UNIVERSE_FILE <- file.path(paths$step01, "universe_genes.txt")
  DE_FILE       <- file.path(paths$step02, "de_genes_rawp05.txt")

  stopifnot(file.exists(UNIVERSE_FILE), file.exists(DE_FILE))

  universe <- readr::read_lines(UNIVERSE_FILE)
  de       <- readr::read_lines(DE_FILE)
  N <- length(universe)
  n <- length(de)
  cat(sprintf("GO-BP ORA — N=%d, n=%d\n", N, n))

  # ---- Outputs ---------------------------------------------------------------
  OUT_ALL  <- file.path(paths$step03_go_bp, "ora_go_bp_all.csv")
  OUT_FILT <- file.path(paths$step03_go_bp, "ora_go_bp_filtered.csv")
  OUT_NOTE <- file.path(paths$step03_go_bp, "run_note.txt")

  if (length(de) == 0) {
    msg <- "DE gene list is empty (step02 produced no genes). Downstream steps will be skipped."
    writeLines(paste0("[GO-BP ORA]\n", msg), OUT_NOTE)
    readr::write_csv(tibble::tibble(), OUT_ALL)
    readr::write_csv(tibble::tibble(), OUT_FILT)
    cat(msg, "\n")
    return(invisible(list(ora_all_file = OUT_ALL, ora_filtered_file = OUT_FILT)))
  }

  de_fraction <- n / N
  ora_mode <- if (!is.null(cfg$pipelineA$ora_min_enrich_ratio_mode) &&
                  nzchar(cfg$pipelineA$ora_min_enrich_ratio_mode))
    cfg$pipelineA$ora_min_enrich_ratio_mode else "adaptive"
  rmin <- if (ora_mode == "adaptive") {
    if      (de_fraction <= 0.10) 1.5
    else if (de_fraction <= 0.30) 1.2
    else                          1.0
  } else {
    fixed_r <- suppressWarnings(as.numeric(cfg$pipelineA$ora_min_enrich_ratio))
    if (!is.na(fixed_r)) fixed_r else 1.5
  }
  cat(sprintf("ORA enrich-ratio threshold: t=%.1f (mode=%s, de_fraction=%.3f)\n",
              rmin, ora_mode, de_fraction))

  # ---- ORA on GO:BP ----------------------------------------------------------
  ego <- enrichGO(
    gene          = de,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1, qvalueCutoff = 1,
    minGSSize     = 1,  maxGSSize = 100000,
    readable      = FALSE
  )

  if (is.null(ego) || nrow(ego@result) == 0) {
    msg <- "No GO:BP terms returned by enrichGO. Check universe/DE IDs or relax filters."
    writeLines(paste0("[GO-BP ORA]\n", msg), OUT_NOTE)
    readr::write_csv(tibble::tibble(), OUT_ALL)
    readr::write_csv(tibble::tibble(), OUT_FILT)
    cat(msg, "\n")
    return(invisible(list(ora_all_file = OUT_ALL, ora_filtered_file = OUT_FILT)))
  }

  # ---- Unpack ratios, add E[k] and enrichment ratio -------------------------
  res <- ego@result |>
    dplyr::arrange(p.adjust) |>
    tidyr::separate_wider_delim(GeneRatio, "/", names = c("k", "n_in_table"), cols_remove = FALSE) |>
    tidyr::separate_wider_delim(BgRatio,   "/", names = c("K", "N_in_table"), cols_remove = FALSE) |>
    dplyr::mutate(dplyr::across(c(k, n_in_table, K, N_in_table), as.numeric),
                  Ek        = (n * K) / N,
                  enr_ratio = k / pmax(Ek, .Machine$double.eps)) |>
    dplyr::relocate(ID, Description, k, K, GeneRatio, BgRatio, Ek, enr_ratio, pvalue, p.adjust, Count)

  readr::write_csv(res, OUT_ALL)

  # ---- Data-driven filters ---------------------------------------------------
  K_q   <- quantile(res$K, probs = c(.10, .50, .90, .95), na.rm = TRUE)
  K_min <- max(15L, as.integer(K_q[1]))
  K_max <- min(500L, as.integer(K_q[4]))

  if (K_min > K_max) {
    K_min <- max(3L, K_max)
    cat(sprintf("  Warning: K_min > K_max after percentile computation; K_min clamped to %d\n", K_min))
  }

  kmin  <- 3L

  # FDR and fallback params (from cfg; safe defaults if keys absent)
  go_fdr_strict  <- if (!is.null(cfg$pipelineA$go_fdr_cutoff_strict))
    suppressWarnings(as.numeric(cfg$pipelineA$go_fdr_cutoff_strict)) else 0.05
  go_fdr_relaxed <- if (!is.null(cfg$pipelineA$go_fdr_cutoff_relaxed))
    suppressWarnings(as.numeric(cfg$pipelineA$go_fdr_cutoff_relaxed)) else 0.20
  go_fallback_mode <- if (!is.null(cfg$pipelineA$go_fallback_mode) &&
                          nzchar(cfg$pipelineA$go_fallback_mode))
    cfg$pipelineA$go_fallback_mode else "relax_fdr_then_topn"
  go_fallback_topn <- if (!is.null(cfg$pipelineA$go_fallback_topn))
    as.integer(cfg$pipelineA$go_fallback_topn) else 50L

  # Strict filter
  res_f <- res |>
    dplyr::filter(p.adjust < go_fdr_strict,
                  dplyr::between(K, K_min, K_max),
                  k >= kmin,
                  enr_ratio >= rmin)

  # ---- Fallback (only when strict returns 0 rows) ----------------------------
  fallback_used <- "none"
  if (nrow(res_f) == 0 && go_fallback_mode != "none") {
    cat(sprintf("Strict GO filter (FDR<%.2f) returned 0 rows — trying fallback mode='%s'\n",
                go_fdr_strict, go_fallback_mode))
    if (grepl("relax_fdr", go_fallback_mode)) {
      res_f <- res |>
        dplyr::filter(p.adjust < go_fdr_relaxed,
                      dplyr::between(K, K_min, K_max),
                      k >= kmin,
                      enr_ratio >= rmin)
      if (nrow(res_f) > 0)
        fallback_used <- sprintf("relaxed_fdr (FDR<%.2f)", go_fdr_relaxed)
    }
    if (nrow(res_f) == 0 && grepl("topn", go_fallback_mode)) {
      res_f <- res |>
        dplyr::arrange(p.adjust, pvalue) |>
        dplyr::slice_head(n = go_fallback_topn)
      if (nrow(res_f) > 0)
        fallback_used <- sprintf("topN (n=%d, sorted by FDR)", go_fallback_topn)
    }
    if (fallback_used != "none")
      cat(sprintf("  Fallback applied: %s -> %d rows kept\n", fallback_used, nrow(res_f)))
    else
      cat("  All fallbacks exhausted — filtered output remains empty.\n")
  }

  readr::write_csv(res_f, OUT_FILT)

  # ---- Run note --------------------------------------------------------------
  de_policy_summary <- tryCatch({
    dp <- yaml::read_yaml(file.path(paths$step02, "de_policy.yml"))
    sprintf("p_col=%s, p_cutoff=%.3f, lfc_min=%s, steps=%s",
            dp$p_col_used, dp$final_p_cutoff, dp$final_lfc_min,
            if (length(dp$tightening_steps) == 0) "none"
            else paste(dp$tightening_steps, collapse = "|"))
  }, error = function(e) "de_policy.yml not found")

  note <- sprintf(paste0(
    "[GO-BP ORA]\n",
    "N=%d, n=%d, de_fraction=%.3f\n",
    "DE policy: %s\n",
    "K_min..K_max=[%d,%d] (from GO-BP K distribution)\n",
    "Strict filters: FDR<%.2f, k>=%d, k/E[k]>=%.1f (mode=%s)\n",
    "Fallback mode: %s  |  Fallback used: %s\n",
    "Tested=%d, Kept=%d (final filtered)\n",
    "Files: %s, %s\n"),
    N, n, de_fraction, de_policy_summary,
    K_min, K_max, go_fdr_strict, kmin, rmin, ora_mode,
    go_fallback_mode, fallback_used,
    nrow(res), nrow(res_f),
    OUT_ALL, OUT_FILT
  )

  writeLines(note, OUT_NOTE)
  cat(note)

  invisible(list(
    ora_all_file      = OUT_ALL,
    ora_filtered_file = OUT_FILT
  ))
}
