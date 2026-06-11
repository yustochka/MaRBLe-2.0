# R/pipelineA/steps/step98_export_final.R
# Pipeline A — Step 98: Export canonical final pathway list to FINAL/
#
# Canonical source priority:
#   1. Tier 1 (stable + consensus)        — when Tier 1 has >= 5 pathways
#   2. Tier 2 (stable only)               — when Tier 1 has < 5 pathways
#   3. All candidates (tiered_all_...)    — when Tier 1 + Tier 2 empty and
#                                           final_fallback_mode == "tier2_then_candidates"
#
# Reads:
#   paths$step07/tier1_final.csv
#   paths$step07/tier2_stable_only.csv      (fallback 1)
#   paths$step07/tiered_all_candidates.csv  (fallback 2)
#
# Config (pipelineA section):
#   final_fallback_mode:   "tier2_then_candidates" (default) | "tier2_only" | "none"
#   final_candidates_topn: null (export all) or integer (top-N by stability desc, p.adjust asc)
#
# Writes into ctx$run_dir/FINAL/:
#   FINAL.csv          chosen source
#   final_summary.txt  pipeline/dataset/run metadata + pathway count + thresholds
#
# Returns (invisibly):
#   final_file    path to FINAL.csv
#   summary_file  path to final_summary.txt

step98_export_final <- function(cfg, paths, ctx) {

  TIER1_FILE      <- file.path(paths$step07, "tier1_final.csv")
  TIER2_FILE      <- file.path(paths$step07, "tier2_stable_only.csv")
  CANDIDATES_FILE <- file.path(paths$step07, "tiered_all_candidates.csv")
  OUT_DIR         <- file.path(ctx$run_dir, "FINAL")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  # Config-driven fallback policy
  fallback_mode   <- if (!is.null(cfg$pipelineA$final_fallback_mode) &&
                         nzchar(cfg$pipelineA$final_fallback_mode))
    cfg$pipelineA$final_fallback_mode else "tier2_then_candidates"
  candidates_topn <- if (!is.null(cfg$pipelineA$final_candidates_topn))
    as.integer(cfg$pipelineA$final_candidates_topn) else NULL

  if (!file.exists(TIER1_FILE)) {
    stop("Canonical source not found: ", TIER1_FILE,
         "\nRun step07_make_tiers() before step98_export_final().")
  }

  tbl          <- readr::read_csv(TIER1_FILE, show_col_types = FALSE)
  final_source <- "Tier 1"
  source_file  <- TIER1_FILE

  # Fallback 1: Tier 2 (when Tier 1 has fewer than 5 pathways)
  tier1_n <- nrow(tbl)
  if (tier1_n < 5 && file.exists(TIER2_FILE)) {
    tbl2 <- readr::read_csv(TIER2_FILE, show_col_types = FALSE)
    if (nrow(tbl2) > 0) {
      tbl          <- tbl2
      final_source <- if (tier1_n == 0)
        "Tier 2 (Tier 1 empty)"
      else
        sprintf("Tier 2 (Tier 1 had only %d pathway%s)", tier1_n, if (tier1_n == 1) "" else "s")
      source_file  <- TIER2_FILE
      cat(if (tier1_n == 0)
        "Tier 1 is empty — falling back to Tier 2 (stable-only).\n"
      else
        sprintf("Tier 1 has only %d pathway%s (< 5) — falling back to Tier 2 (stable-only).\n",
                tier1_n, if (tier1_n == 1) "" else "s"))
    }
  }

  # Fallback 2: all candidates
  if (nrow(tbl) == 0 && fallback_mode == "tier2_then_candidates" &&
      file.exists(CANDIDATES_FILE)) {
    cands <- readr::read_csv(CANDIDATES_FILE, show_col_types = FALSE)
    if (nrow(cands) > 0) {
      if ("stability" %in% names(cands))
        cands <- cands |> dplyr::arrange(dplyr::desc(stability), p.adjust)
      if (!is.null(candidates_topn))
        cands <- dplyr::slice_head(cands, n = candidates_topn)
      tbl          <- cands
      final_source <- sprintf("candidates (Tier1+Tier2 empty; n=%d)", nrow(tbl))
      source_file  <- CANDIDATES_FILE
      cat("Tier 1 and Tier 2 are empty — falling back to all candidates.\n")
    }
  }

  OUT_FILE <- file.path(OUT_DIR, "FINAL.csv")
  readr::write_csv(tbl, OUT_FILE)

  cat(sprintf("FINAL source: %s\nPathways exported: %d\n", final_source, nrow(tbl)))

  # ---- Summary ---------------------------------------------------------------
  topn_line <- if (!is.null(candidates_topn))
    paste0("final_candidates_topn: ", candidates_topn, "\n") else ""

  summary_text <- paste0(
    "Pipeline:            pipelineA\n",
    "Dataset:             ", cfg$dataset$name, "\n",
    "Run ID:              ", ctx$run_id, "\n",
    "Timestamp:           ", ctx$timestamp, "\n",
    "FINAL source:        ", final_source, "\n",
    "Source file:         ", source_file, "\n",
    "Pathways exported:   ", nrow(tbl), "\n",
    "final_fallback_mode: ", fallback_mode, "\n",
    topn_line,
    "\nKey thresholds:\n",
    "  bootstrap B:          ", cfg$bootstrap$B, "\n",
    "  subsample fraction:   ", cfg$bootstrap$SUBSAMPLE_P, "\n",
    "  TAU_STABILITY:        ", cfg$thresholds$TAU_STABILITY, "\n",
    "  JACCARD_MIN:          ", cfg$thresholds$JACCARD_MIN, "\n",
    "  CONS_JACCARD_MIN:     ", cfg$thresholds$CONS_JACCARD_MIN, "\n",
    "  CONS_K_MIN:           ", cfg$thresholds$CONS_K_MIN, "\n"
  )

  SUMMARY_FILE <- file.path(OUT_DIR, "final_summary.txt")
  writeLines(summary_text, SUMMARY_FILE)

  invisible(list(
    final_file   = OUT_FILE,
    summary_file = SUMMARY_FILE
  ))
}
