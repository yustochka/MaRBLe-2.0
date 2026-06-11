# R/pipelineA/steps/step07_make_tiers.R
# Pipeline A — Step 07: Assign tiers to final candidates
#
# Reads:
#   paths$step06/final_shortlist.csv    (written by step06_bootstrap_consensus)
#
# Writes into paths$step07:
#   tier1_final.csv           stable AND consensus
#   tier2_stable_only.csv     stable but NOT in consensus
#   tiered_all_candidates.csv all candidates with tier label
#
# Returns:
#   tier1_file       path to tier1_final.csv
#   tier2_file       path to tier2_stable_only.csv
#   final_table_file path to tiered_all_candidates.csv

step07_make_tiers <- function(cfg, paths, ctx) {

  # ---- Input -----------------------------------------------------------------
  FINAL_FILE <- file.path(paths$step06, "final_shortlist.csv")
  stopifnot(file.exists(FINAL_FILE))

  final <- readr::read_csv(FINAL_FILE, show_col_types = FALSE)

  # Safety: make sure the columns we need exist
  needed  <- c("ID", "Description", "collection", "p.adjust", "k",
               "stability", "selected", "consensus", "keep_final")
  missing <- setdiff(needed, names(final))
  if (length(missing) > 0) {
    stop("final_shortlist.csv is missing columns: ", paste(missing, collapse = ", "))
  }

  # ---- Tier definitions (identical to original) ------------------------------
  # Tier 1 = stable AND consensus
  tier1 <- final |>
    dplyr::filter(keep_final == TRUE) |>
    dplyr::arrange(collection, p.adjust)

  # Tier 2 = stable-only, NOT in the final consensus list
  tier2 <- final |>
    dplyr::filter(selected == TRUE, keep_final == FALSE) |>
    dplyr::arrange(collection, dplyr::desc(stability), p.adjust)

  # Combined with tier label (for plotting / reporting)
  combined <- final |>
    dplyr::mutate(tier = dplyr::case_when(
      keep_final == TRUE ~ "Tier 1 (stable + consensus)",
      selected == TRUE   ~ "Tier 2 (stable only)",
      TRUE               ~ "Not kept"
    )) |>
    dplyr::arrange(
      match(tier, c("Tier 1 (stable + consensus)", "Tier 2 (stable only)", "Not kept")),
      collection, dplyr::desc(stability), p.adjust
    )

  # ---- Outputs ---------------------------------------------------------------
  OUT_TIER1 <- file.path(paths$step07, "tier1_final.csv")
  OUT_TIER2 <- file.path(paths$step07, "tier2_stable_only.csv")
  OUT_ALL   <- file.path(paths$step07, "tiered_all_candidates.csv")

  readr::write_csv(tier1,    OUT_TIER1)
  readr::write_csv(tier2,    OUT_TIER2)
  readr::write_csv(combined, OUT_ALL)

  cat("\n[Tiers written]\n")
  cat("Tier 1 (stable + consensus):  ", nrow(tier1),    " pathways -> ", OUT_TIER1, "\n", sep = "")
  cat("Tier 2 (stable only):         ", nrow(tier2),    " pathways -> ", OUT_TIER2, "\n", sep = "")
  cat("All candidates w/ tier label: ", nrow(combined), " rows     -> ", OUT_ALL,   "\n\n", sep = "")

  invisible(list(
    tier1_file       = OUT_TIER1,
    tier2_file       = OUT_TIER2,
    final_table_file = OUT_ALL
  ))
}
