# R/pipelineA/steps/step03a_ora_wp.R
# Pipeline A — Step 03a: ORA on WikiPathways with data-driven filters
#
# Reads:
#   paths$step01/universe_genes.txt        (written by step01_make_universe)
#   paths$step02/de_genes_rawp05.txt       (written by step02_make_de)
#   cfg$dataset$pathway2gene_file
#
# Writes into paths$step03_wp:
#   ora_wikipathways_all.csv
#   ora_wikipathways_filtered.csv
#   run_note.txt

suppressPackageStartupMessages({
  library(clusterProfiler)
})

step03a_ora_wp <- function(cfg, paths, ctx) {

  # ---- Inputs ----------------------------------------------------------------
  UNIVERSE_FILE <- file.path(paths$step01, "universe_genes.txt")
  DE_FILE       <- file.path(paths$step02, "de_genes_rawp05.txt")
  P2G_FILE      <- cfg$dataset$pathway2gene_file

  stopifnot(file.exists(UNIVERSE_FILE), file.exists(DE_FILE), file.exists(P2G_FILE))

  universe <- readr::read_lines(UNIVERSE_FILE)
  universe <- unique(stringr::str_trim(universe))
  universe <- universe[universe != ""]

  de <- readr::read_lines(DE_FILE)
  de <- unique(stringr::str_trim(de))
  de <- de[de != ""]

  p2g <- readr::read_csv(P2G_FILE, show_col_types = FALSE)

  if (length(universe) == 0) stop("Universe is empty — run step01 first.")

  # ---- Outputs ---------------------------------------------------------------
  OUT_ALL  <- file.path(paths$step03_wp, "ora_wikipathways_all.csv")
  OUT_FILT <- file.path(paths$step03_wp, "ora_wikipathways_filtered.csv")
  OUT_NOTE <- file.path(paths$step03_wp, "run_note.txt")

  if (length(de) == 0) {
    msg <- "DE gene list is empty (step02 produced no genes). Downstream steps will be skipped."
    writeLines(paste0("[ORA run note]\n", msg), OUT_NOTE)
    readr::write_csv(tibble::tibble(), OUT_ALL)
    readr::write_csv(tibble::tibble(), OUT_FILT)
    cat(msg, "\n")
    return(invisible(list(ora_all_file = OUT_ALL, ora_filtered_file = OUT_FILT)))
  }

  # ---- Build term maps -------------------------------------------------------
  term2gene <- p2g |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") |>
    dplyr::select(wpid, SYMBOL) |>
    dplyr::distinct()

  term2name <- p2g |>
    dplyr::filter(!is.na(name), name != "") |>
    dplyr::select(wpid, name) |>
    dplyr::distinct()

  # ---- Data-driven size limits (K range) -------------------------------------
  N <- length(universe)
  n <- length(de)
  cat("Universe N:", N, " | DE n:", n, "\n")

  de_fraction <- n / N
  ora_mode <- if (!is.null(cfg$pipelineA$ora_min_enrich_ratio_mode) &&
                  nzchar(cfg$pipelineA$ora_min_enrich_ratio_mode))
    cfg$pipelineA$ora_min_enrich_ratio_mode else "adaptive"
  r_min <- if (ora_mode == "adaptive") {
    if      (de_fraction <= 0.10) 1.5
    else if (de_fraction <= 0.30) 1.2
    else                          1.0
  } else {
    fixed_r <- suppressWarnings(as.numeric(cfg$pipelineA$ora_min_enrich_ratio))
    if (!is.na(fixed_r)) fixed_r else 1.5
  }
  cat(sprintf("ORA enrich-ratio threshold: t=%.1f (mode=%s, de_fraction=%.3f)\n",
              r_min, ora_mode, de_fraction))

  K_tbl <- term2gene |>
    dplyr::filter(SYMBOL %in% universe) |>
    dplyr::count(wpid, name = "K") |>
    dplyr::inner_join(term2name, by = "wpid")

  if (nrow(K_tbl) == 0) stop("No pathways overlap your universe — check IDs / mapping.")

  K_q   <- quantile(K_tbl$K, probs = c(.10, .50, .90, .95), names = FALSE, type = 7)
  K_min <- max(15L, as.integer(K_q[1]))
  K_max <- min(500L, as.integer(K_q[4]))

  if (K_min > K_max) {
    K_min <- max(3L, K_max)
    cat(sprintf("  Warning: K_min > K_max after percentile computation; K_min clamped to %d\n", K_min))
  }

  cat("Pathway size K percentiles (10/50/90/95%):",
      paste(K_q, collapse = ", "), "\n")
  cat("Using K_min =", K_min, "and K_max =", K_max, "\n")

  # ---- ORA (Fisher/hypergeometric) -------------------------------------------
  ora <- enricher(
    gene          = de,
    universe      = universe,
    TERM2GENE     = term2gene,
    TERM2NAME     = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )

  if (is.null(ora) || nrow(ora@result) == 0) {
    msg <- "No pathways returned by ORA. Try relaxing DE rule or checking IDs/mapping."
    writeLines(paste0("[ORA run note]\n", msg), OUT_NOTE)
    readr::write_csv(tibble::tibble(), OUT_ALL)
    readr::write_csv(tibble::tibble(), OUT_FILT)
    cat(msg, "\n")
    return(invisible(list(ora_all_file = OUT_ALL, ora_filtered_file = OUT_FILT)))
  }

  # ---- Unpack ratios, add E[k] and enrichment ratio -------------------------
  res <- ora@result |>
    dplyr::arrange(p.adjust) |>
    tidyr::separate_wider_delim(GeneRatio, "/", names = c("k", "n_in_table"),
                                cols_remove = FALSE) |>
    tidyr::separate_wider_delim(BgRatio,   "/", names = c("K", "N_in_table"),
                                cols_remove = FALSE) |>
    dplyr::mutate(dplyr::across(c(k, n_in_table, K, N_in_table), as.numeric),
                  Ek        = (n * K) / N,
                  enr_ratio = k / pmax(Ek, .Machine$double.eps)) |>
    dplyr::relocate(ID, Description, k, K, GeneRatio, BgRatio, Ek, enr_ratio,
                    pvalue, p.adjust, Count)

  readr::write_csv(res, OUT_ALL)

  # ---- Practical filters -----------------------------------------------------
  q_fdr <- 0.05
  k_min <- 3L

  res_filt <- res |>
    dplyr::filter(p.adjust < q_fdr,
                  dplyr::between(K, K_min, K_max),
                  k >= k_min,
                  enr_ratio >= r_min)

  readr::write_csv(res_filt, OUT_FILT)

  # ---- Run note --------------------------------------------------------------
  de_policy_summary <- tryCatch({
    dp <- yaml::read_yaml(file.path(paths$step02, "de_policy.yml"))
    sprintf("p_col=%s, p_cutoff=%.3f, lfc_min=%s, steps=%s",
            dp$p_col_used, dp$final_p_cutoff, dp$final_lfc_min,
            if (length(dp$tightening_steps) == 0) "none"
            else paste(dp$tightening_steps, collapse = "|"))
  }, error = function(e) "de_policy.yml not found")

  note <- sprintf(paste0(
    "[ORA run note]\n",
    "N = %d (universe), n = %d (DE), de_fraction = %.3f\n",
    "DE policy: %s\n",
    "K_min..K_max = [%d, %d] (data-driven from K distribution)\n",
    "Filters: FDR < %.2f, k >= %d, k/E[k] >= %.1f (mode=%s)\n",
    "Kept %d / %d pathways after filters.\n",
    "Files: %s, %s\n"),
    N, n, de_fraction, de_policy_summary,
    K_min, K_max, q_fdr, k_min, r_min, ora_mode,
    nrow(res_filt), nrow(res),
    OUT_ALL, OUT_FILT
  )

  writeLines(note, OUT_NOTE)
  cat("\n", note, sep = "")

  invisible(list(
    ora_all_file      = OUT_ALL,
    ora_filtered_file = OUT_FILT
  ))
}
