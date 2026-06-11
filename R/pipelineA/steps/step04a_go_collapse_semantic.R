# R/pipelineA/steps/step04a_go_collapse_semantic.R
# Pipeline A — Step 04a: Collapse redundant GO-BP terms via GO semantic similarity
#
# Reads:
#   paths$step03_go_bp/ora_go_bp_filtered.csv   (written by step03b_ora_go_bp)
#
# Writes into paths$step04/semantic/:
#   reps.csv        representative terms kept at the configured threshold
#   mapping.csv     dropped terms -> their assigned representative
#   run_note.txt    summary (counts + DE gene coverage)
#
# Returns:
#   reps_file      path to reps.csv
#   mapping_file   path to mapping.csv
#   run_note_file  path to run_note.txt

suppressPackageStartupMessages({
  library(GOSemSim)
  library(org.Hs.eg.db)
})

step04a_go_collapse_semantic <- function(cfg, paths, ctx) {

  # ---- Parameters ------------------------------------------------------------
  tau     <- if (!is.null(cfg$thresholds$GO_SEMANTIC_TAU)) cfg$thresholds$GO_SEMANTIC_TAU else 0.60
  MEASURE <- "Wang"

  # ---- Inputs ----------------------------------------------------------------
  INFILE <- file.path(paths$step03_go_bp, "ora_go_bp_filtered.csv")
  stopifnot(file.exists(INFILE))

  tbl <- readr::read_csv(INFILE, show_col_types = FALSE)

  if (nrow(tbl) == 0) {
    stop("Input GO-BP ORA table is empty: ", INFILE)
  }

  req_cols <- c("ID", "Description", "p.adjust", "k", "K", "enr_ratio", "geneID")
  stopifnot(all(req_cols %in% names(tbl)))

  tbl <- tbl |>
    dplyr::transmute(
      ID          = as.character(ID),
      Description = as.character(Description),
      p.adjust    = as.numeric(p.adjust),
      k           = as.integer(k),
      K           = as.integer(K),
      enr_ratio   = as.numeric(enr_ratio),
      geneID      = as.character(geneID)
    ) |>
    dplyr::distinct(ID, .keep_all = TRUE) |>
    dplyr::arrange(p.adjust, dplyr::desc(enr_ratio), dplyr::desc(k))

  # ---- GO semantic data object -----------------------------------------------
  semBP <- GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = "BP")

  # ---- Pairwise GO-term semantic similarity matrix ---------------------------
  go_ids  <- tbl$ID
  sim_mat <- GOSemSim::mgoSim(go_ids, go_ids, semData = semBP,
                               measure = MEASURE, combine = NULL) |>
    as.matrix()
  diag(sim_mat) <- 1
  rownames(sim_mat) <- colnames(sim_mat) <- go_ids

  # ---- Precompute DE hit sets for coverage summaries -------------------------
  hits_list <- tbl |>
    dplyr::mutate(hits = strsplit(geneID, "/")) |>
    dplyr::select(ID, hits)    # dplyr::select to avoid AnnotationDbi::select

  all_hits                    <- unique(unlist(hits_list$hits))
  total_hit_genes_precollapse <- length(all_hits)

  # ---- Greedy collapse given a similarity cutoff tau -------------------------
  collapse_semantic <- function(tbl, sim, cutoff) {
    ord <- order(tbl$p.adjust, -tbl$enr_ratio, -tbl$k)
    dropped <- rep(FALSE, nrow(tbl)); names(dropped) <- tbl$ID
    keep_ids     <- character(0)
    mapping_rows <- list()

    for (i in ord) {
      id <- tbl$ID[i]
      if (dropped[id]) next
      keep_ids <- c(keep_ids, id)

      sims      <- sim[id, ]
      close_ids <- names(sims)[sims >= cutoff & names(sims) != id & !dropped[names(sims)]]
      if (length(close_ids)) {
        dropped[close_ids] <- TRUE
        mapping_rows[[length(mapping_rows) + 1]] <-
          tibble::tibble(
            original_term_id = close_ids,
            kept_term_id     = id,
            pair_similarity  = as.numeric(sims[close_ids])
          )
      }
    }

    reps <- tbl |>
      dplyr::filter(ID %in% keep_ids) |>
      dplyr::arrange(p.adjust, dplyr::desc(enr_ratio), dplyr::desc(k))

    mapping <- if (length(mapping_rows)) dplyr::bind_rows(mapping_rows) else
      tibble::tibble(
        original_term_id = character(),
        kept_term_id     = character(),
        pair_similarity  = numeric()
      )

    list(reps = reps, mapping = mapping)
  }

  # ---- Coverage of DE hits retained by the representatives -------------------
  coverage_pct <- function(reps_ids, hits_tbl) {
    kept_hits <- hits_tbl |>
      dplyr::filter(ID %in% reps_ids) |>
      dplyr::pull(hits) |>
      unlist() |>
      unique()
    100 * length(kept_hits) / total_hit_genes_precollapse
  }

  # ---- Run collapse and write outputs ----------------------------------------
  OUT_DIR <- file.path(paths$step04, "semantic")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  coll    <- collapse_semantic(tbl, sim_mat, cutoff = tau)
  reps    <- coll$reps
  mapping <- coll$mapping
  cov     <- coverage_pct(reps$ID, hits_list)

  reps_path    <- file.path(OUT_DIR, "reps.csv")
  mapping_path <- file.path(OUT_DIR, "mapping.csv")
  note_path    <- file.path(OUT_DIR, "run_note.txt")

  readr::write_csv(reps,    reps_path)
  readr::write_csv(mapping, mapping_path)

  note <- sprintf(paste0(
    "[GO collapse — semantic]\n",
    "measure = %s, threshold tau = %.2f\n",
    "input_terms = %d\n",
    "kept_reps   = %d\n",
    "DE coverage = %.1f%%\n"),
    MEASURE, tau, nrow(tbl), nrow(reps), cov
  )

  writeLines(note, note_path)
  cat(note, "\n")

  invisible(list(
    reps_file     = reps_path,
    mapping_file  = mapping_path,
    run_note_file = note_path
  ))
}
