# R/pipelineB/steps/step04_go_collapse_semantic_agreement.R
# Pipeline B — Step 04: GO:BP semantic redundancy reduction on agreement-filtered GO list.
#
# What this does:
# - Takes the GO terms that passed the CAMERA ∩ fgsea agreement filter (Step03).
# - Collapses semantically similar GO:BP terms using Wang semantic similarity (GOSemSim).
# - Collapse is done separately for Up and Down directions.
# - Runs a single tau (default 0.60; override via cfg$pipelineB$go_semantic_tau).
# - Representative selection: lowest CAMERA FDR → lowest fgsea FDR → shorter name → stronger stats.
#
# Reads (from current run):
#   paths$step03/go_bp_input_agreement.csv
#
# Writes into paths$step04  (results/pipelineB/<run_id>/04_go_collapse/):
#   reps.csv
#   mapping.csv
#   run_note.txt
#
# Returns (invisibly):
#   reps_file      path to reps.csv
#   mapping_file   path to mapping.csv

suppressPackageStartupMessages({
  library(GOSemSim)
  library(org.Hs.eg.db)
})

step04_go_collapse_semantic_agreement <- function(cfg, paths, ctx) {

  # ---- Parameters -----------------------------------------------------------
  INFILE       <- file.path(paths$step03, "go_bp_input_agreement.csv")
  OUT_DIR      <- paths$step04
  REPS_FILE    <- file.path(OUT_DIR, "reps.csv")
  MAPPING_FILE <- file.path(OUT_DIR, "mapping.csv")

  TAU     <- if (!is.null(cfg$pipelineB$go_semantic_tau)) as.numeric(cfg$pipelineB$go_semantic_tau) else 0.60
  MEASURE <- "Wang"

  # Column name candidates (robust detection)
  DIR_CAND     <- c("direction","Direction","dir")
  ID_CAND      <- c("ID","id","go_id","GO","go")
  DESC_CAND    <- c("Description","description","term","name","title")
  CAMFDR_CAND  <- c("cam_FDR","camera_FDR","cam_fdr","FDR_camera","camera_fdr","cam.padj","cam_padj")
  FGFDR_CAND   <- c("fg_FDR","fgsea_FDR","fg_padj","padj_fgsea","fgsea_padj","fg.padj","fg_padj")
  CAMSTAT_CAND <- c("cam_stat","camera_stat","cam_Stat","stat_camera")
  FGNES_CAND   <- c("NES","fg_NES","fgsea_NES","fg_nes")

  stopifnot(file.exists(INFILE))

  # ---- Helpers --------------------------------------------------------------
  getcol <- function(df, candidates) {
    nm <- names(df); nml <- tolower(nm)
    hit <- nm[match(tolower(candidates), nml, nomatch = 0)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }

  stop_if_missing <- function(x, label) {
    if (is.na(x) || !nzchar(x)) stop("Missing required column for: ", label)
  }

  safe_num <- function(x) suppressWarnings(as.numeric(x))

  # Greedy collapse: best-first keep; drop any not-yet-dropped terms with similarity >= τ
  collapse_semantic <- function(tbl, sim_mat, tau) {
    ids     <- tbl$ID
    dropped <- stats::setNames(rep(FALSE, length(ids)), ids)
    keep_ids     <- character(0)
    mapping_rows <- list()

    for (id in ids) {
      if (dropped[[id]]) next

      keep_ids <- c(keep_ids, id)

      sims      <- sim_mat[id, ]
      close_ids <- names(sims)[sims >= tau & names(sims) != id & !dropped[names(sims)]]

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

    reps <- tbl |> dplyr::filter(ID %in% keep_ids)

    mapping <- if (length(mapping_rows)) dplyr::bind_rows(mapping_rows) else
      tibble::tibble(original_term_id = character(),
                     kept_term_id     = character(),
                     pair_similarity  = numeric())

    list(reps = reps, mapping = mapping)
  }

  # ---- 0) Read input --------------------------------------------------------
  raw <- readr::read_csv(INFILE, show_col_types = FALSE)

  id_col       <- getcol(raw, ID_CAND);    stop_if_missing(id_col,     "GO term ID")
  dir_col      <- getcol(raw, DIR_CAND);   stop_if_missing(dir_col,    "direction (Up/Down)")
  desc_col     <- getcol(raw, DESC_CAND)   # optional
  cam_fdr_col  <- getcol(raw, CAMFDR_CAND);  stop_if_missing(cam_fdr_col, "CAMERA FDR")
  fg_fdr_col   <- getcol(raw, FGFDR_CAND);   stop_if_missing(fg_fdr_col,  "fgsea padj/FDR")
  cam_stat_col <- getcol(raw, CAMSTAT_CAND)  # optional
  fg_nes_col   <- getcol(raw, FGNES_CAND)    # optional

  tbl <- raw |>
    dplyr::transmute(
      ID          = as.character(.data[[id_col]]),
      direction   = as.character(.data[[dir_col]]),
      Description = if (!is.na(desc_col)) as.character(.data[[desc_col]]) else NA_character_,
      cam_FDR     = safe_num(.data[[cam_fdr_col]]),
      fg_FDR      = safe_num(.data[[fg_fdr_col]]),
      cam_stat    = if (!is.na(cam_stat_col)) safe_num(.data[[cam_stat_col]]) else NA_real_,
      fg_NES      = if (!is.na(fg_nes_col))  safe_num(.data[[fg_nes_col]])  else NA_real_
    ) |>
    dplyr::filter(!is.na(ID), ID != "", !is.na(cam_FDR), !is.na(fg_FDR)) |>
    dplyr::distinct(ID, direction, .keep_all = TRUE)

  # Normalize direction labels
  tbl <- tbl |>
    dplyr::mutate(direction = dplyr::case_when(
      stringr::str_to_lower(direction) %in% c("up","upregulated","up_reg","positive","+")      ~ "Up",
      stringr::str_to_lower(direction) %in% c("down","downregulated","down_reg","negative","-") ~ "Down",
      TRUE ~ direction
    ))

  # ---- Early exit: 0 input terms --------------------------------------------
  if (nrow(tbl) == 0) {
    cat("GO semantic collapse: 0 input terms — writing empty outputs and returning.\n")
    empty_reps <- tibble::tibble(
      ID=character(), direction=character(), Description=character(),
      cam_FDR=double(), fg_FDR=double(), cam_stat=double(), fg_NES=double()
    )
    empty_map <- tibble::tibble(
      original_term_id=character(), kept_term_id=character(),
      pair_similarity=double(), direction=character()
    )
    readr::write_csv(empty_reps, REPS_FILE)
    readr::write_csv(empty_map,  MAPPING_FILE)
    writeLines("[GO semantic collapse: 0 input terms — skipped]", file.path(OUT_DIR, "run_note.txt"))
    return(invisible(list(reps_file = REPS_FILE, mapping_file = MAPPING_FILE)))
  }

  dirs <- sort(unique(tbl$direction))
  if (!all(c("Up","Down") %in% dirs)) {
    message("Direction levels found: ", paste(dirs, collapse = ", "))
    message("Proceeding, but ideally this file contains 'Up' and/or 'Down'.")
  }

  # ---- 1) Compute semantic similarity matrix --------------------------------
  go_ids <- unique(tbl$ID)

  semBP <- GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = "BP")

  sim_mat <- GOSemSim::mgoSim(
    go_ids, go_ids,
    semData = semBP,
    measure = MEASURE,
    combine = NULL
  ) |> as.matrix()

  rownames(sim_mat) <- colnames(sim_mat) <- go_ids
  diag(sim_mat) <- 1
  sim_mat[is.na(sim_mat)] <- 0

  # ---- 2) Collapse per direction --------------------------------------------
  reps_all    <- list()
  mapping_all <- list()

  for (dir in sort(unique(tbl$direction))) {
    sub <- tbl |> dplyr::filter(direction == dir)

    if (nrow(sub) == 0) next

    # Representative ranking rule (best-first):
    # 1) lowest CAMERA FDR  2) lowest fgsea FDR/padj
    # 3) shorter term name  4) stronger stats (|NES|, |cam_stat|)
    sub_ranked <- sub |>
      dplyr::arrange(
        cam_FDR,
        fg_FDR,
        ifelse(is.na(Description), Inf, nchar(Description)),
        dplyr::desc(abs(fg_NES)),
        dplyr::desc(abs(cam_stat))
      )

    ids_dir <- sub_ranked$ID
    sim_dir <- sim_mat[ids_dir, ids_dir, drop = FALSE]

    coll <- collapse_semantic(sub_ranked, sim_dir, TAU)

    reps_all[[dir]]    <- coll$reps    |> dplyr::mutate(direction = dir)
    mapping_all[[dir]] <- coll$mapping |> dplyr::mutate(direction = dir)
  }

  reps    <- dplyr::bind_rows(reps_all)
  mapping <- dplyr::bind_rows(mapping_all)

  # Fill any missing GO descriptions in representatives
  if (!exists("add_go_term_name", mode = "function"))
    source("R/utils/go_term_names.R")
  reps <- add_go_term_name(reps, id_col = "ID", desc_col = "Description")

  readr::write_csv(reps,    REPS_FILE)
  readr::write_csv(mapping, MAPPING_FILE)

  # ---- 3) Run note ----------------------------------------------------------
  in_total   <- nrow(tbl)
  kept_total <- nrow(reps)
  in_up      <- sum(tbl$direction == "Up")
  in_down    <- sum(tbl$direction == "Down")
  kept_up    <- sum(reps$direction == "Up")
  kept_down  <- sum(reps$direction == "Down")

  note <- glue::glue("
[Pipeline B \u2014 GO semantic collapse (agreement input)]
Measure: {MEASURE}
Threshold \u03c4: {sprintf('%.2f', TAU)}
Collapse scope: within-direction (Up and Down separately)

Input GO terms total: {in_total}   (Up={in_up}, Down={in_down})
Representatives kept: {kept_total} (Up={kept_up}, Down={kept_down})

Representative rule (best-first):
1) lowest CAMERA FDR
2) then lowest fgsea padj/FDR
3) then shorter term name (if available)
4) then stronger stats (|NES|, |cam_stat|) as mild tie-break

Files:
- reps.csv (representatives)
- mapping.csv (collapsed -> representative)
")

  writeLines(note, file.path(OUT_DIR, "run_note.txt"))
  cat(note, "\nDone.\n")

  invisible(list(
    reps_file    = REPS_FILE,
    mapping_file = MAPPING_FILE
  ))
}
