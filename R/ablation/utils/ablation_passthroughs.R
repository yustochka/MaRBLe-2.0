# R/ablation/utils/ablation_passthroughs.R
# Passthrough functions that replace pipeline steps when a component is
# disabled in the ablation config.  Each passthrough:
#   - writes the same files a normal step would write (to keep folder structure identical)
#   - uses the same output schema so downstream steps need no changes
#   - writes a BYPASSED.txt marker so runs are clearly labeled

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

# ---------------------------------------------------------------------------
# Helper: detect a column by name (case-insensitive, first match wins)
# ---------------------------------------------------------------------------
.getcol_pt <- function(df, candidates) {
  nm  <- names(df)
  nml <- tolower(nm)
  hit <- nm[match(tolower(candidates), nml, nomatch = 0)]
  if (length(hit) == 0) NA_character_ else hit[1]
}


# ---------------------------------------------------------------------------
# 1. CAMERA-only step03 passthrough
#    Used when use_fgsea = false.
#    Reads step01 CAMERA output, keeps pathways at FDR <= 0.05, writes the
#    same wp_input_agreement.csv / go_bp_input_agreement.csv that normal
#    step03 produces — but fg_FDR = NA, fgsea_NES = NA.
#
# NOTE: Because fg_FDR = NA, the normal step04 function would drop all GO
# rows (it requires !is.na(fg_FDR)).  The ablation runner therefore always
# uses passthrough_step04_no_semantic when use_fgsea = false.
# ---------------------------------------------------------------------------

passthrough_step03_camera_only <- function(cfg, paths) {
  if (!exists("add_go_term_name", mode = "function"))
    source("R/utils/go_term_names.R")

  CAM_WP_FILE <- file.path(paths$step01, "camera_wikipathways.csv")
  CAM_GO_FILE <- file.path(paths$step01, "camera_go_bp.csv")
  FDR_CUT     <- 0.05

  stopifnot(file.exists(CAM_WP_FILE), file.exists(CAM_GO_FILE))

  # Build agreement-style tibble from a single CAMERA result file
  make_camera_only <- function(cam_df, label) {
    dir_col  <- .getcol_pt(cam_df, c("Direction", "direction"))
    fdr_col  <- .getcol_pt(cam_df, c("FDR", "adj.P.Val"))
    desc_col <- .getcol_pt(cam_df, c("Description", "description", "name"))

    if (is.na(dir_col))
      stop("CAMERA output missing Direction column for ", label)
    if (is.na(fdr_col))
      stop("CAMERA output missing FDR column for ", label)

    cam_df |>
      dplyr::transmute(
        ID          = as.character(ID),
        Description = if (!is.na(desc_col)) as.character(.data[[desc_col]]) else NA_character_,
        direction   = dplyr::case_when(
          tolower(.data[[dir_col]]) == "up"   ~ "Up",
          tolower(.data[[dir_col]]) == "down" ~ "Down",
          TRUE                                ~ NA_character_
        ),
        cam_FDR     = as.numeric(.data[[fdr_col]]),
        fg_FDR      = NA_real_,    # fgsea bypassed
        fgsea_NES   = NA_real_,    # fgsea bypassed
        agreement_source           = "camera_only",
        agreement_mode_used        = "camera_only",
        agreement_relaxed_fdr_used = NA_real_,
        agreement_topn_used        = NA_integer_
      ) |>
      dplyr::filter(!is.na(direction), cam_FDR <= FDR_CUT) |>
      dplyr::arrange(cam_FDR)
  }

  cam_wp <- readr::read_csv(CAM_WP_FILE, show_col_types = FALSE)
  cam_go <- readr::read_csv(CAM_GO_FILE, show_col_types = FALSE)

  wp_out <- make_camera_only(cam_wp, "WikiPathways")
  go_out <- make_camera_only(cam_go, "GO:BP")

  # Fill GO descriptions from GO.db
  go_out <- add_go_term_name(go_out, id_col = "ID", desc_col = "Description")

  wp_path <- file.path(paths$step03, "wp_input_agreement.csv")
  go_path <- file.path(paths$step03, "go_bp_input_agreement.csv")

  readr::write_csv(wp_out, wp_path)
  readr::write_csv(go_out, go_path)

  note <- paste0(
    "[CAMERA-only step03 passthrough — fgsea bypassed]\n",
    "Kept: CAMERA FDR <= ", FDR_CUT, ", same direction.\n",
    "fg_FDR = NA, fgsea_NES = NA for all rows.\n",
    "WP kept: ", nrow(wp_out), "\n",
    "GO kept: ", nrow(go_out), "\n"
  )
  writeLines(note, file.path(paths$step03, "run_note.txt"))
  writeLines("[Step02 skipped: use_fgsea = false]",
             file.path(paths$step02, "BYPASSED.txt"))
  cat(note)

  invisible(list(wp_agree_file = wp_path, go_agree_file = go_path))
}


# ---------------------------------------------------------------------------
# 2. No-semantic step04 passthrough
#    Used when use_semantic_collapse = false (or when use_fgsea = false,
#    because step04 requires fg_FDR non-NA).
#    Copies step03 go_bp_input_agreement.csv to step04/reps.csv with column
#    renaming to match the step04 reps.csv schema.  No GO terms are removed.
# ---------------------------------------------------------------------------

passthrough_step04_no_semantic <- function(paths, reason = "use_semantic_collapse=false") {
  GO_IN_FILE   <- file.path(paths$step03, "go_bp_input_agreement.csv")
  REPS_FILE    <- file.path(paths$step04, "reps.csv")
  MAPPING_FILE <- file.path(paths$step04, "mapping.csv")

  stopifnot(file.exists(GO_IN_FILE))

  raw <- readr::read_csv(GO_IN_FILE, show_col_types = FALSE)

  # step03 GO columns: ID, Description, direction, cam_FDR, fg_FDR, fgsea_NES, ...
  # step04 reps.csv:   ID, direction, Description, cam_FDR, fg_FDR, cam_stat, fg_NES
  reps <- raw |>
    dplyr::transmute(
      ID          = as.character(ID),
      direction   = as.character(direction),
      Description = as.character(Description),
      cam_FDR     = as.numeric(cam_FDR),
      fg_FDR      = as.numeric(fg_FDR),
      cam_stat    = NA_real_,
      fg_NES      = if ("fgsea_NES" %in% names(raw)) as.numeric(fgsea_NES) else NA_real_
    )

  empty_mapping <- tibble::tibble(
    original_term_id = character(),
    kept_term_id     = character(),
    pair_similarity  = numeric(),
    direction        = character()
  )

  readr::write_csv(reps,          REPS_FILE)
  readr::write_csv(empty_mapping, MAPPING_FILE)

  note <- paste0(
    "[No-semantic step04 passthrough — ", reason, "]\n",
    "GO semantic collapse bypassed.\n",
    "All agreement-filtered GO terms passed through unchanged.\n",
    "GO terms in reps.csv: ", nrow(reps), "\n"
  )
  writeLines(note, file.path(paths$step04, "run_note.txt"))
  writeLines(paste0("[Step04 bypassed: ", reason, "]"),
             file.path(paths$step04, "BYPASSED.txt"))
  cat(note)

  invisible(list(reps_file = REPS_FILE, mapping_file = MAPPING_FILE))
}


# ---------------------------------------------------------------------------
# 3. No-overlap step05 passthrough
#    Used when use_overlap_clustering = false.
#    Merges step03 WP + step04 GO (already semantically collapsed or not),
#    computes set_size from step01 CAMERA NGenes (universe-restricted count),
#    splits by direction, and writes representatives_up/down.csv with the
#    same column schema that step06 expects.
#    No Jaccard clustering is applied — every pathway proceeds to step06.
# ---------------------------------------------------------------------------

passthrough_step05_no_overlap <- function(cfg, paths) {
  WP_IN_FILE  <- file.path(paths$step03, "wp_input_agreement.csv")
  GO_IN_FILE  <- file.path(paths$step04, "reps.csv")
  CAM_WP_FILE <- file.path(paths$step01, "camera_wikipathways.csv")
  CAM_GO_FILE <- file.path(paths$step01, "camera_go_bp.csv")

  stopifnot(
    file.exists(WP_IN_FILE), file.exists(GO_IN_FILE),
    file.exists(CAM_WP_FILE), file.exists(CAM_GO_FILE)
  )

  MIN_SET_SIZE <- if (!is.null(cfg$pipelineB$min_set_size))
    as.integer(cfg$pipelineB$min_set_size) else 10L
  MAX_SET_SIZE <- if (!is.null(cfg$pipelineB$max_set_size))
    as.integer(cfg$pipelineB$max_set_size) else 500L

  # ---- Set sizes from CAMERA NGenes (universe-restricted pathway size) ------
  # All terms in step03/step04 passed CAMERA's set-size filter at step01,
  # so NGenes is the authoritative universe-restricted size.
  cam_wp_sizes <- readr::read_csv(CAM_WP_FILE, show_col_types = FALSE) |>
    dplyr::transmute(ID = as.character(ID), set_size = as.integer(NGenes))

  cam_go_sizes <- readr::read_csv(CAM_GO_FILE, show_col_types = FALSE) |>
    dplyr::transmute(ID = as.character(ID), set_size = as.integer(NGenes))

  # ---- Build standardized representative schema for step06 -----------------
  # step06 need_cols: ID, Description, direction, collection,
  #                   set_size, camera_FDR, fgsea_padj, agreement_q, absNES
  # step05 actual output also has: cluster, fgsea_NES

  make_reps <- function(in_df, sizes_df, collection_label, nes_col_name) {
    in_df |>
      dplyr::left_join(sizes_df, by = "ID") |>
      dplyr::filter(!is.na(set_size),
                    set_size >= MIN_SET_SIZE, set_size <= MAX_SET_SIZE) |>
      dplyr::transmute(
        cluster     = NA_integer_,
        ID          = as.character(ID),
        Description = as.character(Description),
        direction   = stringr::str_to_title(as.character(direction)),
        camera_FDR  = as.numeric(cam_FDR),
        fgsea_padj  = as.numeric(fg_FDR),
        fgsea_NES   = if (nes_col_name %in% names(in_df))
                        as.numeric(.data[[nes_col_name]]) else NA_real_,
        collection  = collection_label,
        set_size    = as.integer(set_size),
        agreement_q = pmax(camera_FDR, fgsea_padj, na.rm = TRUE),
        absNES      = ifelse(!is.na(fgsea_NES), abs(fgsea_NES), NA_real_)
      )
  }

  wp_in <- readr::read_csv(WP_IN_FILE, show_col_types = FALSE)
  go_in <- readr::read_csv(GO_IN_FILE, show_col_types = FALSE)

  # WP: fgsea_NES column name is "fgsea_NES" (from step03)
  # GO: NES column name is "fg_NES" (from step04) or "fgsea_NES" (from passthrough step04)
  go_nes_col <- if ("fg_NES" %in% names(go_in)) "fg_NES" else "fgsea_NES"

  wp_reps <- make_reps(wp_in, cam_wp_sizes, "WP", "fgsea_NES")
  go_reps <- make_reps(go_in, cam_go_sizes, "GO", go_nes_col)

  all_reps <- dplyr::bind_rows(wp_reps, go_reps) |>
    dplyr::arrange(direction, agreement_q, dplyr::desc(absNES))

  reps_up   <- dplyr::filter(all_reps, direction == "Up")
  reps_down <- dplyr::filter(all_reps, direction == "Down")

  # ---- Write representative files (step06 input) ---------------------------
  REPS_UP_FILE   <- file.path(paths$step05, "representatives_up.csv")
  REPS_DOWN_FILE <- file.path(paths$step05, "representatives_down.csv")
  readr::write_csv(reps_up,   REPS_UP_FILE)
  readr::write_csv(reps_down, REPS_DOWN_FILE)

  # Write empty edge / cluster files for structural consistency
  empty_edges <- tibble::tibble(
    from = character(), to = character(), jaccard = double(), k = integer()
  )
  empty_reps_schema <- tibble::tibble(
    ID=character(), Description=character(), direction=character(),
    camera_FDR=double(), fgsea_padj=double(), fgsea_NES=double(),
    collection=character(), set_size=integer(),
    agreement_q=double(), absNES=double(), cluster=integer()
  )
  for (fn in c("edges_up.csv", "edges_down.csv"))
    readr::write_csv(empty_edges,      file.path(paths$step05, fn))
  for (fn in c("clusters_up.csv", "clusters_down.csv"))
    readr::write_csv(empty_reps_schema, file.path(paths$step05, fn))

  note <- paste0(
    "[No-overlap step05 passthrough — use_overlap_clustering=false]\n",
    "Jaccard overlap clustering bypassed.\n",
    "All agreement-filtered terms (after optional semantic collapse) passed to step06.\n",
    "set_size sourced from step01 CAMERA NGenes (universe-restricted count).\n",
    "WP terms: ", nrow(wp_reps), "\n",
    "GO terms: ", nrow(go_reps), "\n",
    "Total UP:   ", nrow(reps_up),   "\n",
    "Total DOWN: ", nrow(reps_down), "\n"
  )
  writeLines(note, file.path(paths$step05, "run_note.txt"))
  writeLines("[Step05 bypassed: use_overlap_clustering=false]",
             file.path(paths$step05, "BYPASSED.txt"))
  cat(note)

  invisible(list(
    out_dir              = paths$step05,
    representatives_up   = REPS_UP_FILE,
    representatives_down = REPS_DOWN_FILE
  ))
}
