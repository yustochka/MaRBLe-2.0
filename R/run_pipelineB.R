# run_pipelineB.R
# Entry point for Pipeline B.

library(yaml)

source("R/utils/config_loader.R")
source("R/utils/run_context.R")
source("R/utils/paths_builder.R")
source("R/utils/validate_inputs.R")
source("R/utils/gene_id_convert.R")
source("R/utils/go_term_names.R")
source("R/pipelineB/steps/step01_camera_wp_go_bp.R")
source("R/pipelineB/steps/step02_fgsea_wp_go_bp.R")
source("R/pipelineB/steps/step03_prepare_inputs_agreement_first.R")
source("R/pipelineB/steps/step04_go_collapse_semantic_agreement.R")
source("R/pipelineB/steps/step05_overlap_cluster_B.R")
source("R/pipelineB/steps/step06_consensus_and_themes_A1A2.R")

cfg   <- load_config("config/default.yml")
ctx   <- create_run_context("pipelineB", cfg)
paths <- build_paths_b(ctx)

message("Pipeline B run initialized at ", ctx$run_dir)
message("Dataset: ", cfg$dataset$name)

validate_config_and_inputs(cfg, "pipelineB")

# ---------------------------------------------------------------------------
# Step-toggle helpers
# ---------------------------------------------------------------------------

enabled_steps <- if (!is.null(cfg$pipelineB$run_steps)) {
  as.integer(cfg$pipelineB$run_steps)
} else {
  c(1L, 2L, 3L, 4L, 5L, 6L)
}

step_enabled <- function(n) as.integer(n) %in% enabled_steps

message("Steps enabled: ", paste(enabled_steps, collapse = ", "))

# ---------------------------------------------------------------------------
# Manifest helper — writes manifest.yml into a step folder after it runs
# ---------------------------------------------------------------------------

write_manifest <- function(step_dir, step_num, inputs, outputs) {
  yaml::write_yaml(
    list(
      step      = step_num,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      inputs    = as.list(inputs),
      outputs   = as.list(outputs)
    ),
    file.path(step_dir, "manifest.yml")
  )
}

ran_steps <- integer(0)

# ---------------------------------------------------------------------------
# Step 01: CAMERA enrichment (WP + GO:BP)
# ---------------------------------------------------------------------------

if (step_enabled(1)) {
  result01 <- step01_camera_wp_go_bp(cfg, paths, ctx)
  message("Universe written to:       ", result01$universe_file)
  message("CAMERA WP written to:      ", result01$camera_wp_file)
  message("CAMERA GO:BP written to:   ", result01$camera_go_file)
  write_manifest(paths$step01, 1,
    inputs  = list(expression_file = cfg$dataset$expression_file,
                   metadata_file   = cfg$dataset$metadata_file),
    outputs = list(universe_file   = result01$universe_file,
                   camera_wp_file  = result01$camera_wp_file,
                   camera_go_file  = result01$camera_go_file)
  )
  message("--- Pipeline B Step01 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 1L)
}

# ---------------------------------------------------------------------------
# Step 02: fgsea sensitivity analysis (WP + GO:BP)
# ---------------------------------------------------------------------------

if (step_enabled(2)) {
  result02 <- step02_fgsea_wp_go_bp(cfg, paths, ctx)
  message("fgsea WP written to:       ", result02$fgsea_wp_file)
  message("fgsea GO:BP written to:    ", result02$fgsea_go_file)
  write_manifest(paths$step02, 2,
    inputs  = list(expression_file = cfg$dataset$expression_file,
                   metadata_file   = cfg$dataset$metadata_file,
                   universe_file   = file.path(paths$step01, "universe_genes.txt")),
    outputs = list(fgsea_wp_file   = result02$fgsea_wp_file,
                   fgsea_go_file   = result02$fgsea_go_file)
  )
  message("--- Pipeline B Step02 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 2L)
}

# ---------------------------------------------------------------------------
# Step 03: Agreement filter — keep CAMERA ∩ fgsea (same direction)
# ---------------------------------------------------------------------------

if (step_enabled(3)) {
  result03 <- step03_prepare_inputs_agreement_first(cfg, paths, ctx)
  message("WP agreement written to:    ", result03$wp_agree_file)
  message("GO:BP agreement written to: ", result03$go_agree_file)
  write_manifest(paths$step03, 3,
    inputs  = list(step01_dir = paths$step01,
                   step02_dir = paths$step02),
    outputs = list(wp_agree_file = result03$wp_agree_file,
                   go_agree_file = result03$go_agree_file)
  )
  message("--- Pipeline B Step03 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 3L)
}

# ---------------------------------------------------------------------------
# Step 04: GO semantic collapse on agreement-filtered list
# ---------------------------------------------------------------------------

if (step_enabled(4)) {
  result04 <- step04_go_collapse_semantic_agreement(cfg, paths, ctx)
  message("GO collapse reps written to:    ", result04$reps_file)
  message("GO collapse mapping written to: ", result04$mapping_file)
  write_manifest(paths$step04, 4,
    inputs  = list(step03_dir = paths$step03),
    outputs = list(reps_file    = result04$reps_file,
                   mapping_file = result04$mapping_file)
  )
  message("--- Pipeline B Step04 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 4L)
}

# ---------------------------------------------------------------------------
# Step 05: Overlap-based clustering (Jaccard graph, direction-aware)
# ---------------------------------------------------------------------------

if (step_enabled(5)) {
  result05 <- step05_overlap_cluster_B(cfg, paths, ctx)
  message("Overlap output dir:              ", paths$step05)
  message("Representatives Up written to:   ", result05$representatives_up)
  message("Representatives Down written to: ", result05$representatives_down)
  write_manifest(paths$step05, 5,
    inputs  = list(step01_dir = paths$step01,
                   step03_dir = paths$step03,
                   step04_dir = paths$step04),   # GO input now from step04/reps.csv
    outputs = list(out_dir              = paths$step05,
                   representatives_up   = result05$representatives_up,
                   representatives_down = result05$representatives_down)
  )
  message("--- Pipeline B Step05 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 5L)
}

# ---------------------------------------------------------------------------
# Step 06: Consensus + Themes
# ---------------------------------------------------------------------------

if (step_enabled(6)) {
  result06 <- step06_consensus_and_themes_A1A2(cfg, paths, ctx)
  message("Final UP written to:         ", result06$final_up_file)
  message("Final DOWN written to:       ", result06$final_down_file)
  # manifest intentionally skipped for step06
  message("--- Pipeline B Step06 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 6L)
}

# ---------------------------------------------------------------------------
# Run metadata
# ---------------------------------------------------------------------------

git_hash <- tryCatch(
  trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE)),
  error = function(e) NA_character_
)

run_meta <- list(
  pipeline   = "pipelineB",
  run_id     = ctx$run_id,
  timestamp  = ctx$timestamp,
  dataset    = cfg$dataset$name,
  ran_steps  = as.list(ran_steps),
  git_commit = if (length(git_hash) == 1 && nchar(git_hash) > 0) git_hash else NA_character_,
  run_dir    = ctx$run_dir
)

yaml::write_yaml(run_meta, file.path(ctx$run_dir, "run_meta.yml"))

# ---------------------------------------------------------------------------
# FINAL export — always runs after steps complete
# Applies final selection policy (consensus_first + top-N cap) to step06 full
# tables, then writes FINAL/FINAL.csv + FINAL/final_summary.txt
# ---------------------------------------------------------------------------

FINAL_DIR   <- file.path(ctx$run_dir, "FINAL")
dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)

STEP06_UP   <- file.path(paths$step06, "final_up.csv")
STEP06_DOWN <- file.path(paths$step06, "final_down.csv")

# Final selection policy parameters (with defaults)
FINAL_MODE  <- if (!is.null(cfg$pipelineB$final_mode))
  cfg$pipelineB$final_mode else "consensus_first"
FINAL_TOP_N <- if (!is.null(cfg$pipelineB$final_top_n_per_direction))
  as.integer(cfg$pipelineB$final_top_n_per_direction) else 5L

# Select top-N from a candidate set for one direction
apply_top_n <- function(tbl) {
  tbl |>
    dplyr::arrange(agreement_q, dplyr::desc(absNES), camera_FDR, fgsea_padj) |>
    dplyr::slice_head(n = FINAL_TOP_N)
}

if (file.exists(STEP06_UP) && file.exists(STEP06_DOWN)) {
  read_final_csv <- function(path) {
    readr::read_csv(path, show_col_types = FALSE,
      col_types = readr::cols(
        camera_FDR  = readr::col_double(),
        fgsea_padj  = readr::col_double(),
        fgsea_NES   = readr::col_double(),
        agreement_q = readr::col_double(),
        absNES      = readr::col_double(),
        set_size    = readr::col_integer(),
        cluster     = readr::col_character(),
        .default    = readr::col_character()
      ))
  }
  up   <- read_final_csv(STEP06_UP)
  down <- read_final_csv(STEP06_DOWN)

  # Step03 provenance (for summary + empty note)
  detect_s3_src <- function(fpath) {
    if (!file.exists(fpath)) return(list(rows = 0L, src = "strict_intersection"))
    tmp <- tryCatch(readr::read_csv(fpath, show_col_types = FALSE),
                   error = function(e) tibble::tibble())
    list(rows = nrow(tmp),
         src  = if ("agreement_source" %in% names(tmp) && nrow(tmp) > 0)
                  as.character(tmp$agreement_source[1]) else "strict_intersection")
  }
  s3_wp          <- detect_s3_src(file.path(paths$step03, "wp_input_agreement.csv"))
  s3_go          <- detect_s3_src(file.path(paths$step03, "go_bp_input_agreement.csv"))
  step03_wp_rows <- s3_wp$rows;  step03_wp_src <- s3_wp$src
  step03_go_rows <- s3_go$rows;  step03_go_src <- s3_go$src
  step03_total   <- step03_wp_rows + step03_go_rows

  # Candidate set counts (before top-N cap)
  n_cons_up       <- sum(up$in_consensus   == TRUE, na.rm = TRUE)
  n_solo_up       <- sum(up$solo_high_conf == TRUE, na.rm = TRUE)
  n_fallback_up      <- if ("final_reason" %in% names(up))
    sum(up$final_reason %in% c("fallback_ranked", "fallback_ranked_low_conf"), na.rm = TRUE) else 0L
  n_last_resort_up   <- if ("final_reason" %in% names(up))
    sum(up$final_reason == "fallback_topn_last_resort", na.rm = TRUE) else 0L
  n_cons_down        <- sum(down$in_consensus   == TRUE, na.rm = TRUE)
  n_solo_down        <- sum(down$solo_high_conf == TRUE, na.rm = TRUE)
  n_fallback_down    <- if ("final_reason" %in% names(down))
    sum(down$final_reason %in% c("fallback_ranked", "fallback_ranked_low_conf"), na.rm = TRUE) else 0L
  n_last_resort_down <- if ("final_reason" %in% names(down))
    sum(down$final_reason == "fallback_topn_last_resort", na.rm = TRUE) else 0L

  # consensus_first: consensus > solo > fallback_ranked (Option A)
  candidates_up <- if (FINAL_MODE == "consensus_first") {
    if (n_cons_up >= 1L) dplyr::filter(up, in_consensus   == TRUE)
    else if (n_solo_up >= 1L) dplyr::filter(up, solo_high_conf == TRUE)
    else dplyr::filter(up, keep_final == TRUE)
  } else {
    dplyr::filter(up, keep_final == TRUE)
  }

  candidates_down <- if (FINAL_MODE == "consensus_first") {
    if (n_cons_down >= 1L) dplyr::filter(down, in_consensus   == TRUE)
    else if (n_solo_down >= 1L) dplyr::filter(down, solo_high_conf == TRUE)
    else dplyr::filter(down, keep_final == TRUE)
  } else {
    dplyr::filter(down, keep_final == TRUE)
  }

  source_up <- if (n_cons_up >= 1L) "consensus"
    else if (n_solo_up >= 1L) "solo"
    else if (n_fallback_up >= 1L) "fallback_ranked"
    else if (n_last_resort_up >= 1L) "fallback_topn_last_resort"
    else "none"
  source_down <- if (n_cons_down >= 1L) "consensus"
    else if (n_solo_down >= 1L) "solo"
    else if (n_fallback_down >= 1L) "fallback_ranked"
    else if (n_last_resort_down >= 1L) "fallback_topn_last_resort"
    else "none"

  final_up   <- apply_top_n(candidates_up)
  final_down <- apply_top_n(candidates_down)
  final_tbl  <- dplyr::bind_rows(final_up, final_down)

  # ---- Fill missing Descriptions before writing ----------------------------
  # GO: look up TERM from GO.db using the GOID
  go_na_ids <- final_tbl |>
    dplyr::filter(is.na(Description), collection == "GO", startsWith(ID, "GO:")) |>
    dplyr::pull(ID)
  if (length(go_na_ids) > 0) {
    suppressPackageStartupMessages({ library(AnnotationDbi); library(GO.db) })
    go_raw <- suppressMessages(
      AnnotationDbi::select(GO.db, keys = go_na_ids, keytype = "GOID", columns = "TERM")
    )
    go_lookup <- stats::setNames(go_raw$TERM, go_raw$GOID)
    final_tbl <- final_tbl |>
      dplyr::mutate(
        Description = dplyr::if_else(
          is.na(Description) & collection == "GO" & startsWith(ID, "GO:"),
          go_lookup[ID],
          Description
        )
      )
    message("Filled GO descriptions: ", length(go_na_ids), " terms resolved via GO.db")
  }

  # WP: look up name from step01 CAMERA WP output (has ID + Description)
  wp_na_ids <- final_tbl |>
    dplyr::filter(is.na(Description), collection == "WP", startsWith(ID, "WP")) |>
    dplyr::pull(ID)
  if (length(wp_na_ids) > 0) {
    wp_cam_file <- file.path(paths$step01, "camera_wikipathways.csv")
    if (file.exists(wp_cam_file)) {
      wp_names <- readr::read_csv(wp_cam_file, show_col_types = FALSE) |>
        dplyr::filter(!is.na(Description)) |>
        dplyr::distinct(ID, .keep_all = TRUE)
      wp_lookup <- stats::setNames(wp_names$Description, wp_names$ID)
      final_tbl <- final_tbl |>
        dplyr::mutate(
          Description = dplyr::if_else(
            is.na(Description) & collection == "WP" & startsWith(ID, "WP"),
            wp_lookup[ID],
            Description
          )
        )
      message("Filled WP descriptions: ", length(wp_na_ids), " terms resolved via camera_wikipathways.csv")
    }
    still_na_wp <- final_tbl |>
      dplyr::filter(is.na(Description), collection == "WP") |>
      dplyr::pull(ID)
    if (length(still_na_wp) > 0) {
      warning("No Description found for WP IDs: ", paste(still_na_wp, collapse = ", "))
    }
  }
  # --------------------------------------------------------------------------

  readr::write_csv(final_tbl, file.path(FINAL_DIR, "FINAL.csv"))

  n_up    <- nrow(final_up)
  n_down  <- nrow(final_down)
  n_total <- nrow(final_tbl)

  agree_mode_cfg <- if (!is.null(cfg$pipelineB$agreement_mode) && nzchar(cfg$pipelineB$agreement_mode))
    cfg$pipelineB$agreement_mode else "intersection"

  empty_final_note <- if (n_total == 0) {
    if (step03_total == 0) {
      paste0(
        "\nFINAL empty.\n",
        "Reason: Step03 agreement (mode: '", agree_mode_cfg,
        "') returned 0 pathways under current thresholds.\n",
        "Tip:    Set agreement_mode: 'relax_fdr_then_topn' in pipelineB config to enable fallback.\n"
      )
    } else {
      paste0(
        "\nFINAL empty.\n",
        "Reason: Step03 kept ", step03_total, " pathways",
        " (WP=", step03_wp_rows, " [", step03_wp_src, "]",
        ", GO=", step03_go_rows, " [", step03_go_src, "]) but Step06 kept 0.\n",
        "Tip:    Lower solo_q_max/solo_abs_nes_min, or reduce optionA_quality_floor_q,\n",
        "        or set optionA_allow_low_conf: true in pipelineB config.\n"
      )
    }
  } else ""

  summary_text <- paste0(
    "Pipeline:                pipelineB\n",
    "Dataset:                 ", cfg$dataset$name, "\n",
    "Run ID:                  ", ctx$run_id, "\n",
    "Timestamp:               ", ctx$timestamp, "\n",
    "Steps run:               ", paste(sort(ran_steps), collapse = ", "), "\n",
    "\nStep03 agreement:\n",
    "  WP source:             ", step03_wp_src,  " (", step03_wp_rows,  " pathways)\n",
    "  GO source:             ", step03_go_src,  " (", step03_go_rows,  " pathways)\n",
    "  agreement_mode:        ", agree_mode_cfg, "\n",
    "\nFinal selection policy:\n",
    "  final_mode:            ", FINAL_MODE, "\n",
    "  final_top_n_per_dir:   ", FINAL_TOP_N, "\n",
    "\nStep 06 thresholds (from cfg):\n",
    "  CONS_JACCARD_MIN:      ", cfg$thresholds$CONS_JACCARD_MIN, "\n",
    "  CONS_K_MIN:            ", cfg$thresholds$CONS_K_MIN, "\n",
    "  solo_q_max:            ", cfg$pipelineB$solo_q_max, "\n",
    "  solo_abs_nes_min:      ", cfg$pipelineB$solo_abs_nes_min, "\n",
    "  optionA_enable:        ", isTRUE(cfg$pipelineB$optionA_enable), "\n",
    "  optionA_quality_floor: ", if (!is.null(cfg$pipelineB$optionA_quality_floor_q)) cfg$pipelineB$optionA_quality_floor_q else 0.20, "\n",
    "\nCandidate set counts (step 06 full tables):\n",
    "  UP   consensus rows:         ", n_cons_up,          "\n",
    "  UP   solo rows:              ", n_solo_up,          "\n",
    "  UP   fallback rows:          ", n_fallback_up,      "\n",
    "  UP   last_resort rows:       ", n_last_resort_up,   "\n",
    "  DOWN consensus rows:         ", n_cons_down,        "\n",
    "  DOWN solo rows:              ", n_solo_down,        "\n",
    "  DOWN fallback rows:          ", n_fallback_down,    "\n",
    "  DOWN last_resort rows:       ", n_last_resort_down, "\n",
    "\nFinal selection source used:\n",
    "  UP:                    ", source_up,   "\n",
    "  DOWN:                  ", source_down, "\n",
    "\nFINAL.csv counts (top ", FINAL_TOP_N, " per direction):\n",
    "  n_up:                  ", n_up,    "\n",
    "  n_down:                ", n_down,  "\n",
    "  n_total:               ", n_total, "\n",
    empty_final_note
  )

  writeLines(summary_text, file.path(FINAL_DIR, "final_summary.txt"))
  message("FINAL.csv written to:      ", file.path(FINAL_DIR, "FINAL.csv"))
  message("final_summary.txt written: ", file.path(FINAL_DIR, "final_summary.txt"))
  message("Step03 — WP: ", step03_wp_rows, " [", step03_wp_src, "]",
          "  GO: ", step03_go_rows, " [", step03_go_src, "]")
  message("Candidates — UP consensus: ", n_cons_up, "  solo: ", n_solo_up,
          "  fallback: ", n_fallback_up, "  last_resort: ", n_last_resort_up,
          " | DOWN consensus: ", n_cons_down, "  solo: ", n_solo_down,
          "  fallback: ", n_fallback_down, "  last_resort: ", n_last_resort_down)
  message("FINAL shortlist — UP: ", n_up, "  DOWN: ", n_down, "  total: ", n_total,
          "  (top ", FINAL_TOP_N, " per direction, mode: ", FINAL_MODE, ")")
  if (n_total == 0) message(trimws(empty_final_note))

} else {
  message("FINAL export skipped: step06 output files not found (was step 6 included in run_steps?)")
}

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------

message("")
message("=== Pipeline B run complete ===")
message("Run dir:   ", ctx$run_dir)
message("Dataset:   ", cfg$dataset$name)
message("Steps run: ",
        if (length(ran_steps) > 0) paste(sort(ran_steps), collapse = ", ") else "(none)")
message("Metadata:  ", file.path(ctx$run_dir, "run_meta.yml"))
