# R/ablation/run_ablation_B.R
# Pipeline B ablation runner.
#
# Runs one ablation variant of Pipeline B on a fixed dataset, writing outputs
# to results/ablation/<variant>/<dataset>/ (named, not timestamped).
#
# Usage (from project root):
#   Rscript R/ablation/run_ablation_B.R --variant B_full
#   Rscript R/ablation/run_ablation_B.R --variant B_no_fgsea --overwrite
#
# How it works:
#   1. Loads config/default.yml (dataset, thresholds, all pipelineB settings)
#   2. Merges config/ablation/<variant>.yml (only pipelineB.ablation block)
#   3. Runs steps 01-06 with bypass logic controlled by ablation flags
#   4. Runs FINAL export (top-N cap honoring ablation$top_n_per_direction)
#   5. Writes stage_counts.csv, run_note.txt, config_used.yml, run_meta.yml
#
# Bypass rules:
#   use_fgsea = false
#     -> step02 skipped; step03 uses CAMERA-only passthrough (fg_FDR=NA)
#     -> step04 ALWAYS uses passthrough when use_fgsea=false (fg_FDR=NA
#        breaks step04's !is.na(fg_FDR) filter regardless of use_semantic_collapse)
#
#   use_semantic_collapse = false
#     -> step04 passthrough: copies step03 GO -> step04/reps.csv unchanged
#
#   use_overlap_clustering = false
#     -> step05 passthrough: merges step03 WP + step04 GO, splits by
#        direction, writes representatives_*.csv without clustering.
#        set_size sourced from step01 CAMERA NGenes.
#
#   top_n_per_direction = null
#     -> FINAL export: slice_head() skipped; all keep_final==TRUE rows included.
#
#   use_consensus = false  [NOT YET IMPLEMENTED - reserved for future chunk]
#     -> currently ignored; step06 always applies consensus filter.

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ---------------------------------------------------------------------------
# Parse CLI arguments
# ---------------------------------------------------------------------------

args         <- commandArgs(trailingOnly = TRUE)
variant_name <- NULL
overwrite    <- FALSE
ablation_root <- "results/ablation"   # default; override with --root

i <- 1L
while (i <= length(args)) {
  if (args[i] == "--variant" && i < length(args)) {
    variant_name <- args[i + 1L]; i <- i + 2L
  } else if (args[i] == "--root" && i < length(args)) {
    ablation_root <- args[i + 1L]; i <- i + 2L
  } else if (args[i] == "--overwrite") {
    overwrite <- TRUE; i <- i + 1L
  } else {
    i <- i + 1L
  }
}

if (is.null(variant_name))
  stop("Usage: Rscript R/ablation/run_ablation_B.R --variant <name> [--root <dir>] [--overwrite]\n",
       "Available variants: B_full, B_no_fgsea, B_no_semantic, B_no_overlap, B_no_topcap")

variant_yml <- file.path("config/ablation", paste0(variant_name, ".yml"))
if (!file.exists(variant_yml))
  stop("Variant config not found: ", variant_yml,
       "\nExpected: config/ablation/", variant_name, ".yml")

# ---------------------------------------------------------------------------
# Source utilities and step functions
# ---------------------------------------------------------------------------

source("R/utils/config_loader.R")
source("R/utils/validate_inputs.R")
source("R/utils/gene_id_convert.R")
source("R/utils/go_term_names.R")
source("R/pipelineB/steps/step01_camera_wp_go_bp.R")
source("R/pipelineB/steps/step02_fgsea_wp_go_bp.R")
source("R/pipelineB/steps/step03_prepare_inputs_agreement_first.R")
source("R/pipelineB/steps/step04_go_collapse_semantic_agreement.R")
source("R/pipelineB/steps/step05_overlap_cluster_B.R")
source("R/pipelineB/steps/step06_consensus_and_themes_A1A2.R")
source("R/ablation/utils/ablation_config.R")
source("R/ablation/utils/ablation_passthroughs.R")
source("R/ablation/utils/ablation_stage_counts.R")

# ---------------------------------------------------------------------------
# Load and merge config
# ---------------------------------------------------------------------------

base_cfg <- load_config("config/default.yml")
cfg      <- merge_ablation_config(base_cfg, variant_yml)
abl      <- cfg$pipelineB$ablation   # ablation flags shorthand

dataset_name <- cfg$dataset$name

# ---------------------------------------------------------------------------
# Output directory setup (named, not timestamped)
# ---------------------------------------------------------------------------

run_dir <- file.path(ablation_root, abl$variant_name, dataset_name)

if (dir.exists(run_dir) && !overwrite) {
  existing <- list.files(run_dir, recursive = TRUE)
  if (length(existing) > 0)
    stop("Output directory already exists with content: ", run_dir, "\n",
         "Use --overwrite to replace. This will delete all previous outputs for this variant/dataset.")
}

dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Build per-step output paths (mirrors build_paths_b but uses ablation run_dir)
# ---------------------------------------------------------------------------

build_ablation_paths <- function(base_dir) {
  make_dir <- function(name) {
    p <- file.path(base_dir, name)
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    p
  }
  list(
    step01 = make_dir("01_camera"),
    step02 = make_dir("02_fgsea"),
    step03 = make_dir("03_prepare_inputs"),
    step04 = make_dir("04_go_collapse"),
    step05 = make_dir("05_overlap"),
    step06 = make_dir("06_consensus")
  )
}

paths <- build_ablation_paths(run_dir)

# ---------------------------------------------------------------------------
# Run context (seed + metadata; no timestamped folder needed)
# ---------------------------------------------------------------------------

timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
seed      <- if (!is.null(cfg$run$seed)) cfg$run$seed else 123L
set.seed(seed)

ctx <- list(
  run_id    = paste0("ablation_", abl$variant_name, "_", dataset_name, "_", timestamp),
  run_dir   = run_dir,
  timestamp = timestamp,
  seed      = seed
)

# Save config snapshot (merged base + ablation override)
yaml::write_yaml(cfg, file.path(run_dir, "config_used.yml"))

# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------

message(rep("=", 60))
message("Pipeline B Ablation Runner")
message(rep("=", 60))
message("Variant:   ", abl$variant_name)
message("Dataset:   ", dataset_name)
message("Run dir:   ", run_dir)
message("")
message("Component flags:")
message("  use_fgsea:              ", isTRUE(abl$use_fgsea))
message("  use_semantic_collapse:  ", isTRUE(abl$use_semantic_collapse))
message("  use_overlap_clustering: ", isTRUE(abl$use_overlap_clustering))
message("  use_consensus:          ", isTRUE(abl$use_consensus),
        if (!isTRUE(abl$use_consensus)) "  [NOTE: not yet implemented — step06 always runs]" else "")
message("  top_n_per_direction:    ",
        if (is.null(abl$top_n_per_direction)) "null (no cap)" else abl$top_n_per_direction)
message(rep("-", 60))

validate_config_and_inputs(cfg, "pipelineB")

# ---------------------------------------------------------------------------
# Step 01: CAMERA enrichment (always runs)
# ---------------------------------------------------------------------------

message("\n[Step 01] CAMERA enrichment...")
result01 <- step01_camera_wp_go_bp(cfg, paths, ctx)
message("--- Step01 complete ---")

# ---------------------------------------------------------------------------
# Steps 02 + 03: fgsea and agreement filter (conditional on use_fgsea)
# ---------------------------------------------------------------------------

if (isTRUE(abl$use_fgsea)) {

  message("\n[Step 02] fgsea sensitivity analysis...")
  result02 <- step02_fgsea_wp_go_bp(cfg, paths, ctx)
  message("--- Step02 complete ---")

  message("\n[Step 03] CAMERA \u2229 fgsea agreement filter...")
  result03 <- step03_prepare_inputs_agreement_first(cfg, paths, ctx)
  message("--- Step03 complete ---")

} else {

  message("\n[Step 02] BYPASSED (use_fgsea = false)")
  message("[Step 03] CAMERA-only passthrough...")
  result03 <- passthrough_step03_camera_only(cfg, paths)
  message("--- Step03 passthrough complete ---")

}

# ---------------------------------------------------------------------------
# Step 04: GO semantic collapse (conditional)
#
# IMPORTANT: step04 requires fg_FDR non-NA (filter: !is.na(fg_FDR)).
# When use_fgsea=false, fg_FDR=NA for all rows → step04 would return 0 reps.
# Therefore, the step04 passthrough is ALWAYS used when use_fgsea=false,
# regardless of the use_semantic_collapse flag.
# ---------------------------------------------------------------------------

use_step04_passthrough <- !isTRUE(abl$use_fgsea) || !isTRUE(abl$use_semantic_collapse)

if (!use_step04_passthrough) {

  message("\n[Step 04] GO semantic collapse...")
  result04 <- step04_go_collapse_semantic_agreement(cfg, paths, ctx)
  message("--- Step04 complete ---")

} else {

  pt_reason <- if (!isTRUE(abl$use_fgsea))
    "use_fgsea=false (fg_FDR=NA incompatible with step04)"
  else
    "use_semantic_collapse=false"
  message("\n[Step 04] PASSTHROUGH [", pt_reason, "]")
  result04 <- passthrough_step04_no_semantic(paths, reason = pt_reason)
  message("--- Step04 passthrough complete ---")

}

# ---------------------------------------------------------------------------
# Step 05: Overlap-based clustering (conditional)
# ---------------------------------------------------------------------------

if (isTRUE(abl$use_overlap_clustering)) {

  message("\n[Step 05] Overlap clustering...")
  result05 <- step05_overlap_cluster_B(cfg, paths, ctx)
  message("--- Step05 complete ---")

} else {

  message("\n[Step 05] PASSTHROUGH (use_overlap_clustering = false)")
  result05 <- passthrough_step05_no_overlap(cfg, paths)
  message("--- Step05 passthrough complete ---")

}

# ---------------------------------------------------------------------------
# Step 06: Consensus + Themes (always runs)
# NOTE: use_consensus=false is reserved for a future implementation.
#       Step06 always applies the consensus/solo filter for now.
# ---------------------------------------------------------------------------

message("\n[Step 06] Consensus + Themes...")
result06 <- step06_consensus_and_themes_A1A2(cfg, paths, ctx)
message("--- Step06 complete ---")

# ---------------------------------------------------------------------------
# FINAL export
# Replicates run_pipelineB.R FINAL block with these ablation-specific changes:
#   - FINAL_TOP_N from abl$top_n_per_direction (NULL = no cap)
#   - 'variant' and 'dataset' columns added to FINAL.csv
# ---------------------------------------------------------------------------

FINAL_DIR <- file.path(run_dir, "FINAL")
dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)

STEP06_UP   <- file.path(paths$step06, "final_up.csv")
STEP06_DOWN <- file.path(paths$step06, "final_down.csv")

# Final selection policy
FINAL_MODE  <- if (!is.null(cfg$pipelineB$final_mode))
  cfg$pipelineB$final_mode else "consensus_first"
FINAL_TOP_N <- abl$top_n_per_direction   # NULL = no cap (B_no_topcap)

apply_top_n <- function(tbl) {
  tbl <- tbl |>
    dplyr::arrange(agreement_q, dplyr::desc(absNES), camera_FDR, fgsea_padj)
  if (!is.null(FINAL_TOP_N))
    tbl <- dplyr::slice_head(tbl, n = as.integer(FINAL_TOP_N))
  tbl
}

if (file.exists(STEP06_UP) && file.exists(STEP06_DOWN)) {
  up   <- readr::read_csv(STEP06_UP,   show_col_types = FALSE)
  down <- readr::read_csv(STEP06_DOWN, show_col_types = FALSE)

  # Step03 provenance (for summary)
  detect_s3_src <- function(fpath) {
    if (!file.exists(fpath)) return(list(rows = 0L, src = "unknown"))
    tmp <- tryCatch(readr::read_csv(fpath, show_col_types = FALSE),
                    error = function(e) tibble::tibble())
    list(rows = nrow(tmp),
         src  = if ("agreement_source" %in% names(tmp) && nrow(tmp) > 0)
                  as.character(tmp$agreement_source[1]) else "unknown")
  }
  s3_wp  <- detect_s3_src(file.path(paths$step03, "wp_input_agreement.csv"))
  s3_go  <- detect_s3_src(file.path(paths$step03, "go_bp_input_agreement.csv"))

  # Candidate counts per tier
  count_reason <- function(tbl, reason_vals) {
    if (!"final_reason" %in% names(tbl)) return(0L)
    sum(tbl$final_reason %in% reason_vals, na.rm = TRUE)
  }
  n_cons_up     <- sum(up$in_consensus   == TRUE, na.rm = TRUE)
  n_solo_up     <- sum(up$solo_high_conf == TRUE, na.rm = TRUE)
  n_fb_up       <- count_reason(up, c("fallback_ranked", "fallback_ranked_low_conf"))
  n_lr_up       <- count_reason(up, "fallback_topn_last_resort")
  n_cons_down   <- sum(down$in_consensus   == TRUE, na.rm = TRUE)
  n_solo_down   <- sum(down$solo_high_conf == TRUE, na.rm = TRUE)
  n_fb_down     <- count_reason(down, c("fallback_ranked", "fallback_ranked_low_conf"))
  n_lr_down     <- count_reason(down, "fallback_topn_last_resort")

  # Selection policy: consensus_first
  pick_candidates <- function(tbl, n_cons, n_solo) {
    if (FINAL_MODE == "consensus_first") {
      if (n_cons >= 1L) dplyr::filter(tbl, in_consensus   == TRUE)
      else if (n_solo >= 1L) dplyr::filter(tbl, solo_high_conf == TRUE)
      else dplyr::filter(tbl, keep_final == TRUE)
    } else {
      dplyr::filter(tbl, keep_final == TRUE)
    }
  }

  candidates_up   <- pick_candidates(up,   n_cons_up,   n_solo_up)
  candidates_down <- pick_candidates(down, n_cons_down, n_solo_down)

  source_label <- function(n_cons, n_solo, n_fb, n_lr) {
    if (n_cons >= 1L) "consensus"
    else if (n_solo >= 1L) "solo"
    else if (n_fb >= 1L) "fallback_ranked"
    else if (n_lr >= 1L) "fallback_topn_last_resort"
    else "none"
  }
  source_up   <- source_label(n_cons_up,   n_solo_up,   n_fb_up,   n_lr_up)
  source_down <- source_label(n_cons_down, n_solo_down, n_fb_down, n_lr_down)

  final_up   <- apply_top_n(candidates_up)
  final_down <- apply_top_n(candidates_down)
  final_tbl  <- dplyr::bind_rows(final_up, final_down)

  # ---- Fill missing descriptions -------------------------------------------
  go_na_ids <- final_tbl |>
    dplyr::filter(is.na(Description), collection == "GO", startsWith(ID, "GO:")) |>
    dplyr::pull(ID)
  if (length(go_na_ids) > 0) {
    suppressPackageStartupMessages({ library(AnnotationDbi); library(GO.db) })
    go_raw    <- suppressMessages(
      AnnotationDbi::select(GO.db, keys = go_na_ids, keytype = "GOID", columns = "TERM")
    )
    go_lookup <- stats::setNames(go_raw$TERM, go_raw$GOID)
    final_tbl <- final_tbl |>
      dplyr::mutate(Description = dplyr::if_else(
        is.na(Description) & collection == "GO" & startsWith(ID, "GO:"),
        go_lookup[ID], Description))
  }
  wp_na_ids <- final_tbl |>
    dplyr::filter(is.na(Description), collection == "WP", startsWith(ID, "WP")) |>
    dplyr::pull(ID)
  if (length(wp_na_ids) > 0) {
    wp_cam_file <- file.path(paths$step01, "camera_wikipathways.csv")
    if (file.exists(wp_cam_file)) {
      wp_names  <- readr::read_csv(wp_cam_file, show_col_types = FALSE) |>
        dplyr::filter(!is.na(Description)) |>
        dplyr::distinct(ID, .keep_all = TRUE)
      wp_lookup <- stats::setNames(wp_names$Description, wp_names$ID)
      final_tbl <- final_tbl |>
        dplyr::mutate(Description = dplyr::if_else(
          is.na(Description) & collection == "WP" & startsWith(ID, "WP"),
          wp_lookup[ID], Description))
    }
  }
  # --------------------------------------------------------------------------

  # Add ablation metadata columns
  final_tbl <- final_tbl |>
    dplyr::mutate(variant = abl$variant_name, dataset = dataset_name)

  readr::write_csv(final_tbl, file.path(FINAL_DIR, "FINAL.csv"))

  n_up    <- nrow(final_up)
  n_down  <- nrow(final_down)
  n_total <- nrow(final_tbl)

  # ---- final_summary.txt ---------------------------------------------------
  top_n_str <- if (is.null(FINAL_TOP_N)) "null (no cap)" else as.character(FINAL_TOP_N)
  agree_mode_cfg <- if (!is.null(cfg$pipelineB$agreement_mode))
    cfg$pipelineB$agreement_mode else "intersection"

  summary_text <- paste0(
    "Pipeline B Ablation\n",
    "Variant:                 ", abl$variant_name, "\n",
    "Dataset:                 ", dataset_name, "\n",
    "Timestamp:               ", timestamp, "\n",
    "\nComponent flags:\n",
    "  use_fgsea:             ", isTRUE(abl$use_fgsea), "\n",
    "  use_semantic_collapse: ", isTRUE(abl$use_semantic_collapse), "\n",
    "  use_overlap_clustering:", isTRUE(abl$use_overlap_clustering), "\n",
    "  use_consensus:         ", isTRUE(abl$use_consensus), "\n",
    "  top_n_per_direction:   ", top_n_str, "\n",
    "\nStep04 passthrough applied: ", use_step04_passthrough, "\n",
    if (use_step04_passthrough && !isTRUE(abl$use_fgsea))
      "  (reason: fg_FDR=NA when use_fgsea=false)\n" else "",
    "\nStep03 agreement:\n",
    "  WP source: ", s3_wp$src, " (", s3_wp$rows, " pathways)\n",
    "  GO source: ", s3_go$src, " (", s3_go$rows, " pathways)\n",
    "  agreement_mode (base cfg): ", agree_mode_cfg, "\n",
    "\nFinal selection:\n",
    "  final_mode:            ", FINAL_MODE, "\n",
    "  top_n_per_direction:   ", top_n_str, "\n",
    "\nStep06 candidate counts:\n",
    "  UP   consensus: ", n_cons_up,   "  solo: ", n_solo_up,   "  fallback: ", n_fb_up,   "  last_resort: ", n_lr_up,   "\n",
    "  DOWN consensus: ", n_cons_down, "  solo: ", n_solo_down, "  fallback: ", n_fb_down, "  last_resort: ", n_lr_down, "\n",
    "\nFinal selection source:\n",
    "  UP:   ", source_up,   "\n",
    "  DOWN: ", source_down, "\n",
    "\nFINAL.csv counts:\n",
    "  n_up:    ", n_up,    "\n",
    "  n_down:  ", n_down,  "\n",
    "  n_total: ", n_total, "\n"
  )

  writeLines(summary_text, file.path(FINAL_DIR, "final_summary.txt"))
  message("FINAL.csv:         ", file.path(FINAL_DIR, "FINAL.csv"))
  message("final_summary.txt: ", file.path(FINAL_DIR, "final_summary.txt"))
  message("FINAL shortlist — UP: ", n_up, "  DOWN: ", n_down,
          "  total: ", n_total,
          if (is.null(FINAL_TOP_N)) "  (no cap)" else paste0("  (cap: ", FINAL_TOP_N, " per direction)"))

} else {
  message("FINAL export skipped: step06 output files not found.")
  n_up <- n_down <- n_total <- NA_integer_
  source_up <- source_down <- NA_character_
}

# ---------------------------------------------------------------------------
# run_meta.yml
# ---------------------------------------------------------------------------

git_hash <- tryCatch(
  trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE)),
  error = function(e) NA_character_
)

run_meta <- list(
  pipeline    = "pipelineB_ablation",
  variant     = abl$variant_name,
  dataset     = dataset_name,
  run_id      = ctx$run_id,
  timestamp   = timestamp,
  git_commit  = if (length(git_hash) == 1 && nchar(git_hash) > 0) git_hash else NA_character_,
  run_dir     = run_dir,
  flags = list(
    use_fgsea              = isTRUE(abl$use_fgsea),
    use_semantic_collapse  = isTRUE(abl$use_semantic_collapse),
    use_overlap_clustering = isTRUE(abl$use_overlap_clustering),
    use_consensus          = isTRUE(abl$use_consensus),
    top_n_per_direction    = abl$top_n_per_direction,
    step04_passthrough_applied = use_step04_passthrough
  )
)
yaml::write_yaml(run_meta, file.path(run_dir, "run_meta.yml"))

# ---------------------------------------------------------------------------
# run_note.txt (human-readable per-run record)
# ---------------------------------------------------------------------------

top_n_str <- if (is.null(abl$top_n_per_direction)) "null (no cap)" else abl$top_n_per_direction

run_note_lines <- c(
  paste0("Ablation variant:        ", abl$variant_name),
  paste0("Dataset:                 ", dataset_name),
  paste0("Run date:                ", format(Sys.Date(), "%Y-%m-%d")),
  paste0("Run dir:                 ", run_dir),
  "",
  "Components ON/OFF:",
  paste0("  use_fgsea:              ", isTRUE(abl$use_fgsea)),
  paste0("  use_semantic_collapse:  ", isTRUE(abl$use_semantic_collapse)),
  paste0("  use_overlap_clustering: ", isTRUE(abl$use_overlap_clustering)),
  paste0("  use_consensus:          ", isTRUE(abl$use_consensus)),
  paste0("  top_n_per_direction:    ", top_n_str),
  ""
)

if (use_step04_passthrough && !isTRUE(abl$use_semantic_collapse)) {
  run_note_lines <- c(run_note_lines,
    "Note: step04 passthrough applied (use_semantic_collapse=false).",
    "  All step03 GO terms passed to step05 without semantic collapse.", "")
}
if (use_step04_passthrough && !isTRUE(abl$use_fgsea)) {
  run_note_lines <- c(run_note_lines,
    "Note: step04 passthrough ALSO applied because use_fgsea=false.",
    "  fg_FDR=NA is incompatible with step04's !is.na(fg_FDR) filter.",
    "  Semantic collapse was therefore also bypassed for this variant.", "")
}
if (!isTRUE(abl$use_fgsea)) {
  run_note_lines <- c(run_note_lines,
    "Note: With fgsea bypassed, fgsea_NES=NA for all pathways.",
    "  solo_high_conf requires absNES >= solo_abs_nes_min — this rule",
    "  cannot activate without NES. Final pathways depend on GO<->WP",
    "  consensus only (or fallback mechanisms if agreement_mode triggers them).", "")
}

# Stage funnel (from stage_counts, collected below)
stage_counts_df <- collect_stage_counts(cfg, paths, abl, run_dir)

run_note_lines <- c(run_note_lines, "Stage funnel:")
for (i in seq_len(nrow(stage_counts_df))) {
  r    <- stage_counts_df[i, ]
  byp  <- if (r$bypassed) " [BYPASSED]" else ""
  note <- if (nzchar(r$note)) paste0(" (", r$note, ")") else ""
  cnt  <- if (is.na(r$count)) "NA" else as.character(r$count)
  run_note_lines <- c(run_note_lines,
    sprintf("  %-28s %-5s: %s%s%s",
            r$stage, r$direction, cnt, byp, note))
}

if (!is.na(n_total)) {
  run_note_lines <- c(run_note_lines, "",
    paste0("Final selection source:"),
    paste0("  UP:   ", source_up),
    paste0("  DOWN: ", source_down))
}

writeLines(run_note_lines, file.path(run_dir, "run_note.txt"))

# ---------------------------------------------------------------------------
# Final banner
# ---------------------------------------------------------------------------

message("")
message(rep("=", 60))
message("Pipeline B Ablation — COMPLETE")
message("Variant:   ", abl$variant_name)
message("Dataset:   ", dataset_name)
message("Run dir:   ", run_dir)
message("Metadata:  ", file.path(run_dir, "run_meta.yml"))
message(rep("=", 60))
