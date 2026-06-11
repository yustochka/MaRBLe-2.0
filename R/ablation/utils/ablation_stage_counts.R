# R/ablation/utils/ablation_stage_counts.R
# Collect pathway counts at each pipeline stage from completed step output files
# and write stage_counts.csv to the run directory.
#
# Columns: variant, dataset, stage, direction, count, bypassed, note
#
# Stage IDs (in pipeline order):
#   camera_wp_sig      CAMERA WP at FDR <= 0.05
#   camera_go_sig      CAMERA GO:BP at FDR <= 0.05
#   fgsea_wp_sig       fgsea WP at padj <= 0.05   (NA/bypassed if use_fgsea=false)
#   fgsea_go_sig       fgsea GO:BP at padj <= 0.05
#   after_agreement_wp WP terms passing step03
#   after_agreement_go GO terms passing step03
#   after_semantic_go  GO reps after step04 (= agreement count if bypassed)
#   after_overlap_reps step05 reps, split by direction
#   final              FINAL.csv shortlist, split by direction

collect_stage_counts <- function(cfg, paths, ablation_cfg, run_dir) {

  variant <- ablation_cfg$variant_name
  dataset <- cfg$dataset$name

  rows <- list()

  add_row <- function(stage, direction, count, bypassed = FALSE, note = "") {
    rows[[length(rows) + 1]] <<- tibble::tibble(
      variant   = variant,
      dataset   = dataset,
      stage     = stage,
      direction = direction,
      count     = as.integer(count),
      bypassed  = as.logical(bypassed),
      note      = as.character(note)
    )
  }

  # Safe CSV reader — returns NULL on error
  safe_read <- function(fpath) {
    if (!file.exists(fpath)) return(NULL)
    tryCatch(readr::read_csv(fpath, show_col_types = FALSE),
             error = function(e) NULL)
  }

  # Count rows in a file matching a column condition
  count_rows <- function(df, col, op, val) {
    if (is.null(df) || !col %in% names(df)) return(NA_integer_)
    sum(op(df[[col]], val), na.rm = TRUE)
  }

  # ---- Step 01: CAMERA significant (FDR <= 0.05) ---------------------------
  cam_wp <- safe_read(file.path(paths$step01, "camera_wikipathways.csv"))
  cam_go <- safe_read(file.path(paths$step01, "camera_go_bp.csv"))

  cam_wp_fdr_col <- if (!is.null(cam_wp) && "FDR" %in% names(cam_wp)) "FDR" else "adj.P.Val"
  cam_go_fdr_col <- if (!is.null(cam_go) && "FDR" %in% names(cam_go)) "FDR" else "adj.P.Val"

  add_row("camera_wp_sig", "all",
          count_rows(cam_wp, cam_wp_fdr_col, `<=`, 0.05))
  add_row("camera_go_sig", "all",
          count_rows(cam_go, cam_go_fdr_col, `<=`, 0.05))

  # ---- Step 02: fgsea significant (padj <= 0.05) ---------------------------
  fgsea_bypassed <- !isTRUE(ablation_cfg$use_fgsea)

  if (fgsea_bypassed) {
    add_row("fgsea_wp_sig", "all", NA_integer_, bypassed = TRUE,
            note = "use_fgsea=false")
    add_row("fgsea_go_sig", "all", NA_integer_, bypassed = TRUE,
            note = "use_fgsea=false")
  } else {
    fg_wp <- safe_read(file.path(paths$step02, "fgsea_wikipathways.csv"))
    fg_go <- safe_read(file.path(paths$step02, "fgsea_go_bp.csv"))
    add_row("fgsea_wp_sig", "all", count_rows(fg_wp, "padj", `<=`, 0.05))
    add_row("fgsea_go_sig", "all", count_rows(fg_go, "padj", `<=`, 0.05))
  }

  # ---- Step 03: agreement filter -------------------------------------------
  wp_agree <- safe_read(file.path(paths$step03, "wp_input_agreement.csv"))
  go_agree <- safe_read(file.path(paths$step03, "go_bp_input_agreement.csv"))

  wp_agree_src <- if (!is.null(wp_agree) && "agreement_source" %in% names(wp_agree) && nrow(wp_agree) > 0)
    as.character(wp_agree$agreement_source[1]) else "unknown"
  go_agree_src <- if (!is.null(go_agree) && "agreement_source" %in% names(go_agree) && nrow(go_agree) > 0)
    as.character(go_agree$agreement_source[1]) else "unknown"

  add_row("after_agreement_wp", "all",
          if (is.null(wp_agree)) NA_integer_ else nrow(wp_agree),
          note = wp_agree_src)
  add_row("after_agreement_go", "all",
          if (is.null(go_agree)) NA_integer_ else nrow(go_agree),
          note = go_agree_src)

  # ---- Step 04: semantic collapse ------------------------------------------
  # Bypassed when use_semantic_collapse=false OR when use_fgsea=false
  # (fg_FDR=NA incompatible with step04's filter).
  semantic_bypassed <- !isTRUE(ablation_cfg$use_semantic_collapse) || fgsea_bypassed
  reps04 <- safe_read(file.path(paths$step04, "reps.csv"))

  add_row("after_semantic_go", "all",
          if (is.null(reps04)) NA_integer_ else nrow(reps04),
          bypassed = semantic_bypassed,
          note = dplyr::case_when(
            !isTRUE(ablation_cfg$use_semantic_collapse) ~ "use_semantic_collapse=false",
            fgsea_bypassed                              ~ "forced passthrough (fg_FDR=NA)",
            TRUE                                        ~ ""
          ))

  # ---- Step 05: overlap clustering representatives -------------------------
  overlap_bypassed <- !isTRUE(ablation_cfg$use_overlap_clustering)
  reps_up   <- safe_read(file.path(paths$step05, "representatives_up.csv"))
  reps_down <- safe_read(file.path(paths$step05, "representatives_down.csv"))

  overlap_note <- if (overlap_bypassed) "use_overlap_clustering=false" else ""
  add_row("after_overlap_reps", "Up",
          if (is.null(reps_up))   NA_integer_ else nrow(reps_up),
          bypassed = overlap_bypassed, note = overlap_note)
  add_row("after_overlap_reps", "Down",
          if (is.null(reps_down)) NA_integer_ else nrow(reps_down),
          bypassed = overlap_bypassed, note = overlap_note)

  # ---- FINAL shortlist -----------------------------------------------------
  final_file <- file.path(run_dir, "FINAL", "FINAL.csv")
  final_df   <- safe_read(final_file)

  add_row("final", "Up",
          if (is.null(final_df)) NA_integer_
          else sum(toupper(final_df$direction) == "UP", na.rm = TRUE))
  add_row("final", "Down",
          if (is.null(final_df)) NA_integer_
          else sum(toupper(final_df$direction) == "DOWN", na.rm = TRUE))

  # ---- Write ---------------------------------------------------------------
  out      <- dplyr::bind_rows(rows)
  out_path <- file.path(run_dir, "stage_counts.csv")
  readr::write_csv(out, out_path)
  message("stage_counts.csv written: ", out_path)
  invisible(out)
}
