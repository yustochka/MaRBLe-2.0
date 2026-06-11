# R/pipelineB/steps/step03_prepare_inputs_agreement_first.R
# Pipeline B — Step 03: Prepare high-confidence pathway inputs
# Keep pathways significant in BOTH CAMERA and fgsea AND with matching direction.
#
# Reads (from previous steps in the current run):
#   paths$step01/camera_wikipathways.csv
#   paths$step01/camera_go_bp.csv
#   paths$step02/fgsea_wikipathways.csv
#   paths$step02/fgsea_go_bp.csv
#
# Writes into paths$step03  (results/pipelineB/<run_id>/03_prepare_inputs/):
#   wp_input_agreement.csv    columns: ID, Description, direction, cam_FDR, fg_FDR, fgsea_NES
#   go_bp_input_agreement.csv columns: ID, Description, direction, cam_FDR, fg_FDR, fgsea_NES
#   run_note.txt
#
# Returns (invisibly):
#   wp_agree_file    path to wp_input_agreement.csv
#   go_agree_file    path to go_bp_input_agreement.csv
#
# agreement_mode (cfg$pipelineB$agreement_mode):
#   "intersection"        — strict CAMERA∩fgsea at FDR/padj <= 0.05 (default)
#   "relax_fdr_then_topn" — if strict yields 0, retry with agreement_relaxed_fdr;
#                           if still 0, take top agreement_topn per direction

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(glue)
})

step03_prepare_inputs_agreement_first <- function(cfg, paths, ctx) {

  # ---- Input paths (from previous steps) ------------------------------------
  CAM_WP_FILE <- file.path(paths$step01, "camera_wikipathways.csv")
  CAM_GO_FILE <- file.path(paths$step01, "camera_go_bp.csv")
  FG_WP_FILE  <- file.path(paths$step02, "fgsea_wikipathways.csv")
  FG_GO_FILE  <- file.path(paths$step02, "fgsea_go_bp.csv")

  stopifnot(
    file.exists(CAM_WP_FILE), file.exists(CAM_GO_FILE),
    file.exists(FG_WP_FILE),  file.exists(FG_GO_FILE)
  )

  OUT_DIR <- paths$step03

  # ---- Agreement mode parameters (Task 2) ------------------------------------
  AGREE_MODE  <- if (!is.null(cfg$pipelineB$agreement_mode) && nzchar(cfg$pipelineB$agreement_mode))
    cfg$pipelineB$agreement_mode else "intersection"
  AGREE_RELAX <- if (!is.null(cfg$pipelineB$agreement_relaxed_fdr))
    as.numeric(cfg$pipelineB$agreement_relaxed_fdr) else 0.20
  AGREE_TOPN  <- if (!is.null(cfg$pipelineB$agreement_topn))
    as.integer(cfg$pipelineB$agreement_topn) else 50L

  # ---- Helpers ---------------------------------------------------------------

  normalize_cam <- function(df) {
    has_fdr <- "FDR" %in% names(df)
    df |>
      dplyr::rename_with(~ "ID", dplyr::matches("^ID$|^Pathway$|^pathway$")) |>
      dplyr::mutate(
        cam_FDR = if (has_fdr) FDR else adj.P.Val,
        cam_dir = dplyr::case_when(
          Direction %in% c("Up", "UP")    ~  1,
          Direction %in% c("Down","DOWN") ~ -1,
          TRUE ~ NA_real_
        )
      ) |>
      dplyr::select(ID, cam_FDR, cam_dir, dplyr::everything())
  }

  normalize_fgsea <- function(df) {
    df |>
      dplyr::rename_with(~ "ID", dplyr::matches("^ID$|^pathway$")) |>
      dplyr::mutate(
        fg_FDR = padj,
        fg_dir = sign(NES)
      ) |>
      dplyr::select(ID, fg_FDR, fg_dir, dplyr::everything())
  }

  agreement_filter <- function(cam_df, fg_df, out_path, label) {

    cam_df <- normalize_cam(cam_df)
    fg_df  <- normalize_fgsea(fg_df)

    # Base: inner-join on ID, require same direction (no FDR cut yet)
    merged_base <- cam_df |>
      dplyr::inner_join(fg_df, by = "ID", suffix = c("_cam", "_fg")) |>
      dplyr::filter(!is.na(cam_dir), !is.na(fg_dir), cam_dir == fg_dir)

    # Strict filter (always applied first)
    merged    <- merged_base |> dplyr::filter(cam_FDR <= 0.05, fg_FDR <= 0.05)
    agree_src <- "strict_intersection"
    fallback_note <- ""

    # Fallback (only when strict yields 0 and mode != "intersection")
    if (nrow(merged) == 0 && AGREE_MODE != "intersection") {
      merged_relax <- merged_base |>
        dplyr::filter(cam_FDR <= AGREE_RELAX, fg_FDR <= AGREE_RELAX)

      if (nrow(merged_relax) > 0) {
        merged        <- merged_relax
        agree_src     <- "relaxed_fdr"
        fallback_note <- sprintf(" [fallback: relaxed FDR \u2264 %.2f]", AGREE_RELAX)
      } else if (AGREE_MODE == "relax_fdr_then_topn" && nrow(merged_base) > 0) {
        merged <- merged_base |>
          dplyr::mutate(.aq = pmax(cam_FDR, fg_FDR, na.rm = TRUE)) |>
          dplyr::group_by(cam_dir) |>
          dplyr::slice_min(.aq, n = AGREE_TOPN, with_ties = FALSE) |>
          dplyr::ungroup() |>
          dplyr::select(-.aq)
        agree_src     <- "topn"
        fallback_note <- sprintf(" [fallback: top-%d per direction]", AGREE_TOPN)
      }
    }

    # Description resolution (safe on 0-row tibbles)
    has_desc_cam <- "Description_cam" %in% names(merged)
    has_desc_fg  <- "Description_fg"  %in% names(merged)

    merged <- merged |>
      dplyr::mutate(
        Description = dplyr::coalesce(
          if (has_desc_cam) .data[["Description_cam"]] else NA_character_,
          if (has_desc_fg)  .data[["Description_fg"]]  else NA_character_
        )
      )

    out <- merged |>
      dplyr::transmute(
        ID,
        Description,
        direction = ifelse(cam_dir == 1, "Up", "Down"),
        cam_FDR,
        fg_FDR,
        fgsea_NES                  = if ("NES" %in% names(merged)) .data[["NES"]] else NA_real_,
        agreement_source           = agree_src,
        agreement_mode_used        = AGREE_MODE,
        agreement_relaxed_fdr_used = AGREE_RELAX,
        agreement_topn_used        = as.integer(AGREE_TOPN)
      ) |>
      dplyr::arrange(cam_FDR, fg_FDR)

    # Fill GO term descriptions before writing (no-op for WP IDs)
    if (!exists("add_go_term_name", mode = "function"))
      source("R/utils/go_term_names.R")
    out <- add_go_term_name(out, id_col = "ID", desc_col = "Description")

    readr::write_csv(out, out_path)

    cat(glue::glue(
      "[Agreement filter \u2014 {label}]\n",
      "CAMERA \u2229 fgsea (FDR \u2264 0.05, same direction){fallback_note}\n",
      "Agreement source used: {agree_src}\n",
      "Kept pathways: {nrow(out)}\n",
      "Output: {out_path}\n\n"
    ))

    out
  }

  # ---- Run -------------------------------------------------------------------
  cam_wp <- readr::read_csv(CAM_WP_FILE, show_col_types = FALSE)
  fg_wp  <- readr::read_csv(FG_WP_FILE,  show_col_types = FALSE)

  cam_go <- readr::read_csv(CAM_GO_FILE, show_col_types = FALSE)
  fg_go  <- readr::read_csv(FG_GO_FILE,  show_col_types = FALSE)

  wp_out_path <- file.path(OUT_DIR, "wp_input_agreement.csv")
  go_out_path <- file.path(OUT_DIR, "go_bp_input_agreement.csv")

  wp_agree <- agreement_filter(cam_wp, fg_wp, wp_out_path, "WikiPathways")
  go_agree <- agreement_filter(cam_go, fg_go, go_out_path, "GO:BP")

  # ---- Run note --------------------------------------------------------------
  note <- glue::glue("
[Pipeline B \u2014 Agreement inputs]

agreement_mode: {AGREE_MODE}
Strict criteria:
- CAMERA FDR \u2264 0.05
- fgsea padj \u2264 0.05
- Same enrichment direction

WikiPathways kept: {nrow(wp_agree)}
GO:BP kept:        {nrow(go_agree)}

Outputs:
- {wp_out_path}
- {go_out_path}
")

  writeLines(note, file.path(OUT_DIR, "run_note.txt"))
  cat(note, "\nDone.\n")

  invisible(list(
    wp_agree_file = wp_out_path,
    go_agree_file = go_out_path
  ))
}
