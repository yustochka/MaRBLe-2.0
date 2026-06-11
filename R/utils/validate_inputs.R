# R/utils/validate_inputs.R
# Pre-flight checks run at the start of both pipeline runners.
# Validates: required input files exist + critical packages available +
#            (Pipeline A) statistics file gene/p-value columns are resolvable.
# Does NOT install anything.

validate_config_and_inputs <- function(cfg, pipeline = c("pipelineA", "pipelineB")) {
  pipeline <- match.arg(pipeline)
  errors   <- character(0)

  message("Validating dataset: ", cfg$dataset$name, " (", pipeline, ")")

  # ---- 1. Required dataset files (all from cfg$dataset) ----------------------

  req_files <- list(
    expression_file   = cfg$dataset$expression_file,
    pathway2gene_file = cfg$dataset$pathway2gene_file
  )

  if (pipeline == "pipelineB") {
    req_files[["metadata_file"]] <- cfg$dataset$metadata_file
  }

  if (pipeline == "pipelineA") {
    req_files[["statistics_file"]] <- cfg$dataset$statistics_file
  }

  for (nm in names(req_files)) {
    f <- req_files[[nm]]
    if (is.null(f) || !nzchar(f)) {
      errors <- c(errors, sprintf("[file] %s: (null or empty — check config)", nm))
    } else if (!file.exists(f)) {
      errors <- c(errors, sprintf("[file] %s: not found — %s", nm, f))
    }
  }

  # ---- 2. Statistics file column pre-check (Pipeline A only) -----------------
  # Peeks at the statistics file header and reports which gene/p-value columns
  # will be used by step01 and step02.  Errors if a configured column is absent.

  if (pipeline == "pipelineA") {
    stats_path <- cfg$dataset$statistics_file
    if (!is.null(stats_path) && nzchar(stats_path) && file.exists(stats_path)) {

      stats_hdr <- tryCatch(
        readr::read_csv(stats_path, show_col_types = FALSE, n_max = 0),
        error = function(e) NULL
      )

      if (!is.null(stats_hdr)) {
        nms  <- names(stats_hdr)
        norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))

        # --- Gene column -------------------------------------------------------
        gc <- cfg$dataset$gene_col
        if (!is.null(gc) && nzchar(gc)) {
          # Config override: must exist in statistics file OR expression file
          if (gc %in% nms) {
            message("  gene_col: '", gc, "' (configured — found in statistics_file)")
          } else {
            expr_path <- cfg$dataset$expression_file
            found_in_expr <- !is.null(expr_path) && file.exists(expr_path) &&
              tryCatch({
                eh <- readr::read_csv(expr_path, show_col_types = FALSE, n_max = 0)
                gc %in% names(eh)
              }, error = function(e) FALSE)
            if (found_in_expr) {
              message("  gene_col: '", gc, "' (configured — found in expression_file)")
            } else {
              errors <- c(errors, sprintf(
                "[columns] gene_col '%s' not found in statistics_file or expression_file.\n    stats cols: %s",
                gc, paste(nms, collapse = ", ")))
            }
          }
        } else {
          # Auto-detect preview (mirrors step02 logic)
          g_exact <- c("gene", "genes", "symbol", "geneid", "gene_id",
                       "external_gene_name", "hgnc_symbol")
          g_fuzzy <- c("genesymbol", "gene_symbol", "hgncsymbol")
          detected_gc <- nms[tolower(nms) %in% g_exact]
          if (length(detected_gc) == 0)
            detected_gc <- nms[norm(nms) %in% norm(g_fuzzy)]
          if (length(detected_gc) > 0)
            message("  gene_col: '", detected_gc[1], "' (auto-detected)")
          else
            message("  gene_col: NOT DETECTED — consider setting cfg$dataset$gene_col")
        }

        # --- P-value column ----------------------------------------------------
        pc <- cfg$dataset$p_col
        if (!is.null(pc) && nzchar(pc)) {
          if (pc %in% nms)
            message("  p_col:    '", pc, "' (configured — found in statistics_file)")
          else
            errors <- c(errors, sprintf(
              "[columns] p_col '%s' not found in statistics_file.\n    stats cols: %s",
              pc, paste(nms, collapse = ", ")))
        } else {
          p_exact <- c("p", "pvalue", "p.value", "p-value", "raw_p", "p_val",
                       "pval", "pr(>|t|)", "padj", "p.adjust", "fdr", "adj.p.value")
          detected_pc <- nms[tolower(nms) %in% p_exact]
          if (length(detected_pc) == 0)
            detected_pc <- nms[norm(nms) %in% norm(c("adjpvalue", "pvalueadj"))]
          if (length(detected_pc) > 0)
            message("  p_col:    '", detected_pc[1], "' (auto-detected)")
          else
            message("  p_col:    NOT DETECTED — consider setting cfg$dataset$p_col")
        }
      }
    }
  }

  # ---- 2b. Metadata column check (Pipeline B only) ---------------------------
  # Verifies that sample_id_col and group_col exist in metadata_file before
  # step01 runs, so the user gets a clear error with the available columns.

  if (pipeline == "pipelineB") {
    meta_path <- cfg$dataset$metadata_file
    if (!is.null(meta_path) && nzchar(meta_path) && file.exists(meta_path)) {
      meta_hdr <- tryCatch(
        readr::read_csv(meta_path, show_col_types = FALSE, n_max = 0),
        error = function(e) NULL
      )
      if (!is.null(meta_hdr)) {
        nms <- names(meta_hdr)
        sid <- if (!is.null(cfg$dataset$sample_id_col) && nzchar(cfg$dataset$sample_id_col))
          cfg$dataset$sample_id_col else "SampleID"
        grp <- if (!is.null(cfg$dataset$group_col) && nzchar(cfg$dataset$group_col))
          cfg$dataset$group_col else "Genotype"

        if (sid %in% nms) {
          message("  sample_id_col: '", sid, "' (found in metadata_file)")
        } else {
          errors <- c(errors, sprintf(
            "[columns] sample_id_col '%s' not found in metadata_file.\n    Available: %s\n    Set cfg$dataset$sample_id_col in config.",
            sid, paste(nms, collapse = ", ")))
        }

        if (grp %in% nms) {
          message("  group_col:     '", grp, "' (found in metadata_file)")
        } else {
          errors <- c(errors, sprintf(
            "[columns] group_col '%s' not found in metadata_file.\n    Available: %s\n    Set cfg$dataset$group_col in config.",
            grp, paste(nms, collapse = ", ")))
        }
      }
    }

    # gene_id_type check: report value; for non-SYMBOL warn if rownames non-numeric
    gid <- if (!is.null(cfg$dataset$gene_id_type) && nzchar(cfg$dataset$gene_id_type))
      toupper(cfg$dataset$gene_id_type) else "SYMBOL"

    if (gid == "SYMBOL") {
      message("  gene_id_type:  '", gid, "' (no conversion)")
    } else if (gid %in% c("ENTREZ", "ENTREZID", "ENSEMBL")) {
      message("  gene_id_type:  '", gid, "' (->SYMBOL conversion enabled via org.Hs.eg.db)")
      # Quick sanity check: peek at a few row IDs
      expr_path <- cfg$dataset$expression_file
      if (!is.null(expr_path) && nzchar(expr_path) && file.exists(expr_path)) {
        peek <- tryCatch(
          readr::read_csv(expr_path, show_col_types = FALSE, n_max = 5),
          error = function(e) NULL
        )
        if (!is.null(peek) && ncol(peek) > 0) {
          sample_ids <- as.character(peek[[1]])
          if (gid %in% c("ENTREZ", "ENTREZID") &&
              !all(grepl("^[0-9]+$", sample_ids))) {
            message("  Warning: gene_id_type=ENTREZ but some row IDs are not purely numeric.",
                    " (first few: ", paste(head(sample_ids, 3), collapse = ", "), ")")
          }
        }
      }
    } else {
      errors <- c(errors, sprintf(
        "[gene_id_type] '%s' is not supported. Use SYMBOL, ENTREZ, or ENSEMBL.", gid))
    }
  }

  # ---- 3. Critical packages (requireNamespace only — no installs) ------------

  pkgs <- c(
    "yaml", "readr", "dplyr", "stringr", "tibble",
    "purrr", "glue", "limma", "org.Hs.eg.db", "GOSemSim", "igraph"
  )

  if (pipeline == "pipelineB") {
    pkgs <- c(pkgs, "fgsea", "AnnotationDbi", "GO.db")
  }

  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    errors <- c(errors,
                sprintf("[packages] not installed: %s",
                        paste(missing_pkgs, collapse = ", ")))
  }

  # ---- 4. Report -------------------------------------------------------------

  if (length(errors) > 0) {
    message("Pre-flight validation FAILED (", pipeline, "):")
    for (e in errors) message("  - ", e)
    stop("Fix the above issues before running the pipeline.", call. = FALSE)
  }

  message("Pre-flight validation passed (", pipeline, ").")
  invisible(TRUE)
}
