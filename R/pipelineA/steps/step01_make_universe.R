# R/pipelineA/steps/step01_make_universe.R
# Pipeline A — Step 01: define measured∩mapped gene universe (SYMBOL)
#
# Reads:
#   cfg$dataset$expression_file
#   cfg$dataset$pathway2gene_file
#   cfg$dataset$statistics_file   (fallback for gene symbols when expression file
#                                   uses non-symbol IDs; only read if gene_col
#                                   override is set and not found in expr file)
#
# Column overrides (optional, in config under dataset:):
#   gene_col   exact column name for gene symbols (default: auto-detect)
#
# Writes into paths$step01:
#   universe_genes.txt
#   debug_measured_not_in_mapping.csv
#   debug_mapping_not_measured.csv
#   run_note.txt

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

step01_make_universe <- function(cfg, paths, ctx) {

  # ---- Inputs ----------------------------------------------------------------
  EXPR_FILE <- cfg$dataset$expression_file
  MAP_FILE  <- cfg$dataset$pathway2gene_file

  # ---- Outputs (all inside the current run folder) ---------------------------
  OUT_UNIVERSE_TXT <- file.path(paths$step01, "universe_genes.txt")
  OUT_DEBUG_MEAS   <- file.path(paths$step01, "debug_measured_not_in_mapping.csv")
  OUT_DEBUG_MAP    <- file.path(paths$step01, "debug_mapping_not_measured.csv")
  OUT_NOTE         <- file.path(paths$step01, "run_note.txt")

  # ---- Read inputs -----------------------------------------------------------
  expr <- suppressMessages(readr::read_csv(EXPR_FILE, show_col_types = FALSE))
  map  <- suppressMessages(readr::read_csv(MAP_FILE,  show_col_types = FALSE))

  # ---- Gene symbol column detection ------------------------------------------
  # Normalises names (strips punctuation/spaces) for fuzzy matching.
  # Catches: SYMBOL, Gene Symbol, gene_symbol, GeneSymbol, hgnc_symbol, etc.
  find_sym_col <- function(nms) {
    patterns <- c("symbol", "genesymbol", "hgncsymbol", "externalgene", "hugo")
    nml <- gsub("[^a-z0-9]", "", tolower(nms))
    hit <- nms[nml %in% patterns]
    if (length(hit) > 0) hit[1] else NA_character_
  }

  # ---- Measured genes --------------------------------------------------------
  # Priority: cfg$dataset$gene_col override  >  auto-detect in expression file.
  # If the override column is absent from the expression file, fall back to the
  # statistics file (Pipeline A always has cfg$dataset$statistics_file).
  gene_col_cfg  <- if (!is.null(cfg$dataset$gene_col) && nzchar(cfg$dataset$gene_col))
                     cfg$dataset$gene_col else NA_character_
  measured_file <- EXPR_FILE   # updated below if we fall back to stats file

  if (!is.na(gene_col_cfg)) {
    if (gene_col_cfg %in% names(expr)) {
      expr_sym_col <- gene_col_cfg
      measured <- expr[[expr_sym_col]] %>% as.character() %>% str_trim() %>%
                  discard(~ is.na(.x) || .x == "") %>% unique()
    } else {
      # Column not in expression file — look in statistics file
      stats_path <- cfg$dataset$statistics_file
      if (is.null(stats_path) || !file.exists(stats_path))
        stop("gene_col '", gene_col_cfg, "' not found in expression file and ",
             "cfg$dataset$statistics_file is missing or not found.\n",
             "  Expression file columns: ", paste(names(expr), collapse = ", "))
      stats_tmp <- suppressMessages(readr::read_csv(stats_path, show_col_types = FALSE))
      if (!gene_col_cfg %in% names(stats_tmp))
        stop("gene_col '", gene_col_cfg, "' not found in expression file or statistics file.\n",
             "  Expression file columns:  ", paste(names(expr), collapse = ", "), "\n",
             "  Statistics file columns:  ", paste(names(stats_tmp), collapse = ", "))
      expr_sym_col  <- gene_col_cfg
      measured_file <- stats_path
      measured <- stats_tmp[[gene_col_cfg]] %>% as.character() %>% str_trim() %>%
                  discard(~ is.na(.x) || .x == "") %>% unique()
    }
  } else {
    # Auto-detect a SYMBOL-like column in the expression file
    expr_sym_col <- find_sym_col(names(expr))
    if (is.na(expr_sym_col)) expr_sym_col <- names(expr)[1]
    measured <- expr[[expr_sym_col]] %>% as.character() %>% str_trim() %>%
                discard(~ is.na(.x) || .x == "") %>% unique()
    if (length(measured) > 0 &&
        all(grepl("^[0-9]+$", measured[!is.na(measured) & nchar(measured) > 0]))) {
      warning("Column '", expr_sym_col, "' in the expression file contains numeric gene IDs, ",
              "not gene symbols. Universe will likely be empty. ",
              "Set cfg$dataset$gene_col to a SYMBOL column (e.g. 'Gene Symbol').")
    }
  }

  # ---- Mapping genes (Pathway2Gene) ------------------------------------------
  map_sym_col <- find_sym_col(names(map))
  if (is.na(map_sym_col)) map_sym_col <- names(map)[1]

  mapped <- map[[map_sym_col]] %>%
    as.character() %>%
    str_trim() %>%
    discard(~ is.na(.x) || .x == "") %>%
    unique()

  # ---- Compute universe + debug sets -----------------------------------------
  universe            <- sort(intersect(measured, mapped))
  measured_not_mapped <- sort(setdiff(measured, mapped))
  mapped_not_measured <- sort(setdiff(mapped, measured))

  # Coverage
  cov_measured <- if (length(measured) > 0) length(universe) / length(measured) else NA_real_
  cov_mapped   <- if (length(mapped)   > 0) length(universe) / length(mapped)   else NA_real_

  # ---- Write outputs ---------------------------------------------------------
  readr::write_lines(universe, OUT_UNIVERSE_TXT)

  readr::write_csv(tibble(SYMBOL = measured_not_mapped), OUT_DEBUG_MEAS)
  readr::write_csv(tibble(SYMBOL = mapped_not_measured), OUT_DEBUG_MAP)

  measured_src_label <- if (measured_file == EXPR_FILE) "expression_file" else "statistics_file"

  note <- c(
    "Pipeline A — Universe (measured ∩ mapped)",
    paste0("dataset:      ", cfg$dataset$name),
    paste0("EXPR_FILE:    ", EXPR_FILE),
    paste0("MAP_FILE:     ", MAP_FILE),
    paste0("gene column:  '", expr_sym_col, "' (from ", measured_src_label, ")"),
    paste0("map  SYMBOL column used: ", map_sym_col),
    "",
    paste0("Measured genes: ", length(measured)),
    paste0("Mapped genes:   ", length(mapped)),
    paste0("Universe (N):   ", length(universe)),
    "",
    paste0("Measured not mapped: ", length(measured_not_mapped)),
    paste0("Mapped not measured: ", length(mapped_not_measured)),
    "",
    paste0("Coverage (Universe / measured): ", sprintf("%.1f%%", 100 * cov_measured)),
    paste0("Coverage (Universe / mapped):   ", sprintf("%.1f%%", 100 * cov_mapped)),
    "",
    paste0("Wrote: ", OUT_UNIVERSE_TXT),
    paste0("Wrote: ", OUT_DEBUG_MEAS),
    paste0("Wrote: ", OUT_DEBUG_MAP)
  )
  writeLines(note, OUT_NOTE)

  # ---- Console summary -------------------------------------------------------
  cat("\n=== Pipeline A Step 01: Universe ===\n")
  cat("Dataset:  ", cfg$dataset$name, "\n", sep = "")
  cat("Gene col: '", expr_sym_col, "' (from ", measured_src_label, ")\n", sep = "")
  cat("Measured: ", length(measured), "\n", sep = "")
  cat("Mapped:   ", length(mapped), "\n", sep = "")
  cat("Universe: ", length(universe), "\n", sep = "")
  cat("Wrote:    ", OUT_UNIVERSE_TXT, "\n\n", sep = "")

  if (length(universe) < 500) {
    warning("Universe looks very small (<500). Check that both inputs use the same ID type (SYMBOL).")
  }

  invisible(list(universe_file = OUT_UNIVERSE_TXT))
}
