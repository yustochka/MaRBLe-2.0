# shiny/R/data_loader.R
# ------------------------------------------------------------------
# Data-loading helpers for the MaRBLe Pathway/GO Shortlist Explorer.
#
# Path safety: all functions resolve paths relative to the project
# root, so they work whether the calling R session's working directory
# is the project root OR the shiny/ subdirectory (Shiny's default
# when launched via shiny::runApp('shiny')).
#
# No external package dependencies — base R only.
# ------------------------------------------------------------------

# ── Internal helpers ────────────────────────────────────────────────────────

# Locate the project root by searching for the presence of data/ + results/.
# Checks the current directory first, then one level up (for shiny/ context).
.marble_root <- function() {
  for (d in c(".", "..")) {
    if (dir.exists(file.path(d, "data")) &&
        dir.exists(file.path(d, "results"))) {
      return(normalizePath(d))
    }
  }
  normalizePath(".")  # fallback: return current dir
}

# Resolve a path that may be relative to the project root.
# If the path already exists as-is, normalise and return it.
# Otherwise try prepending .marble_root(); return constructed path
# even when it doesn't exist (callers check existence themselves).
.abs <- function(rel_path) {
  if (is.na(rel_path) || !nzchar(rel_path)) return(rel_path)
  if (file.exists(rel_path) || dir.exists(rel_path)) {
    return(normalizePath(rel_path))
  }
  candidate <- file.path(.marble_root(), rel_path)
  if (file.exists(candidate) || dir.exists(candidate)) {
    return(normalizePath(candidate))
  }
  candidate  # return even if missing — caller decides what to do
}

# Safe single-value extraction from a data frame row.
.get <- function(df, col) {
  if (nrow(df) > 0 && col %in% names(df)) df[[col]][1L] else NA
}


# ── 1. detect_datasets ──────────────────────────────────────────────────────

#' Scan data/ for dataset subfolders.
#' Returns a sorted character vector of folder names; character(0) if none.
detect_datasets <- function(data_dir = "data") {
  abs_dir <- .abs(data_dir)
  if (!dir.exists(abs_dir)) {
    warning("detect_datasets: data directory not found at: ", abs_dir)
    return(character(0))
  }
  dirs <- list.dirs(abs_dir, recursive = FALSE, full.names = FALSE)
  sort(dirs[nzchar(dirs)])
}


# ── 2. load_dataset_summary ─────────────────────────────────────────────────

#' Read results/analysis/dataset_summary.csv.
#' Returns the data frame, or an empty data frame (with attribute "missing"=TRUE)
#' if the file is absent.
load_dataset_summary <- function(
    summary_path = "results/analysis/dataset_summary.csv") {

  abs_path <- .abs(summary_path)
  if (!file.exists(abs_path)) {
    warning("load_dataset_summary: file not found: ", abs_path)
    empty <- data.frame(
      dataset         = character(0),
      pipeline        = character(0),
      run_dir         = character(0),
      n_universe      = integer(0),
      n_de            = integer(0),
      n_final         = integer(0),
      final_source    = character(0),
      consensus_pairs = integer(0)
    )
    attr(empty, "missing") <- TRUE
    return(empty)
  }

  tryCatch(
    read.csv(abs_path, stringsAsFactors = FALSE),
    error = function(e) {
      warning("load_dataset_summary: read error: ", conditionMessage(e))
      structure(data.frame(), missing = TRUE)
    }
  )
}


# ── 3. load_threshold_policy ────────────────────────────────────────────────

#' Read results/analysis/threshold_policy.csv.
#' Returns the data frame, or an empty data frame (with attribute "missing"=TRUE)
#' if the file is absent.
load_threshold_policy <- function(
    policy_path = "results/analysis/threshold_policy.csv") {

  abs_path <- .abs(policy_path)
  if (!file.exists(abs_path)) {
    warning("load_threshold_policy: file not found: ", abs_path)
    empty <- data.frame(
      dataset                    = character(0),
      classification             = character(0),
      recommended_TAU_STABILITY  = numeric(0),
      recommended_CONS_JACCARD_MIN = numeric(0),
      confidence_flag            = character(0)
    )
    attr(empty, "missing") <- TRUE
    return(empty)
  }

  tryCatch(
    read.csv(abs_path, stringsAsFactors = FALSE),
    error = function(e) {
      warning("load_threshold_policy: read error: ", conditionMessage(e))
      structure(data.frame(), missing = TRUE)
    }
  )
}


# ── 4. find_run_dir ─────────────────────────────────────────────────────────

#' Locate the run output directory for a given dataset + pipeline.
#'
#' Resolution order:
#'   1. results/run_registry.csv  (Pipeline A only; has pipeline_a_dir column)
#'   2. results/analysis/dataset_summary.csv  (both pipelines; has run_dir column)
#'   3. Folder scan of results/pipelineA/ or results/pipelineB/
#'
#' Returns an absolute path string, or NA_character_ if nothing found.
find_run_dir <- function(dataset, pipeline, results_dir = "results") {
  root     <- .marble_root()
  abs_rdir <- .abs(results_dir)

  # Resolve a run_dir value that may be relative to the project root.
  try_dir <- function(val) {
    if (is.na(val) || !nzchar(as.character(val))) return(NULL)
    for (p in c(as.character(val), file.path(root, as.character(val)))) {
      if (dir.exists(p)) return(normalizePath(p))
    }
    NULL
  }

  # -- Strategy 1: run_registry.csv (Pipeline A only) ----------------------
  if (pipeline == "A") {
    reg_path <- file.path(abs_rdir, "run_registry.csv")
    if (file.exists(reg_path)) {
      reg <- tryCatch(read.csv(reg_path, stringsAsFactors = FALSE),
                      error = function(e) NULL)
      if (!is.null(reg) &&
          all(c("dataset", "pipeline_a_dir") %in% names(reg))) {
        rows <- reg[reg$dataset == dataset, , drop = FALSE]
        if (nrow(rows) > 0) {
          result <- try_dir(rows$pipeline_a_dir[1L])
          if (!is.null(result)) return(result)
        }
      }
    }
  }

  # -- Strategy 2: dataset_summary.csv (both pipelines) -------------------
  summ_path <- file.path(abs_rdir, "analysis", "dataset_summary.csv")
  if (file.exists(summ_path)) {
    summ <- tryCatch(read.csv(summ_path, stringsAsFactors = FALSE),
                     error = function(e) NULL)
    if (!is.null(summ) &&
        all(c("dataset", "pipeline", "run_dir") %in% names(summ))) {
      pipeline_name <- if (pipeline == "A") "pipelineA" else "pipelineB"
      rows <- summ[summ$dataset == dataset &
                   summ$pipeline == pipeline_name, , drop = FALSE]
      if (nrow(rows) > 0) {
        result <- try_dir(rows$run_dir[1L])
        if (!is.null(result)) return(result)
      }
    }
  }

  # -- Strategy 3: folder scan --------------------------------------------
  subdir   <- if (pipeline == "A") "pipelineA" else "pipelineB"
  scan_dir <- file.path(abs_rdir, subdir)
  if (dir.exists(scan_dir)) {
    candidates <- list.dirs(scan_dir, recursive = FALSE, full.names = TRUE)
    matches <- candidates[
      grepl(dataset, basename(candidates), fixed = TRUE) &
      !grepl("archive", basename(candidates), ignore.case = TRUE)
    ]
    if (length(matches) > 0) {
      # Pick newest by lexicographic sort (folder names embed timestamps)
      return(normalizePath(sort(matches, decreasing = TRUE)[1L]))
    }
  }

  NA_character_
}


# ── 5. find_final_file ──────────────────────────────────────────────────────

#' Resolve the path to FINAL.csv for the requested dataset/pipeline/preset.
#'
#' Returns a named list:
#'   $found            TRUE / FALSE
#'   $path             absolute path or NA_character_
#'   $message          human-readable status string
#'   $effective_preset the preset actually used (may differ from requested)
find_final_file <- function(dataset, pipeline, preset = "Default") {

  ok <- function(path, msg, eff = preset) {
    list(found = TRUE,  path = path,          message = msg, effective_preset = eff)
  }
  fail <- function(msg, eff = preset) {
    list(found = FALSE, path = NA_character_, message = msg, effective_preset = eff)
  }

  # ---- Recommended -------------------------------------------------------
  if (preset == "Recommended") {
    if (pipeline == "B") {
      # No policy run available for Pipeline B — fall back to Default silently
      res <- find_final_file(dataset, pipeline, "Default")
      res$effective_preset <- "Default (Recommended not available for Pipeline B)"
      res$message <- paste0(
        "Recommended preset is only available for Pipeline A in v1. ",
        "Showing Default results for Pipeline B instead."
      )
      return(res)
    }
    # Pipeline A: look in policy_runs/
    policy_csv <- .abs(file.path(
      "results", "analysis", "policy_runs", dataset, "FINAL.csv"
    ))
    if (file.exists(policy_csv)) {
      return(ok(policy_csv,
                paste0("Recommended (policy-optimised) results loaded for ", dataset, "."),
                "Recommended"))
    }
    return(fail(paste0(
      "Recommended policy run not found for ", dataset, ". ",
      "Run results/analysis/04_run_with_policy.R first, ",
      "or switch to the Default preset."
    )))
  }

  # ---- Default -----------------------------------------------------------
  run_dir <- find_run_dir(dataset, pipeline)
  if (is.na(run_dir)) {
    return(fail(paste0(
      "No results directory found for ", dataset, " (Pipeline ", pipeline, "). ",
      "Run Rscript R/run_pipeline", pipeline, ".R from the project root to generate results, ",
      "then reload this dashboard."
    )))
  }

  final_csv <- file.path(run_dir, "FINAL", "FINAL.csv")
  if (!file.exists(final_csv)) {
    return(fail(paste0(
      "Run directory found (", basename(run_dir), ") ",
      "but FINAL/FINAL.csv is missing. ",
      "The pipeline run may be incomplete."
    )))
  }

  ok(normalizePath(final_csv),
     paste0("Results loaded from: ", basename(run_dir), "."),
     "Default")
}


# ── 6. load_shortlist ───────────────────────────────────────────────────────

#' Load the final shortlist and add standardised column aliases.
#'
#' Standardised columns added (originals are kept unchanged):
#'   term_id, term_name, collection, q_value, gene_count,
#'   stability, direction, consensus, rank, source_pipeline
#'
#' Returns a named list:
#'   $found            TRUE / FALSE
#'   $data             data frame (empty if not found)
#'   $message          human-readable status string
#'   $final_file       path used (or NA)
#'   $effective_preset preset actually used
load_shortlist <- function(dataset, pipeline, preset = "Default") {

  ff <- find_final_file(dataset, pipeline, preset)

  if (!ff$found) {
    return(list(
      found            = FALSE,
      data             = data.frame(),
      message          = ff$message,
      final_file       = NA_character_,
      effective_preset = ff$effective_preset
    ))
  }

  df <- tryCatch(
    read.csv(ff$path, stringsAsFactors = FALSE),
    error = function(e) {
      warning("load_shortlist: could not read ", ff$path, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(df) || nrow(df) == 0) {
    return(list(
      found            = FALSE,
      data             = data.frame(),
      message          = if (is.null(df))
                           paste0("Could not read file: ", ff$path)
                         else
                           paste0("FINAL.csv is empty for ", dataset,
                                  " Pipeline ", pipeline, "."),
      final_file       = ff$path,
      effective_preset = ff$effective_preset
    ))
  }

  df <- .standardise_shortlist(df, pipeline)

  list(
    found            = TRUE,
    data             = df,
    message          = ff$message,
    final_file       = ff$path,
    effective_preset = ff$effective_preset
  )
}

# Internal: add standardised column aliases without modifying original columns.
.standardise_shortlist <- function(df, pipeline) {
  has <- function(col) col %in% names(df)

  # Identifiers
  df$term_id   <- if (has("ID"))          df$ID          else NA_character_
  df$term_name <- if (has("Description")) df$Description else NA_character_

  # Collection (already present in both pipelines; add if missing)
  if (!has("collection")) df$collection <- NA_character_

  # q_value: significance measure
  if (pipeline == "A") {
    df$q_value <- if (has("p.adjust"))   as.numeric(df$p.adjust)   else NA_real_
  } else {
    # Pipeline B: prefer agreement_q; fall back to camera_FDR then fgsea_padj
    df$q_value <- if (has("agreement_q") && !all(is.na(df$agreement_q))) {
      as.numeric(df$agreement_q)
    } else if (has("camera_FDR")) {
      as.numeric(df$camera_FDR)
    } else if (has("fgsea_padj")) {
      as.numeric(df$fgsea_padj)
    } else {
      NA_real_
    }
  }

  # Gene count
  df$gene_count <- if (pipeline == "A" && has("k")) {
    as.integer(df$k)
  } else if (pipeline == "B" && has("set_size")) {
    as.integer(df$set_size)
  } else {
    NA_integer_
  }

  # Stability: Pipeline A has it; Pipeline B does not
  df$stability <- if (pipeline == "A" && has("stability")) {
    as.numeric(df$stability)
  } else {
    NA_real_
  }

  # Direction: Pipeline B only
  df$direction <- if (pipeline == "B" && has("direction")) {
    df$direction
  } else {
    NA_character_
  }

  # Consensus
  df$consensus <- if (pipeline == "A" && has("consensus")) {
    as.logical(df$consensus)
  } else if (pipeline == "B" && has("in_consensus")) {
    as.logical(df$in_consensus)
  } else {
    NA
  }

  # Source pipeline label
  df$source_pipeline <- pipeline

  # Rank by q_value ascending (NA values go last)
  ord     <- order(df$q_value, na.last = TRUE)
  df      <- df[ord, , drop = FALSE]
  df$rank <- seq_len(nrow(df))
  rownames(df) <- NULL

  df
}


# ── 7. load_summary_cards ───────────────────────────────────────────────────

#' Assemble the five summary-card values for the dashboard.
#'
#' Draws from dataset_summary.csv and threshold_policy.csv.
#' Returns NA for any value that cannot be found — never crashes.
#'
#' Returns a named list:
#'   $dataset, $pipeline, $final_terms, $universe_size,
#'   $policy_class, $consensus_pairs, $preset
load_summary_cards <- function(dataset, pipeline, preset = "Default") {

  pipeline_name <- if (pipeline == "A") "pipelineA" else "pipelineB"

  # dataset_summary.csv row for this dataset + pipeline
  summ  <- load_dataset_summary()
  row_s <- if (nrow(summ) > 0 &&
               all(c("dataset", "pipeline") %in% names(summ))) {
    summ[summ$dataset == dataset & summ$pipeline == pipeline_name,
         , drop = FALSE]
  } else {
    data.frame()
  }

  # threshold_policy.csv row for this dataset
  pol   <- load_threshold_policy()
  row_p <- if (nrow(pol) > 0 && "dataset" %in% names(pol)) {
    pol[pol$dataset == dataset, , drop = FALSE]
  } else {
    data.frame()
  }

  # For the Recommended preset, final_terms comes from the policy-run FINAL.csv
  final_terms <- if (preset == "Recommended" && pipeline == "A") {
    ff <- find_final_file(dataset, pipeline, "Recommended")
    if (ff$found) {
      df_tmp <- tryCatch(read.csv(ff$path, stringsAsFactors = FALSE),
                         error = function(e) NULL)
      if (!is.null(df_tmp)) nrow(df_tmp) else .get(row_s, "n_final")
    } else {
      .get(row_s, "n_final")
    }
  } else {
    .get(row_s, "n_final")
  }

  list(
    dataset         = dataset,
    pipeline        = pipeline,
    final_terms     = final_terms,
    universe_size   = .get(row_s, "n_universe"),
    policy_class    = .get(row_p, "classification"),
    consensus_pairs = .get(row_s, "consensus_pairs"),
    preset          = preset
  )
}


# ── 8. Gene map helpers ─────────────────────────────────────────────────────

#' Find the Pathway2Gene CSV for a dataset.
#' Searches in dataset-specific locations; falls back to the shared
#' data/dataset_0/processed/Pathway2Gene.csv used by all datasets.
#' Returns NA_character_ if nothing is found.
find_pathway2gene_file <- function(dataset, data_dir = "data") {
  abs_dir <- .abs(data_dir)
  candidates <- c(
    file.path(abs_dir, dataset,    "processed", "Pathway2Gene.csv"),
    file.path(abs_dir, dataset,    "processed", "pathway2gene.csv"),
    file.path(abs_dir, dataset,    "Pathway2Gene.csv"),
    file.path(abs_dir, dataset,    "pathway2gene.csv"),
    # Shared file used by all datasets (from config)
    file.path(abs_dir, "dataset_0", "processed", "Pathway2Gene.csv")
  )
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p))
  }
  NA_character_
}

#' Load WikiPathways gene sets from Pathway2Gene.csv.
#' Returns a named list: names = pathway IDs (WP*), values = sorted character vectors of SYMBOL.
#' Returns an empty named list (with attribute "message") if the file or columns are missing.
#' Base R only.
load_pathway_gene_map <- function(dataset, data_dir = "data") {
  empty <- function(msg) { x <- list(); attr(x, "message") <- msg; x }

  p2g_path <- find_pathway2gene_file(dataset, data_dir)
  if (is.na(p2g_path))
    return(empty(paste0("Pathway2Gene.csv not found for dataset: ", dataset)))

  df <- tryCatch(
    read.csv(p2g_path, stringsAsFactors = FALSE),
    error = function(e) { warning("load_pathway_gene_map: ", conditionMessage(e)); NULL }
  )
  if (is.null(df) || nrow(df) == 0)
    return(empty("Pathway2Gene.csv is empty or could not be read."))

  # Identify pathway ID column
  id_col <- c("wpid", "pathway_id", "ID", "id", "GO_ID", "go_id", "GOID")
  id_col <- id_col[id_col %in% names(df)][1L]

  # Identify gene symbol column
  sym_col <- c("SYMBOL", "symbol", "gene_symbol", "gene", "Gene")
  sym_col <- sym_col[sym_col %in% names(df)][1L]

  if (is.na(id_col) || is.na(sym_col))
    return(empty(paste0(
      "Cannot find ID/symbol columns in Pathway2Gene.csv. Columns: ",
      paste(names(df), collapse = ", ")
    )))

  ids  <- df[[id_col]]
  syms <- df[[sym_col]]
  keep <- !is.na(ids) & nzchar(ids) & !is.na(syms) & nzchar(syms)

  groups <- split(syms[keep], ids[keep])
  result <- lapply(groups, function(x) sort(unique(x)))
  attr(result, "source") <- p2g_path
  result
}

#' Look up gene symbols for GO term IDs using org.Hs.eg.db (if available).
#' Returns a named list: names = GO IDs, values = sorted character vectors of SYMBOL.
#' Only processes IDs starting with "GO:".
#' Returns empty vectors for any ID not found; fails gracefully if packages unavailable.
load_go_gene_map <- function(term_ids) {
  go_ids <- unique(term_ids[grepl("^GO:", as.character(term_ids))])
  if (length(go_ids) == 0L) return(list())

  needs <- c("AnnotationDbi", "org.Hs.eg.db")
  if (any(!vapply(needs, requireNamespace, logical(1L), quietly = TRUE))) {
    warning("load_go_gene_map: AnnotationDbi/org.Hs.eg.db not available.")
    return(setNames(vector("list", length(go_ids)), go_ids))
  }

  tbl <- tryCatch(
    suppressMessages(suppressWarnings(
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys    = go_ids,
        columns = "SYMBOL",
        keytype = "GOALL"
      )
    )),
    error = function(e) { warning("load_go_gene_map: ", conditionMessage(e)); NULL }
  )

  # Build result: one entry per GO ID
  result <- setNames(vector("list", length(go_ids)), go_ids)
  if (!is.null(tbl) && nrow(tbl) > 0) {
    valid  <- !is.na(tbl$SYMBOL) & nzchar(tbl$SYMBOL) &
              !is.na(tbl$GOALL)  & nzchar(tbl$GOALL)
    groups <- split(tbl$SYMBOL[valid], tbl$GOALL[valid])
    cleaned <- lapply(groups, function(x) sort(unique(x)))
    for (nm in names(cleaned)) {
      if (nm %in% names(result)) result[[nm]] <- cleaned[[nm]]
    }
  }
  result
}


# ── 9. Pipeline funnel counts ───────────────────────────────────────────────

# Internal: safely read nrow of a CSV; returns NA_integer_ if missing/unreadable.
.csv_nrow <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_integer_)
  df <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) NA_integer_ else nrow(df)
}

# Internal: add two integers; NA-safe (returns available value if one is NA, NA if both are).
.add_safe <- function(a, b) {
  if (is.na(a) && is.na(b)) NA_integer_
  else if (is.na(a)) b
  else if (is.na(b)) a
  else a + b
}

# Internal: build funnel from dataset_summary.csv when stage files are unavailable.
.funnel_from_summary <- function(dataset, pipeline) {
  summ <- load_dataset_summary()
  pname <- if (pipeline == "A") "pipelineA" else "pipelineB"
  row   <- if (nrow(summ) > 0) summ[summ$dataset == dataset & summ$pipeline == pname, , drop = FALSE]
           else data.frame()

  if (nrow(row) == 0) {
    return(data.frame(stage  = "No data available",
                      count  = NA_integer_,
                      source = "not available",
                      stringsAsFactors = FALSE))
  }

  if (pipeline == "A") {
    nc <- as.integer(.get(row, "n_candidates"))
    n2 <- as.integer(.get(row, "n_tier2"))
    n1 <- as.integer(.get(row, "n_tier1"))
    nf <- as.integer(.get(row, "n_final"))
    data.frame(
      stage  = c("All candidates",
                 if (!is.na(n2) && n2 > 0) "Stable (Tier 2)" else NULL,
                 if (!is.na(n1) && n1 > 0) "Tier 1 (consensus)" else NULL,
                 "Final shortlist"),
      count  = c(nc,
                 if (!is.na(n2) && n2 > 0) n2 else NULL,
                 if (!is.na(n1) && n1 > 0) n1 else NULL,
                 nf),
      source = "dataset_summary (partial)",
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      stage  = "Final shortlist",
      count  = as.integer(.get(row, "n_final")),
      source = "dataset_summary (partial)",
      stringsAsFactors = FALSE
    )
  }
}

#' Load per-stage term counts for the pipeline funnel chart.
#'
#' Reads count data from per-step output CSVs in the run directory.
#' Falls back to dataset_summary.csv if stage files are unavailable.
#'
#' Returns a data.frame with columns: stage, count, source.
#' Rows are ordered first → last pipeline stage.
load_funnel_counts <- function(dataset, pipeline, preset = "Default") {
  rd <- find_run_dir(dataset, pipeline)

  # ── Pipeline A ────────────────────────────────────────────────────────────
  if (pipeline == "A") {
    if (is.na(rd)) return(.funnel_from_summary(dataset, "A"))

    n_ora <- .add_safe(
      .csv_nrow(file.path(rd, "03_ora", "go_bp", "ora_go_bp_filtered.csv")),
      .csv_nrow(file.path(rd, "03_ora", "wp",    "ora_wikipathways_filtered.csv"))
    )

    n_overlap <- .add_safe(
      .csv_nrow(file.path(rd, "05_overlap", "go_bp", "representatives.csv")),
      .csv_nrow(file.path(rd, "05_overlap", "wp",    "representatives.csv"))
    )

    n_stable <- .csv_nrow(file.path(rd, "07_tiers", "tier2_stable_only.csv"))
    if (is.na(n_stable)) {
      # fallback: count selected == TRUE in final_shortlist.csv
      fsl <- file.path(rd, "06_bootstrap_consensus", "final_shortlist.csv")
      if (file.exists(fsl)) {
        df <- tryCatch(read.csv(fsl, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(df) && "selected" %in% names(df))
          n_stable <- sum(as.logical(df$selected), na.rm = TRUE)
      }
    }

    n_final <- .csv_nrow(file.path(rd, "FINAL", "FINAL.csv"))

    # Build data frame (skip any all-NA stages)
    labels  <- c("ORA candidates", "After redundancy reduction",
                 "Stable candidates",  "Final shortlist")
    counts  <- c(n_ora, n_overlap, n_stable, n_final)
    sources <- ifelse(!is.na(counts), "stage file", "not available")
    sources[4] <- if (!is.na(n_final)) "final file" else "not available"

    df <- data.frame(stage  = labels,
                     count  = as.integer(counts),
                     source = sources,
                     stringsAsFactors = FALSE)
    df <- df[!is.na(df$count), , drop = FALSE]

    if (nrow(df) == 0) return(.funnel_from_summary(dataset, "A"))
    return(df)
  }

  # ── Pipeline B ────────────────────────────────────────────────────────────
  if (is.na(rd)) return(.funnel_from_summary(dataset, "B"))

  n_agree <- .add_safe(
    .csv_nrow(file.path(rd, "03_prepare_inputs", "go_bp_input_agreement.csv")),
    .csv_nrow(file.path(rd, "03_prepare_inputs", "wp_input_agreement.csv"))
  )

  n_overlap <- .add_safe(
    .csv_nrow(file.path(rd, "05_overlap", "representatives_up.csv")),
    .csv_nrow(file.path(rd, "05_overlap", "representatives_down.csv"))
  )

  n_final <- .csv_nrow(file.path(rd, "FINAL", "FINAL.csv"))

  labels  <- c("Agreement filter (CAMERA ∩ fgsea)",
               "After redundancy reduction",
               "Final shortlist (top-N cap)")
  counts  <- c(n_agree, n_overlap, n_final)
  sources <- ifelse(!is.na(counts), "stage file", "not available")
  sources[3] <- if (!is.na(n_final)) "final file" else "not available"

  df <- data.frame(stage  = labels,
                   count  = as.integer(counts),
                   source = sources,
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$count), , drop = FALSE]

  if (nrow(df) == 0) return(.funnel_from_summary(dataset, "B"))
  df
}


# ── 10. Network edges ───────────────────────────────────────────────────────

#' Load pairwise edges for the redundancy network chart.
#'
#' Pipeline A: uses cross_collection_pairs.csv (GO-WP consensus) and
#'             05_overlap/wp/edges.csv (WP-WP Jaccard), filtered to shortlist.
#' Pipeline B: uses 05_overlap/edges_up.csv + edges_down.csv, filtered to shortlist.
#' Fallback: cluster-based edges from the cluster column in shortlist_df.
#'
#' @param shortlist_df  standardised data frame from load_shortlist()$data, or NULL
#' @return list(found, nodes, edges, message, source)
load_network_edges <- function(dataset, pipeline,
                               preset      = "Default",
                               shortlist_df = NULL) {

  # ── Skeleton data frames ──────────────────────────────────────────────────
  empty_nodes <- data.frame(
    id = character(), label = character(), collection = character(),
    direction = character(), cluster = integer(), gene_count = integer(),
    stringsAsFactors = FALSE
  )
  empty_edges <- data.frame(
    from = character(), to = character(), weight = numeric(),
    edge_type = character(), shared_genes = integer(),
    stringsAsFactors = FALSE
  )
  fail <- function(msg) list(found = FALSE, nodes = empty_nodes,
                              edges = empty_edges, message = msg,
                              source = "none")

  # ── Build nodes from shortlist_df ─────────────────────────────────────────
  nodes <- if (!is.null(shortlist_df) && nrow(shortlist_df) > 0L) {
    data.frame(
      id         = as.character(shortlist_df$term_id),
      label      = as.character(shortlist_df$term_name),
      collection = as.character(shortlist_df$collection),
      direction  = if ("direction"  %in% names(shortlist_df))
                     as.character(shortlist_df$direction) else NA_character_,
      cluster    = if ("cluster"    %in% names(shortlist_df))
                     suppressWarnings(as.integer(shortlist_df$cluster)) else NA_integer_,
      gene_count = if ("gene_count" %in% names(shortlist_df))
                     as.integer(shortlist_df$gene_count) else NA_integer_,
      stringsAsFactors = FALSE
    )
  } else empty_nodes

  shortlist_ids <- nodes$id
  rd            <- find_run_dir(dataset, pipeline)

  if (nrow(nodes) == 0L) return(fail("No shortlist terms to display."))
  if (is.na(rd)) {
    return(list(found = TRUE, nodes = nodes, edges = empty_edges,
                message = "No run directory — showing terms without edges.",
                source = "shortlist only"))
  }

  edges       <- empty_edges
  source_used <- ""

  # ── Safe CSV nrow helper (returns NULL if file missing/unreadable) ─────────
  .read_csv_safe <- function(path) {
    if (!file.exists(path)) return(NULL)
    tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  }

  if (pipeline == "A") {
    # ── Priority 1: cross_collection_pairs.csv (GO ↔ WP consensus) ──────────
    ccp <- .read_csv_safe(
      file.path(rd, "06_bootstrap_consensus", "cross_collection_pairs.csv")
    )
    if (!is.null(ccp) && nrow(ccp) > 0L &&
        all(c("A", "B", "jaccard", "k") %in% names(ccp))) {
      ccp_f <- ccp[ccp$A %in% shortlist_ids & ccp$B %in% shortlist_ids, ]
      if (nrow(ccp_f) > 0L) {
        edges <- rbind(edges, data.frame(
          from = ccp_f$A, to = ccp_f$B, weight = ccp_f$jaccard,
          edge_type = "consensus", shared_genes = as.integer(ccp_f$k),
          stringsAsFactors = FALSE
        ))
        source_used <- "cross_collection_pairs.csv"
      }
    }

    # ── Priority 2: 05_overlap/wp/edges.csv (WP-WP Jaccard) ─────────────────
    wp_e <- .read_csv_safe(file.path(rd, "05_overlap", "wp", "edges.csv"))
    if (!is.null(wp_e) && nrow(wp_e) > 0L &&
        all(c("from", "to", "jaccard", "k") %in% names(wp_e))) {
      wp_f <- wp_e[wp_e$from %in% shortlist_ids & wp_e$to %in% shortlist_ids, ]
      if (nrow(wp_f) > 0L) {
        edges <- rbind(edges, data.frame(
          from = wp_f$from, to = wp_f$to, weight = wp_f$jaccard,
          edge_type = "overlap", shared_genes = as.integer(wp_f$k),
          stringsAsFactors = FALSE
        ))
        source_used <- if (nzchar(source_used))
          paste0(source_used, " + wp/edges.csv") else "05_overlap/wp/edges.csv"
      }
    }

  } else {
    # ── Pipeline B: edges_up.csv + edges_down.csv ────────────────────────────
    for (sfx in c("up", "down")) {
      epath <- file.path(rd, "05_overlap", paste0("edges_", sfx, ".csv"))
      e <- .read_csv_safe(epath)
      if (is.null(e) || nrow(e) == 0L) next
      if (!all(c("from", "to", "jaccard", "k") %in% names(e))) next
      ef <- e[e$from %in% shortlist_ids & e$to %in% shortlist_ids, ]
      if (nrow(ef) == 0L) next
      edges <- rbind(edges, data.frame(
        from = ef$from, to = ef$to, weight = ef$jaccard,
        edge_type = "overlap", shared_genes = as.integer(ef$k),
        stringsAsFactors = FALSE
      ))
      nm <- paste0("edges_", sfx, ".csv")
      source_used <- if (nzchar(source_used)) paste0(source_used, " + ", nm) else nm
    }
  }

  # ── Fallback: cluster-based edges ─────────────────────────────────────────
  if (nrow(edges) == 0L && any(!is.na(nodes$cluster))) {
    grps <- split(nodes$id, nodes$cluster)
    cl_edges <- do.call(rbind, lapply(grps, function(ids) {
      ids <- ids[!is.na(ids)]
      if (length(ids) < 2L) return(NULL)
      combns <- utils::combn(ids, 2L, simplify = FALSE)
      do.call(rbind, lapply(combns, function(p) {
        data.frame(from = p[1L], to = p[2L], weight = 0.5,
                   edge_type = "cluster", shared_genes = NA_integer_,
                   stringsAsFactors = FALSE)
      }))
    }))
    if (!is.null(cl_edges) && nrow(cl_edges) > 0L) {
      edges       <- cl_edges
      source_used <- "cluster assignments (fallback)"
    }
  }

  # ── Compose result ─────────────────────────────────────────────────────────
  msg <- if (nrow(edges) == 0L) {
    "No edges found between shortlisted terms."
  } else {
    paste0(nrow(edges), " edge(s) from: ", source_used, ".")
  }

  list(found = TRUE, nodes = nodes, edges = edges,
       message = msg, source = source_used)
}


# ── 11. load_dataset_meta ───────────────────────────────────────────────────

#' Load all display-level metadata for a single dataset.
#'
#' Pulls from threshold_policy.csv and dataset_summary.csv.
#' Accepts pre-loaded data frames to avoid repeated file reads.
#'
#' Returns a named list with keys:
#'   dataset, n_universe, de_fraction,
#'   classification, tau, cons_j, confidence,
#'   n_final_a, final_source_a, n_final_b, final_source_b
#' All values are NA when not found — never crashes.
load_dataset_meta <- function(dataset,
                               policy_df  = NULL,
                               summary_df = NULL) {
  pol  <- if (!is.null(policy_df))  policy_df  else load_threshold_policy()
  summ <- if (!is.null(summary_df)) summary_df else load_dataset_summary()

  pol_row <- if (nrow(pol) > 0 && "dataset" %in% names(pol))
               pol[pol$dataset == dataset, , drop = FALSE]
             else data.frame()

  summ_a  <- if (nrow(summ) > 0 &&
                 all(c("dataset", "pipeline") %in% names(summ)))
               summ[summ$dataset == dataset & summ$pipeline == "pipelineA",
                    , drop = FALSE]
             else data.frame()

  summ_b  <- if (nrow(summ) > 0 &&
                 all(c("dataset", "pipeline") %in% names(summ)))
               summ[summ$dataset == dataset & summ$pipeline == "pipelineB",
                    , drop = FALSE]
             else data.frame()

  list(
    dataset        = dataset,
    n_universe     = .get(summ_a, "n_universe"),
    de_fraction    = .get(summ_a, "de_fraction"),
    classification = .get(pol_row, "classification"),
    tau            = .get(pol_row, "recommended_TAU_STABILITY"),
    cons_j         = .get(pol_row, "recommended_CONS_JACCARD_MIN"),
    confidence     = .get(pol_row, "confidence_flag"),
    n_final_a      = .get(summ_a, "n_final"),
    final_source_a = .get(summ_a, "final_source"),
    n_final_b      = .get(summ_b, "n_final"),
    final_source_b = .get(summ_b, "final_source")
  )
}


# ── 12. load_comparison_summary ─────────────────────────────────────────────

#' Load one-row cross-pipeline comparison summary for a dataset.
#' Returns a named list or NULL if the dataset is not found.
load_comparison_summary <- function(ds) {
  f <- .abs("results/analysis/cross_pipeline/comparison_table.csv")
  if (!file.exists(f)) return(NULL)
  tbl <- tryCatch(
    read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(tbl) || nrow(tbl) == 0L || !"dataset" %in% names(tbl)) return(NULL)
  row <- tbl[tbl$dataset == ds, , drop = FALSE]
  if (nrow(row) == 0L) return(NULL)
  as.list(row[1L, ])
}


# ── 13. load_overlap_table ───────────────────────────────────────────────────

#' Load exact shared pathway IDs between Pipeline A and B for a dataset.
#' Returns a data frame (0 rows when none found).
load_overlap_table <- function(ds) {
  f <- .abs("results/analysis/cross_pipeline/pathway_overlap.csv")
  empty <- data.frame(ID = character(), Description = character(),
                      stability = numeric(), direction = character(),
                      absNES = numeric(), stringsAsFactors = FALSE)
  if (!file.exists(f)) return(empty)
  tbl <- tryCatch(
    read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) empty
  )
  if (nrow(tbl) == 0L || !"dataset" %in% names(tbl)) return(empty)
  tbl[tbl$dataset == ds, , drop = FALSE]
}
