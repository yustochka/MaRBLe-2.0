# R/analysis/00_diagnostics.R
#
# Scans existing baseline results for all 4 datasets × 2 pipelines.
# Reads only; never re-runs any pipeline step.
#
# Output: results/analysis/dataset_summary.csv
#
# Run from project root:
#   Rscript R/analysis/00_diagnostics.R

suppressPackageStartupMessages(library(yaml))

# ---- Configuration -----------------------------------------------------------

ROOT <- getwd()  # must be run from project root

# Verified baseline run directories (relative to ROOT).
# Note: older runs use "DATASE_" typo (missing T) in folder names.
PIPELINEA_RUNS <- list(
  list(dataset = "dataset_0", dir = "results/pipelineA/DATASE_0_BASELINE_2026-03-01_18-00-42"),
  list(dataset = "dataset_1", dir = "results/pipelineA/DATASE_1_BASELINE_2026-03-03_15-00-39"),
  list(dataset = "dataset_2", dir = "results/pipelineA/DATASET_2_BASELINE_2026-03-04_11-29-07"),
  list(dataset = "dataset_3", dir = "results/pipelineA/DATASET_3_BASELINE_2026-03-04_14-54-56")
)

PIPELINEB_RUNS <- list(
  list(dataset = "dataset_0", dir = "results/pipelineB/DATASE_0_BASELINE_2026-03-01_18-04-29"),
  list(dataset = "dataset_1", dir = "results/pipelineB/DATASE_1_BASELINE_2026-03-04_09-54-24"),
  list(dataset = "dataset_2", dir = "results/pipelineB/DATASET_2_BASELINE_2026-03-04_14-43-36"),
  list(dataset = "dataset_3", dir = "results/pipelineB/DATASET_3_BASELINE_2026-03-04_12-32-56")
)

OUT_DIR  <- file.path(ROOT, "results", "analysis")
OUT_FILE <- file.path(OUT_DIR, "dataset_summary.csv")

# ---- Helper functions --------------------------------------------------------

# Count data rows in a CSV (total lines minus 1 header row).
# Returns NA if file is missing; 0 if file has only a header.
csv_nrow <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  n <- length(readLines(path, warn = FALSE))
  max(0L, as.integer(n) - 1L)
}

# Count lines in a plain text file (no header assumed; one item per line).
txt_nlines <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]  # drop blank lines
  length(lines)
}

# Read a single column from a CSV and return its mean.
# Returns NA if file is missing or column not found.
csv_col_mean <- function(path, col) {
  if (!file.exists(path)) return(NA_real_)
  df <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !col %in% names(df)) return(NA_real_)
  mean(as.numeric(df[[col]]), na.rm = TRUE)
}

# Count rows in a CSV where a logical/character column equals TRUE / "TRUE".
csv_count_true <- function(path, col) {
  if (!file.exists(path)) return(NA_integer_)
  df <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !col %in% names(df)) return(NA_integer_)
  sum(as.character(df[[col]]) == "TRUE", na.rm = TRUE)
}

# Read the unique values of one column from a CSV, collapsed to a string.
# Returns NA if file or column missing.
csv_col_unique_str <- function(path, col) {
  if (!file.exists(path)) return(NA_character_)
  df <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !col %in% names(df)) return(NA_character_)
  vals <- unique(trimws(as.character(df[[col]])))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0) return(NA_character_)
  paste(vals, collapse = "; ")
}

# Parse Pipeline A final_summary.txt for FINAL source and n_final.
# Returns a named list: final_source (character).
parse_A_summary <- function(path) {
  out <- list(final_source = NA_character_)
  if (!file.exists(path)) return(out)
  lines <- readLines(path, warn = FALSE)

  # dataset_0 format: "Pathways (Tier 1):  11" (no FINAL source line)
  tier1_line <- grep("^Pathways \\(Tier 1\\):", lines, value = TRUE)
  if (length(tier1_line) > 0) {
    out$final_source <- "tier1"
    return(out)
  }

  # datasets 1-3 format: "FINAL source:  Tier 2 ..."
  src_line <- grep("^FINAL source:", lines, value = TRUE)
  if (length(src_line) > 0) {
    raw <- trimws(sub("^FINAL source:", "", src_line[1]))
    # Normalise to short tag.
    # Use ^ anchor so "Tier 2 (Tier 1 empty)" does NOT match "tier1".
    if (grepl("^Tier 1", raw))             out$final_source <- "tier1"
    else if (grepl("^Tier 2", raw))        out$final_source <- "tier2"
    else if (grepl("(?i)candidate", raw))  out$final_source <- "candidates"
    else                                   out$final_source <- raw
  }
  out
}

# Parse Pipeline B final_summary.txt for final_source and agreement_source.
# Returns a named list.
parse_B_summary <- function(path) {
  out <- list(final_source = NA_character_, agreement_source_summary = NA_character_)
  if (!file.exists(path)) return(out)
  lines <- readLines(path, warn = FALSE)

  # --- final_source: read "Final selection source used:" block
  idx <- grep("Final selection source used:", lines)
  if (length(idx) > 0) {
    # Next lines are "  UP:  xxx" and "  DOWN:  yyy"
    block <- lines[seq(idx[1] + 1, min(idx[1] + 3, length(lines)))]
    up_line   <- grep("UP:", block, value = TRUE)
    down_line <- grep("DOWN:", block, value = TRUE)
    up_src   <- if (length(up_line))   trimws(sub(".*UP:", "", up_line[1]))   else NA_character_
    down_src <- if (length(down_line)) trimws(sub(".*DOWN:", "", down_line[1])) else NA_character_
    if (!is.na(up_src) && !is.na(down_src) && up_src == down_src) {
      out$final_source <- up_src
    } else {
      out$final_source <- paste0("UP:", up_src, "/DOWN:", down_src)
    }
  }

  # --- agreement_source: from "Step03 agreement:" block (only newer runs)
  wp_line <- grep("WP source:", lines, value = TRUE)
  go_line <- grep("GO source:", lines, value = TRUE)
  if (length(wp_line) > 0 || length(go_line) > 0) {
    wp_src <- if (length(wp_line)) trimws(sub(".*WP source:", "", wp_line[1])) else NA_character_
    go_src <- if (length(go_line)) trimws(sub(".*GO source:", "", go_line[1])) else NA_character_
    # Extract just the source tag before the parenthetical, e.g. "strict_intersection (14 pathways)" → "strict_intersection"
    clean <- function(s) if (is.na(s)) NA_character_ else trimws(sub("\\s*\\(.*", "", s))
    wp_tag <- clean(wp_src)
    go_tag <- clean(go_src)
    if (!is.na(wp_tag) && !is.na(go_tag) && wp_tag == go_tag) {
      out$agreement_source_summary <- wp_tag
    } else {
      out$agreement_source_summary <- paste0("WP:", wp_tag, "/GO:", go_tag)
    }
  }
  out
}

# ---- Extract one Pipeline A row ---------------------------------------------

extract_pipelineA <- function(dataset, rel_dir) {
  dir    <- file.path(ROOT, rel_dir)
  notes  <- character(0)

  # --- Universe + DE (from de_policy.yml when available) ---
  policy_path <- file.path(dir, "02_de", "de_policy.yml")
  policy      <- NULL
  if (file.exists(policy_path)) {
    policy <- tryCatch(yaml::read_yaml(policy_path), error = function(e) NULL)
  } else {
    notes <- c(notes, "de_policy.yml missing (older run; DE fields derived from text files)")
  }

  n_universe <- if (!is.null(policy$universe_N))
    as.integer(policy$universe_N)
  else
    txt_nlines(file.path(dir, "01_universe", "universe_genes.txt"))

  n_de <- if (!is.null(policy$de_count_n))
    as.integer(policy$de_count_n)
  else
    txt_nlines(file.path(dir, "02_de", "de_genes_rawp05.txt"))

  de_fraction <- if (!is.null(policy$de_fraction))
    round(as.numeric(policy$de_fraction), 4)
  else if (!is.na(n_de) && !is.na(n_universe) && n_universe > 0)
    round(n_de / n_universe, 4)
  else
    NA_real_

  # Ladder steps: list may be empty (no tightening) or have items like "L1: |log2FC| >= 0.5"
  de_ladder <- if (!is.null(policy)) {
    steps <- policy$tightening_steps
    if (is.null(steps) || length(steps) == 0) "none" else paste(steps, collapse = "; ")
  } else {
    NA_character_
  }

  # --- Tier counts ---
  n_tier1      <- csv_nrow(file.path(dir, "07_tiers", "tier1_final.csv"))
  n_tier2      <- csv_nrow(file.path(dir, "07_tiers", "tier2_stable_only.csv"))
  n_candidates <- csv_nrow(file.path(dir, "07_tiers", "tiered_all_candidates.csv"))

  # --- Final ---
  n_final      <- csv_nrow(file.path(dir, "FINAL", "FINAL.csv"))
  summary_info <- parse_A_summary(file.path(dir, "FINAL", "final_summary.txt"))
  final_source <- summary_info$final_source

  # --- Stability (mean across all representative pathways) ---
  mean_stab_go <- csv_col_mean(
    file.path(dir, "06_bootstrap_consensus", "go_bp_stability.csv"), "stability"
  )
  mean_stab_wp <- csv_col_mean(
    file.path(dir, "06_bootstrap_consensus", "wp_stability.csv"), "stability"
  )

  # --- Cross-collection consensus pairs ---
  consensus_pairs <- csv_nrow(
    file.path(dir, "06_bootstrap_consensus", "cross_collection_pairs.csv")
  )

  data.frame(
    dataset                   = dataset,
    pipeline                  = "pipelineA",
    run_dir                   = rel_dir,
    n_universe                = n_universe,
    n_de                      = n_de,
    de_fraction               = de_fraction,
    de_policy_level_or_ladder = de_ladder,
    n_tier1                   = n_tier1,
    n_tier2                   = n_tier2,
    n_candidates              = n_candidates,
    n_final                   = n_final,
    final_source              = final_source,
    agreement_source          = NA_character_,   # N/A for Pipeline A
    mean_stability_go         = round(mean_stab_go, 4),
    mean_stability_wp         = round(mean_stab_wp, 4),
    consensus_pairs           = consensus_pairs,
    notes                     = paste(notes, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

# ---- Extract one Pipeline B row ---------------------------------------------

extract_pipelineB <- function(dataset, rel_dir) {
  dir   <- file.path(ROOT, rel_dir)
  notes <- character(0)

  # --- Universe (count lines; no header in universe_genes.txt) ---
  n_universe <- txt_nlines(file.path(dir, "01_camera", "universe_genes.txt"))

  # --- Final counts ---
  n_final      <- csv_nrow(file.path(dir, "FINAL", "FINAL.csv"))
  summary_info <- parse_B_summary(file.path(dir, "FINAL", "final_summary.txt"))
  final_source <- summary_info$final_source

  # --- Agreement source ---
  # Available only in newer runs (dataset_2, dataset_3).
  # Older runs (dataset_0, dataset_1) lack the agreement_source column in step03 CSVs.
  agreement_src <- csv_col_unique_str(
    file.path(dir, "03_prepare_inputs", "wp_input_agreement.csv"),
    "agreement_source"
  )
  if (is.na(agreement_src)) {
    # Try GO file as fallback
    agreement_src <- csv_col_unique_str(
      file.path(dir, "03_prepare_inputs", "go_bp_input_agreement.csv"),
      "agreement_source"
    )
  }
  if (is.na(agreement_src)) {
    # Use parsed value from summary.txt if available (newer runs have "Step03 agreement" block)
    agreement_src <- summary_info$agreement_source_summary
  }
  if (is.na(agreement_src)) {
    notes <- c(notes, "agreement_source not available (older run format)")
  }

  # --- Consensus pairs (count in_consensus==TRUE across UP + DOWN step06 outputs) ---
  n_cons_up   <- csv_count_true(file.path(dir, "06_consensus", "final_up.csv"),   "in_consensus")
  n_cons_down <- csv_count_true(file.path(dir, "06_consensus", "final_down.csv"), "in_consensus")
  consensus_pairs <- if (!is.na(n_cons_up) && !is.na(n_cons_down)) {
    n_cons_up + n_cons_down
  } else {
    NA_integer_
  }

  data.frame(
    dataset                   = dataset,
    pipeline                  = "pipelineB",
    run_dir                   = rel_dir,
    n_universe                = n_universe,
    n_de                      = NA_integer_,    # Pipeline B uses all ranked genes
    de_fraction               = NA_real_,
    de_policy_level_or_ladder = NA_character_,
    n_tier1                   = NA_integer_,    # Pipeline B has no tier system
    n_tier2                   = NA_integer_,
    n_candidates              = NA_integer_,
    n_final                   = n_final,
    final_source              = final_source,
    agreement_source          = agreement_src,
    mean_stability_go         = NA_real_,       # Pipeline B has no bootstrap stability
    mean_stability_wp         = NA_real_,
    consensus_pairs           = consensus_pairs,
    notes                     = paste(notes, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

# ---- Main: collect all rows --------------------------------------------------

cat("=== 00_diagnostics.R: scanning baseline results ===\n\n")

rows      <- list()
warnings_ <- character(0)

for (entry in PIPELINEA_RUNS) {
  full_dir <- file.path(ROOT, entry$dir)
  if (!dir.exists(full_dir)) {
    msg <- sprintf("MISSING: %s (skipping)", entry$dir)
    cat(msg, "\n")
    warnings_ <- c(warnings_, msg)
    next
  }
  cat(sprintf("Pipeline A | %s ...\n", entry$dataset))
  row <- tryCatch(
    extract_pipelineA(entry$dataset, entry$dir),
    error = function(e) {
      msg <- sprintf("ERROR in %s pipelineA: %s", entry$dataset, conditionMessage(e))
      cat(msg, "\n")
      warnings_ <<- c(warnings_, msg)
      NULL
    }
  )
  if (!is.null(row)) rows[[length(rows) + 1]] <- row
}

for (entry in PIPELINEB_RUNS) {
  full_dir <- file.path(ROOT, entry$dir)
  if (!dir.exists(full_dir)) {
    msg <- sprintf("MISSING: %s (skipping)", entry$dir)
    cat(msg, "\n")
    warnings_ <- c(warnings_, msg)
    next
  }
  cat(sprintf("Pipeline B | %s ...\n", entry$dataset))
  row <- tryCatch(
    extract_pipelineB(entry$dataset, entry$dir),
    error = function(e) {
      msg <- sprintf("ERROR in %s pipelineB: %s", entry$dataset, conditionMessage(e))
      cat(msg, "\n")
      warnings_ <<- c(warnings_, msg)
      NULL
    }
  )
  if (!is.null(row)) rows[[length(rows) + 1]] <- row
}

if (length(rows) == 0) stop("No rows collected — check run directories.")

summary_df <- do.call(rbind, rows)
rownames(summary_df) <- NULL

# ---- Write output ------------------------------------------------------------

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, OUT_FILE, row.names = FALSE)
cat(sprintf("\nWrote: %s\n", OUT_FILE))

# ---- Print summary -----------------------------------------------------------

cat("\n")
cat("=== SUMMARY ===\n")
cat(sprintf("Rows collected: %d / %d expected\n",
            nrow(summary_df),
            length(PIPELINEA_RUNS) + length(PIPELINEB_RUNS)))
cat("\n")

# Print a readable table
fmt_val <- function(x) if (is.na(x)) "NA" else as.character(x)

cat(sprintf("%-12s %-12s %8s %6s %7s %-10s %7s %8s %8s %-30s %-26s\n",
            "dataset", "pipeline", "n_univ", "n_de", "de_frac",
            "de_ladder", "n_final", "stab_go", "stab_wp",
            "final_source", "agreement_source"))
cat(strrep("-", 130), "\n")

for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%-12s %-12s %8s %6s %7s %-10s %7s %8s %8s %-30s %-26s\n",
              r$dataset,
              r$pipeline,
              fmt_val(r$n_universe),
              fmt_val(r$n_de),
              fmt_val(r$de_fraction),
              fmt_val(r$de_policy_level_or_ladder),
              fmt_val(r$n_final),
              fmt_val(r$mean_stability_go),
              fmt_val(r$mean_stability_wp),
              fmt_val(r$final_source),
              fmt_val(r$agreement_source)
  ))
}

cat("\n")
cat("Tier / consensus counts (Pipeline A only):\n")
cat(sprintf("%-12s %8s %8s %12s %16s\n",
            "dataset", "n_tier1", "n_tier2", "n_candidates", "consensus_pairs"))
cat(strrep("-", 62), "\n")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  if (r$pipeline == "pipelineA") {
    cat(sprintf("%-12s %8s %8s %12s %16s\n",
                r$dataset,
                fmt_val(r$n_tier1),
                fmt_val(r$n_tier2),
                fmt_val(r$n_candidates),
                fmt_val(r$consensus_pairs)
    ))
  }
}

cat("\n")
if (length(warnings_) > 0) {
  cat("Warnings / missing files:\n")
  for (w in warnings_) cat(" -", w, "\n")
} else {
  cat("No warnings.\n")
}

# Report notes from individual rows
notes_present <- summary_df[nzchar(trimws(summary_df$notes)), c("dataset", "pipeline", "notes")]
if (nrow(notes_present) > 0) {
  cat("\nPer-row notes:\n")
  for (i in seq_len(nrow(notes_present))) {
    cat(sprintf("  [%s / %s] %s\n",
                notes_present$dataset[i], notes_present$pipeline[i], notes_present$notes[i]))
  }
}

cat("\nDone.\n")
