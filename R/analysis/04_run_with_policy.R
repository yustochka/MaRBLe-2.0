# R/analysis/04_run_with_policy.R
#
# Generalized threshold automation for Pipeline A.
# Works for any dataset listed in results/run_registry.csv — no code changes needed
# to add new datasets.
#
# Classification and threshold selection are computed entirely from data
# (bootstrap matrices, universe, gene sets) — no reliance on precomputed
# threshold_policy.csv or hardcoded dataset names.
#
# Usage:
#   Rscript R/analysis/04_run_with_policy.R --dataset dataset_0
#   Rscript R/analysis/04_run_with_policy.R --dataset dataset_3 --tau 0.70 --cons_j 0.10
#
# To add a new dataset: append a row to results/run_registry.csv.
#
# Required inputs (pre-existing, never regenerated):
#   results/run_registry.csv                            <- dataset -> run dir mapping
#   <run_dir>/06_bootstrap_consensus/go_bp_boot_matrix.csv
#   <run_dir>/06_bootstrap_consensus/wp_boot_matrix.csv
#   <run_dir>/06_bootstrap_consensus/go_bp_stability.csv
#   <run_dir>/06_bootstrap_consensus/wp_stability.csv
#   <run_dir>/05_overlap/go_bp/representatives.csv
#   <run_dir>/05_overlap/wp/representatives.csv
#   <run_dir>/01_universe/universe_genes.txt
#   <run_dir>/07_tiers/tiered_all_candidates.csv
#   <wp_file>                                           <- from registry
#
# Optional (used for LOW_COVERAGE frac_in_matrix check only):
#   results/analysis/expressed_fraction.csv
#
# Outputs:
#   results/analysis/policy_runs/<dataset>/FINAL.csv
#   results/analysis/policy_runs/<dataset>/run_summary.txt
#   results/analysis/policy_runs/<dataset>/comparison.csv
#
# Heavy steps NOT rerun: ORA (step03), bootstrap sampling (step06 loop),
#   GO semantic collapse (step04), Jaccard clustering (step05).
# Pairwise GO x WP Jaccard IS recomputed (fast; needed for CONS_J < 0.15).

suppressPackageStartupMessages({
  library(readr)
  library(org.Hs.eg.db)   # load before anything that might mask select()
  library(AnnotationDbi)
  # dplyr is NOT loaded globally: it masks AnnotationDbi::select()
})

ROOT <- getwd()

# =============================================================================
# Constants
# =============================================================================

CONS_K_MIN          <- 5L
LOW_UNIVERSE        <- 2000L
LOW_COVERAGE_FRAC   <- 0.20
BASELINE_TAU        <- 0.70
BASELINE_CONS_J     <- 0.15
SWEEP_TAU_VALUES    <- c(0.50, 0.60, 0.70)
SWEEP_CONS_J_VALUES <- c(0.10, 0.15, 0.20)

# Threshold rules per classification (BORDERLINE cons_j is data-driven)
THRESHOLD_RULES <- list(
  HEALTHY      = list(tau = 0.70, cons_j = 0.15),
  NO_OVERLAP   = list(tau = 0.70, cons_j = 0.10),
  LOW_COVERAGE = list(tau = 0.50, cons_j = 0.10),
  BORDERLINE   = list(tau = 0.70, cons_j = NULL)   # cons_j filled from sweep
)

# =============================================================================
# Registry
# =============================================================================

load_run_registry <- function(registry_path) {
  if (!file.exists(registry_path))
    stop("Run registry not found: ", registry_path,
         "\nCreate it at results/run_registry.csv with columns:",
         " dataset, pipeline_a_dir, wp_file")
  df <- read_csv(registry_path, show_col_types = FALSE)
  required <- c("dataset", "pipeline_a_dir", "wp_file")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0)
    stop("Run registry missing columns: ", paste(missing, collapse = ", "))
  df
}

# =============================================================================
# Data loading helpers
# =============================================================================

# Bootstrap matrix -> named stability vector (colMeans).
# readr preserves "GO:xxxxx" column names; base read.csv mangles them.
load_stability <- function(path) {
  if (!file.exists(path)) stop("Bootstrap matrix not found: ", path)
  mat_raw <- read_csv(path, show_col_types = FALSE)
  ids  <- names(mat_raw)[-1]
  mat  <- as.matrix(mat_raw[, -1, drop = FALSE])
  stab <- colMeans(mat, na.rm = TRUE)
  names(stab) <- ids
  stab
}

# Stability CSV -> data.frame with ID, Description, p.adjust, k, stability.
load_stability_meta <- function(path) {
  if (!file.exists(path)) stop("Stability CSV not found: ", path)
  df <- read_csv(path, show_col_types = FALSE)
  df$ID          <- as.character(df$ID)
  df$Description <- as.character(df$Description)
  df
}

# Plain gene-per-line text file -> character vector.
read_genes <- function(path) {
  if (!file.exists(path)) stop("Gene list not found: ", path)
  x <- readLines(path, warn = FALSE)
  unique(trimws(x[nzchar(trimws(x))]))
}

# Try to load median frac_in_matrix from expressed_fraction.csv (optional).
# Returns NA_real_ if file absent or dataset has no pipelineA rows.
load_frac_coverage <- function(expr_csv, dataset) {
  if (!file.exists(expr_csv)) {
    message("  [info] expressed_fraction.csv not found",
            " -- LOW_COVERAGE check uses n_universe only")
    return(NA_real_)
  }
  df  <- read_csv(expr_csv, show_col_types = FALSE)
  sub <- df[df$dataset == dataset & df$pipeline == "pipelineA", ]
  if (nrow(sub) == 0) {
    message("  [info] No pipelineA rows for '", dataset,
            "' in expressed_fraction.csv -- using n_universe only")
    return(NA_real_)
  }
  median(sub$frac_in_matrix, na.rm = TRUE)
}

# =============================================================================
# Gene set builders
# =============================================================================

# GO-BP gene sets filtered to universe (base R, no dplyr).
build_go_sets <- function(go_ids, universe) {
  if (length(go_ids) == 0) return(list())
  raw <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = unique(go_ids),
      keytype = "GO",
      columns = c("SYMBOL", "ONTOLOGY")
    )
  )
  id_col <- intersect(c("GO", "GOID", "GOALL"), names(raw))
  if (length(id_col) == 0) stop("org.Hs.eg.db returned no GO id column")
  if (!"ONTOLOGY" %in% names(raw)) raw$ONTOLOGY <- NA_character_
  raw$ID <- as.character(raw[[id_col[1]]])
  raw <- raw[!is.na(raw$ID) & !is.na(raw$SYMBOL), ]
  raw <- raw[is.na(raw$ONTOLOGY) | raw$ONTOLOGY == "BP", ]
  raw <- raw[raw$SYMBOL %in% universe, ]
  raw <- raw[!duplicated(raw[, c("ID", "SYMBOL")]), ]
  gene_list <- split(raw$SYMBOL, raw$ID)
  lapply(gene_list, unique)
}

# WP gene sets filtered to universe (base R, no dplyr).
build_wp_sets <- function(wp_ids, universe, wp_map_path) {
  if (!file.exists(wp_map_path)) stop("Pathway2Gene.csv not found: ", wp_map_path)
  if (length(wp_ids) == 0) return(list())
  raw <- read_csv(wp_map_path, show_col_types = FALSE)
  id_cands  <- c("wpid","ID","id","Pathway","pathway","pathway_id","WPID","WP_ID","WP")
  sym_cands <- c("SYMBOL","symbol","GeneSymbol","gene_symbol","Gene","gene",
                 "HGNC_symbol","HUGO")
  id_col  <- names(raw)[tolower(names(raw)) %in% tolower(id_cands)]
  sym_col <- names(raw)[tolower(names(raw)) %in% tolower(sym_cands)]
  if (!length(id_col) || !length(sym_col))
    stop("Pathway2Gene.csv missing ID/SYMBOL columns. Found: ",
         paste(names(raw), collapse = ", "))
  tbl <- data.frame(
    ID     = as.character(raw[[id_col[1]]]),
    SYMBOL = as.character(raw[[sym_col[1]]]),
    stringsAsFactors = FALSE
  )
  tbl <- tbl[tbl$ID %in% wp_ids & tbl$SYMBOL %in% universe, ]
  tbl <- tbl[!duplicated(tbl), ]
  gene_list <- split(tbl$SYMBOL, tbl$ID)
  lapply(gene_list, unique)
}

# All pairwise GO x WP Jaccard (no threshold filter).
# Returns data.frame: A (GO id), B (WP id), jaccard, k.
compute_all_jaccard <- function(go_sets, wp_sets) {
  if (!length(go_sets) || !length(wp_sets))
    return(data.frame(A=character(), B=character(),
                      jaccard=numeric(), k=integer(), stringsAsFactors=FALSE))
  a_ids <- names(go_sets)
  b_ids <- names(wp_sets)
  rows  <- vector("list", length(a_ids) * length(b_ids))
  idx   <- 0L
  for (a in a_ids) {
    for (b in b_ids) {
      inter <- length(intersect(go_sets[[a]], wp_sets[[b]]))
      uni   <- length(union(go_sets[[a]], wp_sets[[b]]))
      if (uni == 0L) next
      idx <- idx + 1L
      rows[[idx]] <- data.frame(A=a, B=b, jaccard=inter/uni, k=inter,
                                stringsAsFactors=FALSE)
    }
  }
  if (idx == 0L) return(data.frame(A=character(), B=character(),
                                    jaccard=numeric(), k=integer(), stringsAsFactors=FALSE))
  do.call(rbind, rows[seq_len(idx)])
}

# =============================================================================
# Mini-sweep: count-only threshold application for classification
# =============================================================================

apply_thresholds_counts <- function(stability_go, stability_wp, all_pairs,
                                     tau, cons_j, cons_k) {
  go_sel <- names(stability_go)[stability_go >= tau]
  wp_sel <- names(stability_wp)[stability_wp >= tau]
  pf     <- all_pairs[all_pairs$jaccard >= cons_j & all_pairs$k >= cons_k, ,
                       drop = FALSE]
  gc     <- unique(pf$A)
  wc     <- unique(pf$B)
  t1     <- sum(go_sel %in% gc) + sum(wp_sel %in% wc)
  t2     <- length(go_sel) + length(wp_sel) - t1
  src    <- if (t1 >= 5L) "tier1" else if (t2 > 0L) "tier2" else "candidates"
  data.frame(tau=tau, cons_j=cons_j, n_tier1=t1, n_tier2=t2,
             n_consensus_pairs=nrow(pf), final_source=src,
             stringsAsFactors=FALSE)
}

run_mini_sweep <- function(stability_go, stability_wp, all_pairs,
                            tau_vals, cons_j_vals, cons_k) {
  rows <- vector("list", length(tau_vals) * length(cons_j_vals))
  i <- 0L
  for (tau in tau_vals) {
    for (cj in cons_j_vals) {
      i <- i + 1L
      rows[[i]] <- apply_thresholds_counts(stability_go, stability_wp,
                                            all_pairs, tau, cj, cons_k)
    }
  }
  do.call(rbind, rows)
}

# =============================================================================
# Classification and threshold selection
# =============================================================================

classify_dataset <- function(n_universe, med_frac_in_matrix, sweep) {
  # --- LOW_COVERAGE ---
  n_ok    <- !is.na(n_universe) && n_universe >= LOW_UNIVERSE
  frac_ok <- is.na(med_frac_in_matrix) || med_frac_in_matrix >= LOW_COVERAGE_FRAC
  if (!n_ok || !frac_ok) {
    return(list(
      classification        = "LOW_COVERAGE",
      min_cons_j_tier1      = NA_real_,
      max_consensus_pairs   = max(sweep$n_consensus_pairs, na.rm = TRUE),
      baseline_n_tier1      = NA_integer_,
      baseline_final_source = NA_character_
    ))
  }

  max_pairs    <- max(sweep$n_consensus_pairs, na.rm = TRUE)
  tier1_rows   <- sweep[sweep$final_source == "tier1", ]
  tier1_achiev <- nrow(tier1_rows) > 0
  min_cj_tier1 <- if (tier1_achiev) min(tier1_rows$cons_j) else NA_real_

  baseline <- sweep[
    abs(sweep$tau    - BASELINE_TAU)    < 1e-9 &
    abs(sweep$cons_j - BASELINE_CONS_J) < 1e-9, ]
  bl_n_tier1 <- if (nrow(baseline) == 1) baseline$n_tier1     else NA_integer_
  bl_src     <- if (nrow(baseline) == 1) baseline$final_source else NA_character_
  bl_is_tier1 <- isTRUE(bl_src == "tier1")

  # --- NO_OVERLAP ---
  if (!tier1_achiev && max_pairs == 0) {
    return(list(
      classification        = "NO_OVERLAP",
      min_cons_j_tier1      = NA_real_,
      max_consensus_pairs   = max_pairs,
      baseline_n_tier1      = bl_n_tier1,
      baseline_final_source = bl_src
    ))
  }

  # --- BORDERLINE ---
  if (!bl_is_tier1 && tier1_achiev) {
    return(list(
      classification        = "BORDERLINE",
      min_cons_j_tier1      = min_cj_tier1,
      max_consensus_pairs   = max_pairs,
      baseline_n_tier1      = bl_n_tier1,
      baseline_final_source = bl_src
    ))
  }

  # --- HEALTHY ---
  list(
    classification        = "HEALTHY",
    min_cons_j_tier1      = min_cj_tier1,
    max_consensus_pairs   = max_pairs,
    baseline_n_tier1      = bl_n_tier1,
    baseline_final_source = bl_src
  )
}

select_thresholds <- function(classification, min_cons_j_tier1) {
  rule <- THRESHOLD_RULES[[classification]]
  if (is.null(rule)) stop("Unknown classification: ", classification)
  if (classification == "BORDERLINE")
    rule$cons_j <- if (!is.na(min_cons_j_tier1)) min_cons_j_tier1 else 0.10
  rule
}

# =============================================================================
# Full threshold application (returns IDs + counts for final run)
# =============================================================================

apply_thresholds_full <- function(stability_go, stability_wp, all_pairs,
                                   n_candidates_fixed, tau, cons_j, cons_k) {
  go_selected <- names(stability_go)[stability_go >= tau]
  wp_selected <- names(stability_wp)[stability_wp >= tau]
  pairs_filt  <- all_pairs[all_pairs$jaccard >= cons_j & all_pairs$k >= cons_k, ,
                            drop = FALSE]
  go_cons_ids <- unique(pairs_filt$A)
  wp_cons_ids <- unique(pairs_filt$B)
  tier1_go <- go_selected[go_selected %in% go_cons_ids]
  tier1_wp <- wp_selected[wp_selected %in% wp_cons_ids]
  tier2_go <- go_selected[!go_selected %in% go_cons_ids]
  tier2_wp <- wp_selected[!wp_selected %in% wp_cons_ids]
  n_tier1  <- length(tier1_go) + length(tier1_wp)
  n_tier2  <- length(tier2_go) + length(tier2_wp)
  final_source <- if (n_tier1 >= 5L) "tier1" else if (n_tier2 > 0L) "tier2" else "candidates"
  list(
    tier1_ids         = c(tier1_go, tier1_wp),
    tier2_ids         = c(tier2_go, tier2_wp),
    n_tier1           = n_tier1,
    n_tier2           = n_tier2,
    n_consensus_pairs = nrow(pairs_filt),
    final_source      = final_source,
    n_candidates      = n_candidates_fixed
  )
}

# =============================================================================
# Build FINAL.csv data frame
# =============================================================================

attach_meta <- function(ids, go_meta, wp_meta, tau, consensus_flag) {
  if (length(ids) == 0) return(data.frame())
  rows <- lapply(ids, function(id) {
    is_go <- startsWith(id, "GO:")
    meta  <- if (is_go) go_meta[go_meta$ID == id, ] else wp_meta[wp_meta$ID == id, ]
    if (nrow(meta) == 0) {
      message("  WARNING: no metadata found for ", id, " -- skipping")
      return(NULL)
    }
    data.frame(
      ID          = id,
      Description = as.character(meta$Description[1]),
      p.adjust    = meta$p.adjust[1],
      k           = as.integer(meta$k[1]),
      stability   = meta$stability[1],
      selected    = meta$stability[1] >= tau,
      collection  = if (is_go) "GO" else "WP",
      consensus   = consensus_flag,
      keep_final  = TRUE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

build_final_df <- function(result, go_meta, wp_meta, tiered_candidates, tau) {
  if (result$final_source == "tier1") {
    attach_meta(result$tier1_ids, go_meta, wp_meta, tau, consensus_flag = TRUE)
  } else if (result$final_source == "tier2") {
    attach_meta(result$tier2_ids, go_meta, wp_meta, tau, consensus_flag = FALSE)
  } else {
    cands <- tiered_candidates
    cands$selected   <- cands$stability >= tau
    cands$consensus  <- FALSE
    cands$keep_final <- TRUE
    keep_cols <- c("ID","Description","p.adjust","k","stability",
                   "selected","collection","consensus","keep_final")
    cands[, intersect(keep_cols, names(cands)), drop = FALSE]
  }
}

# =============================================================================
# Output writers
# =============================================================================

write_run_summary <- function(out_dir, dataset, classification, signals,
                               tau, cons_j, override_used, result,
                               n_universe, med_frac_in_matrix) {
  confidence <- switch(classification,
    HEALTHY      = "high",
    BORDERLINE   = "medium",
    NO_OVERLAP   = "low",
    LOW_COVERAGE = "unreliable"
  )
  frac_str <- if (!is.na(med_frac_in_matrix)) sprintf("%.3f", med_frac_in_matrix) else "NA (not computed)"
  txt <- paste0(
    "=== run_summary.txt ===\n",
    sprintf("dataset:                    %s\n",  dataset),
    sprintf("n_universe:                 %d\n",  n_universe),
    sprintf("med_frac_in_matrix:         %s\n",  frac_str),
    sprintf("max_consensus_pairs:        %d\n",  signals$max_consensus_pairs),
    sprintf("baseline_n_tier1:           %s\n",
            if (!is.na(signals$baseline_n_tier1)) as.character(signals$baseline_n_tier1) else "NA"),
    sprintf("min_cons_j_achieving_tier1: %s\n",
            if (!is.na(signals$min_cons_j_tier1)) as.character(signals$min_cons_j_tier1) else "NA"),
    sprintf("classification:             %s\n",  classification),
    sprintf("confidence_flag:            %s\n",  confidence),
    sprintf("TAU_STABILITY:              %.2f\n", tau),
    sprintf("CONS_JACCARD_MIN:           %.2f\n", cons_j),
    sprintf("CONS_K_MIN:                 %d\n",  CONS_K_MIN),
    sprintf("overrides_used:             %s\n",
            if (override_used) "YES" else "NO (rule-based)"),
    sprintf("final_source:               %s\n",  result$final_source),
    sprintf("n_tier1:                    %d\n",  result$n_tier1),
    sprintf("n_tier2:                    %d\n",  result$n_tier2),
    sprintf("n_consensus_pairs:          %d\n",  result$n_consensus_pairs),
    sprintf("n_candidates_pool:          %s\n",
            if (!is.na(result$n_candidates)) as.character(result$n_candidates) else "NA")
  )
  writeLines(txt, file.path(out_dir, "run_summary.txt"))
  cat(txt)
}

write_comparison <- function(out_dir, dataset, tau_rec, cons_j_rec,
                              result_rec, mini_sweep) {
  baseline <- mini_sweep[
    abs(mini_sweep$tau    - BASELINE_TAU)    < 1e-9 &
    abs(mini_sweep$cons_j - BASELINE_CONS_J) < 1e-9, ]

  comp <- data.frame(
    dataset                    = dataset,
    baseline_TAU               = BASELINE_TAU,
    baseline_CONS_J            = BASELINE_CONS_J,
    baseline_n_tier1           = if (nrow(baseline)==1) baseline$n_tier1           else NA_integer_,
    baseline_n_tier2           = if (nrow(baseline)==1) baseline$n_tier2           else NA_integer_,
    baseline_consensus_pairs   = if (nrow(baseline)==1) baseline$n_consensus_pairs else NA_integer_,
    baseline_final_source      = if (nrow(baseline)==1) baseline$final_source      else NA_character_,
    recommended_TAU            = tau_rec,
    recommended_CONS_J         = cons_j_rec,
    recommended_n_tier1        = result_rec$n_tier1,
    recommended_n_tier2        = result_rec$n_tier2,
    recommended_consensus_pairs = result_rec$n_consensus_pairs,
    recommended_final_source   = result_rec$final_source,
    stringsAsFactors = FALSE
  )
  write_csv(comp, file.path(out_dir, "comparison.csv"))
  cat(sprintf("  Wrote: %s\n", file.path(out_dir, "comparison.csv")))

  # Print brief comparison table
  cat("\n  --- Comparison: baseline vs recommended ---\n")
  cat(sprintf("  %26s  %12s  %12s\n", "", "baseline", "recommended"))
  bl_t1  <- if (nrow(baseline)==1) baseline$n_tier1           else NA
  bl_t2  <- if (nrow(baseline)==1) baseline$n_tier2           else NA
  bl_cp  <- if (nrow(baseline)==1) baseline$n_consensus_pairs else NA
  bl_src <- if (nrow(baseline)==1) baseline$final_source      else "NA"
  cat(sprintf("  %26s  %12.2f  %12.2f\n", "TAU_STABILITY",    BASELINE_TAU,    tau_rec))
  cat(sprintf("  %26s  %12.2f  %12.2f\n", "CONS_JACCARD_MIN", BASELINE_CONS_J, cons_j_rec))
  cat(sprintf("  %26s  %12s  %12d\n", "n_tier1",
              if (!is.na(bl_t1)) as.character(bl_t1) else "NA", result_rec$n_tier1))
  cat(sprintf("  %26s  %12s  %12d\n", "n_tier2",
              if (!is.na(bl_t2)) as.character(bl_t2) else "NA", result_rec$n_tier2))
  cat(sprintf("  %26s  %12s  %12d\n", "cons_pairs",
              if (!is.na(bl_cp)) as.character(bl_cp) else "NA", result_rec$n_consensus_pairs))
  cat(sprintf("  %26s  %12s  %12s\n", "final_source", bl_src, result_rec$final_source))
  cat("\n")
}

# =============================================================================
# CLI argument parsing
# =============================================================================

parse_args <- function() {
  args   <- commandArgs(trailingOnly = TRUE)
  result <- list(dataset = NULL, tau = NULL, cons_j = NULL)
  i <- 1L
  while (i <= length(args)) {
    switch(args[i],
      "--dataset" = { result$dataset <- args[i+1L]; i <- i + 2L },
      "--tau"     = { result$tau     <- as.numeric(args[i+1L]); i <- i + 2L },
      "--cons_j"  = { result$cons_j  <- as.numeric(args[i+1L]); i <- i + 2L },
      stop("Unknown argument: ", args[i])
    )
  }
  if (is.null(result$dataset) || !nzchar(result$dataset))
    stop("--dataset is required.\nUsage: Rscript 04_run_with_policy.R --dataset <name>")
  result
}

# =============================================================================
# Main
# =============================================================================

args    <- parse_args()
dataset <- args$dataset

cat("=== 04_run_with_policy.R ===\n")
cat(sprintf("dataset: %s\n\n", dataset))

# ---- Load run registry -------------------------------------------------------
registry_path <- file.path(ROOT, "results", "run_registry.csv")
registry      <- load_run_registry(registry_path)

entry <- registry[registry$dataset == dataset, ]
if (nrow(entry) == 0)
  stop("Dataset '", dataset, "' not found in run registry: ", registry_path,
       "\nKnown datasets: ", paste(registry$dataset, collapse = ", "))

run_dir <- file.path(ROOT, entry$pipeline_a_dir[1])
wp_file <- file.path(ROOT, entry$wp_file[1])

if (!dir.exists(run_dir))
  stop("Run directory not found: ", run_dir)

# ---- Load bootstrap matrices -> stability ------------------------------------
cat("Loading bootstrap matrices...\n")
go_stab <- load_stability(
  file.path(run_dir, "06_bootstrap_consensus", "go_bp_boot_matrix.csv"))
wp_stab <- load_stability(
  file.path(run_dir, "06_bootstrap_consensus", "wp_boot_matrix.csv"))
cat(sprintf("  GO terms: %d  |  WP terms: %d\n", length(go_stab), length(wp_stab)))

go_meta <- load_stability_meta(
  file.path(run_dir, "06_bootstrap_consensus", "go_bp_stability.csv"))
wp_meta <- load_stability_meta(
  file.path(run_dir, "06_bootstrap_consensus", "wp_stability.csv"))

# ---- Restrict stability to step05 representatives ----------------------------
go_reps_path <- file.path(run_dir, "05_overlap", "go_bp", "representatives.csv")
wp_reps_path <- file.path(run_dir, "05_overlap", "wp",    "representatives.csv")
if (!file.exists(go_reps_path)) stop("GO representatives.csv not found: ", go_reps_path)
if (!file.exists(wp_reps_path)) stop("WP representatives.csv not found: ", wp_reps_path)

go_rep_ids <- as.character(read_csv(go_reps_path, show_col_types = FALSE)$ID)
wp_rep_ids <- as.character(read_csv(wp_reps_path, show_col_types = FALSE)$ID)

go_stab <- go_stab[names(go_stab) %in% go_rep_ids]
wp_stab <- wp_stab[names(wp_stab) %in% wp_rep_ids]
go_meta <- go_meta[go_meta$ID %in% go_rep_ids, ]
wp_meta <- wp_meta[wp_meta$ID %in% wp_rep_ids, ]
cat(sprintf("  After rep filter: GO=%d  WP=%d\n", length(go_stab), length(wp_stab)))

# ---- Load universe and candidates pool ---------------------------------------
universe   <- read_genes(file.path(run_dir, "01_universe", "universe_genes.txt"))
n_universe <- length(universe)
cat(sprintf("  Universe: %d genes\n", n_universe))

cands_path         <- file.path(run_dir, "07_tiers", "tiered_all_candidates.csv")
if (!file.exists(cands_path)) stop("tiered_all_candidates.csv not found: ", cands_path)
tiered_candidates  <- read_csv(cands_path, show_col_types = FALSE)
n_candidates_fixed <- nrow(tiered_candidates)
cat(sprintf("  Candidates pool: %d\n\n", n_candidates_fixed))

# ---- Build gene sets -> compute pairwise Jaccard -----------------------------
cat("Building gene sets...\n")
go_sets <- tryCatch(
  build_go_sets(names(go_stab), universe),
  error = function(e) { message("  ERROR building GO sets: ", e$message); list() }
)
wp_sets <- tryCatch(
  build_wp_sets(names(wp_stab), universe, wp_file),
  error = function(e) { message("  ERROR building WP sets: ", e$message); list() }
)
cat(sprintf("  GO: %d sets  |  WP: %d sets\n", length(go_sets), length(wp_sets)))

# Drop reps without gene sets (mirrors step06 behaviour)
go_stab <- go_stab[names(go_stab) %in% names(go_sets)]
wp_stab <- wp_stab[names(wp_stab) %in% names(wp_sets)]

cat("Computing all pairwise Jaccard (GO x WP)...\n")
all_pairs <- compute_all_jaccard(go_sets, wp_sets)
cat(sprintf("  Total GO x WP pairs: %d\n\n", nrow(all_pairs)))

# ---- Optional: frac_in_matrix for LOW_COVERAGE detection ---------------------
expr_csv <- file.path(ROOT, "results", "analysis", "expressed_fraction.csv")
med_frac_in_matrix <- load_frac_coverage(expr_csv, dataset)
if (!is.na(med_frac_in_matrix))
  cat(sprintf("  med_frac_in_matrix (pipelineA): %.3f\n\n", med_frac_in_matrix))

# ---- Mini-sweep for classification -------------------------------------------
cat("Running mini-sweep for classification...\n")
mini_sweep <- run_mini_sweep(go_stab, wp_stab, all_pairs,
                              SWEEP_TAU_VALUES, SWEEP_CONS_J_VALUES, CONS_K_MIN)

cat(sprintf("  %-5s  %-6s  %7s  %7s  %8s  %12s\n",
            "TAU", "CONS_J", "n_tier1", "n_tier2", "n_pairs", "final_source"))
for (i in seq_len(nrow(mini_sweep))) {
  r      <- mini_sweep[i, ]
  marker <- if (abs(r$tau - BASELINE_TAU) < 1e-9 &
                  abs(r$cons_j - BASELINE_CONS_J) < 1e-9) " <- baseline" else ""
  cat(sprintf("  %-5.2f  %-6.2f  %7d  %7d  %8d  %12s%s\n",
              r$tau, r$cons_j, r$n_tier1, r$n_tier2, r$n_consensus_pairs,
              r$final_source, marker))
}
cat("\n")

# ---- Classify ----------------------------------------------------------------
signals        <- classify_dataset(n_universe, med_frac_in_matrix, mini_sweep)
classification <- signals$classification
cat(sprintf("Classification: %s\n", classification))
cat(sprintf("  n_universe=%d  med_frac=%s  max_pairs=%d  tier1_achievable=%s\n\n",
            n_universe,
            if (!is.na(med_frac_in_matrix)) sprintf("%.3f", med_frac_in_matrix) else "NA",
            signals$max_consensus_pairs,
            if (!is.na(signals$min_cons_j_tier1)) "YES" else "NO"))

# ---- Select thresholds (rule-based; CLI overrides take precedence) -----------
thresholds <- select_thresholds(classification, signals$min_cons_j_tier1)
rec_tau    <- thresholds$tau
rec_cons_j <- thresholds$cons_j

override_used <- FALSE
if (!is.null(args$tau)) {
  cat(sprintf("  Override TAU: %.2f -> %.2f\n", rec_tau, args$tau))
  rec_tau <- args$tau;  override_used <- TRUE
}
if (!is.null(args$cons_j)) {
  cat(sprintf("  Override CONS_J: %.2f -> %.2f\n", rec_cons_j, args$cons_j))
  rec_cons_j <- args$cons_j;  override_used <- TRUE
}
cat(sprintf("  TAU_STABILITY:    %.2f%s\n", rec_tau,
            if (!override_used) " (rule-based)" else " (override)"))
cat(sprintf("  CONS_JACCARD_MIN: %.2f%s\n\n", rec_cons_j,
            if (!override_used) " (rule-based)" else " (override)"))

# ---- Apply chosen thresholds -------------------------------------------------
cat(sprintf("Applying thresholds: TAU=%.2f, CONS_J=%.2f, CONS_K_MIN=%d\n",
            rec_tau, rec_cons_j, CONS_K_MIN))
result <- apply_thresholds_full(go_stab, wp_stab, all_pairs,
                                 n_candidates_fixed,
                                 tau    = rec_tau,
                                 cons_j = rec_cons_j,
                                 cons_k = CONS_K_MIN)
cat(sprintf("  n_stable_go:       %d\n", sum(go_stab >= rec_tau)))
cat(sprintf("  n_stable_wp:       %d\n", sum(wp_stab >= rec_tau)))
cat(sprintf("  n_consensus_pairs: %d\n", result$n_consensus_pairs))
cat(sprintf("  n_tier1:           %d\n", result$n_tier1))
cat(sprintf("  n_tier2:           %d\n", result$n_tier2))
cat(sprintf("  final_source:      %s\n\n", result$final_source))

# ---- Build FINAL data frame --------------------------------------------------
final_df <- build_final_df(result, go_meta, wp_meta, tiered_candidates, rec_tau)
if (is.null(final_df) || nrow(final_df) == 0) {
  warning("FINAL.csv will be empty -- no pathways selected at the applied thresholds.")
  final_df <- data.frame(
    ID=character(), Description=character(), p.adjust=numeric(), k=integer(),
    stability=numeric(), selected=logical(), collection=character(),
    consensus=logical(), keep_final=logical(), stringsAsFactors=FALSE
  )
}
cat(sprintf("FINAL shortlist: %d pathways (%s)\n\n", nrow(final_df), result$final_source))

# ---- Write outputs -----------------------------------------------------------
out_dir <- file.path(ROOT, "results", "analysis", "policy_runs", dataset)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_csv(final_df, file.path(out_dir, "FINAL.csv"))
cat(sprintf("Wrote: %s\n\n", file.path(out_dir, "FINAL.csv")))

write_run_summary(out_dir, dataset, classification, signals,
                  rec_tau, rec_cons_j, override_used, result,
                  n_universe, med_frac_in_matrix)

write_comparison(out_dir, dataset, rec_tau, rec_cons_j, result, mini_sweep)

cat(sprintf("Done. Output: %s\n", out_dir))
