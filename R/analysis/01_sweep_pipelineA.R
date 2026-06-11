# R/analysis/01_sweep_pipelineA.R
#
# Threshold sensitivity sweep for Pipeline A using existing bootstrap matrices.
#
# IMPORTANT: Does NOT rerun bootstrap (expensive).
#   Reads existing go_bp_boot_matrix.csv and wp_boot_matrix.csv.
#   Loads gene sets once per dataset (org.Hs.eg.db + Pathway2Gene.csv) to
#   compute Jaccard for all CONS_JACCARD_MIN values including 0.10.
#
# Logic mirrors step06/step07 exactly:
#   stability[i]  = colMeans(boot_matrix)[i]
#   selected      = stability >= TAU_STABILITY
#   consensus     = term appears in any (GO,WP) pair with jaccard >= CONS_J & k >= CONS_K_MIN
#   Tier 1        = selected AND consensus   (keep_final in step06)
#   Tier 2        = selected AND NOT consensus
#   n_final       = step98 export logic: tier1>=5 -> tier1; tier2>0 -> tier2; else candidates
#
# CONS_K_MIN is fixed at 5 (baseline) throughout.
# CONS_JACCARD_MIN and TAU_STABILITY are swept.
#
# Output: results/analysis/sweep_pipelineA.csv
#
# Run from project root:
#   Rscript R/analysis/01_sweep_pipelineA.R

suppressPackageStartupMessages({
  library(readr)
  library(org.Hs.eg.db)   # load before anything that might mask select()
  library(AnnotationDbi)
  # dplyr is NOT loaded globally: it masks AnnotationDbi::select() and breaks
  # S4 dispatch even when calling AnnotationDbi::select() explicitly.
  # Use :: prefix if dplyr operations are needed.
})

# ---- Configuration -----------------------------------------------------------

ROOT <- getwd()  # must be run from project root

TAU_VALUES    <- c(0.50, 0.60, 0.70)   # TAU_STABILITY sweep
CONS_J_VALUES <- c(0.10, 0.15, 0.20)   # CONS_JACCARD_MIN sweep
CONS_K_MIN    <- 5L                     # fixed (baseline value; not swept)

# Baseline run directories (relative to ROOT).
RUNS <- list(
  list(dataset = "dataset_0", dir = "results/pipelineA/DATASE_0_BASELINE_2026-03-01_18-00-42"),
  list(dataset = "dataset_1", dir = "results/pipelineA/DATASE_1_BASELINE_2026-03-03_15-00-39"),
  list(dataset = "dataset_2", dir = "results/pipelineA/DATASET_2_BASELINE_2026-03-04_11-29-07"),
  list(dataset = "dataset_3", dir = "results/pipelineA/DATASET_3_BASELINE_2026-03-04_14-54-56")
)

OUT_DIR  <- file.path(ROOT, "results", "analysis")
OUT_FILE <- file.path(OUT_DIR, "sweep_pipelineA.csv")

# ---- Helper: load the bootstrap matrix ---------------------------------------
# Returns a named numeric vector: stability[pathway_id] = colMean.
# Uses readr to preserve "GO:xxxxx" colnames (base read.csv replaces : with .).
load_stability <- function(path) {
  if (!file.exists(path)) return(NULL)
  mat_raw <- read_csv(path, show_col_types = FALSE)
  # First column is "bootstrap" (row labels b1..bB); rest are pathway IDs
  ids  <- names(mat_raw)[-1]
  mat  <- as.matrix(mat_raw[, -1])
  stab <- colMeans(mat, na.rm = TRUE)
  names(stab) <- ids
  stab
}

# ---- Helper: read a plain gene-per-line text file ----------------------------
read_genes <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- readLines(path, warn = FALSE)
  unique(trimws(x[nzchar(trimws(x))]))
}

# ---- Helper: build GO gene sets from org.Hs.eg.db, filtered to universe -----
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

  # Base R only — avoids dplyr masking AnnotationDbi::select
  raw$ID <- as.character(raw[[id_col[1]]])
  raw <- raw[!is.na(raw$ID) & !is.na(raw$SYMBOL), ]
  raw <- raw[is.na(raw$ONTOLOGY) | raw$ONTOLOGY == "BP", ]
  raw <- raw[raw$SYMBOL %in% universe, ]
  raw <- raw[!duplicated(raw[, c("ID", "SYMBOL")]), ]

  # Build named list: ID -> character vector of gene symbols
  gene_list <- split(raw$SYMBOL, raw$ID)
  lapply(gene_list, unique)
}

# ---- Helper: build WP gene sets from Pathway2Gene.csv -----------------------
build_wp_sets <- function(wp_ids, universe, wp_map_path) {
  if (!file.exists(wp_map_path)) stop("Pathway2Gene.csv not found: ", wp_map_path)
  if (length(wp_ids) == 0) return(list())

  raw <- read_csv(wp_map_path, show_col_types = FALSE)

  # Detect id and symbol columns flexibly (same logic as step06)
  id_candidates  <- c("wpid", "ID", "id", "Pathway", "pathway", "pathway_id", "WPID", "WP_ID", "WP")
  sym_candidates <- c("SYMBOL", "symbol", "GeneSymbol", "gene_symbol", "Gene", "gene",
                      "HGNC_symbol", "HUGO")

  id_col  <- names(raw)[tolower(names(raw)) %in% tolower(id_candidates)]
  sym_col <- names(raw)[tolower(names(raw)) %in% tolower(sym_candidates)]

  if (!length(id_col) || !length(sym_col)) {
    stop("Pathway2Gene.csv missing expected ID/SYMBOL columns. Found: ",
         paste(names(raw), collapse = ", "))
  }

  # Base R only
  tbl <- data.frame(
    ID     = as.character(raw[[id_col[1]]]),
    SYMBOL = as.character(raw[[sym_col[1]]]),
    stringsAsFactors = FALSE
  )
  tbl <- tbl[tbl$ID %in% wp_ids & tbl$SYMBOL %in% universe, ]
  tbl <- tbl[!duplicated(tbl), ]

  # Build named list: ID -> character vector of gene symbols
  gene_list <- split(tbl$SYMBOL, tbl$ID)
  lapply(gene_list, unique)
}

# ---- Helper: compute ALL pairwise Jaccard (no threshold) --------------------
# Returns a data frame: A (GO id), B (WP id), jaccard, k (intersection size)
compute_all_jaccard <- function(go_sets, wp_sets) {
  if (!length(go_sets) || !length(wp_sets)) {
    return(data.frame(A = character(), B = character(),
                      jaccard = numeric(), k = integer(),
                      stringsAsFactors = FALSE))
  }
  a_ids <- names(go_sets)
  b_ids <- names(wp_sets)

  rows <- vector("list", length(a_ids) * length(b_ids))
  idx  <- 0L
  for (a in a_ids) {
    for (b in b_ids) {
      inter <- length(intersect(go_sets[[a]], wp_sets[[b]]))
      uni   <- length(union(go_sets[[a]], wp_sets[[b]]))
      if (uni == 0L) next
      idx <- idx + 1L
      rows[[idx]] <- data.frame(A = a, B = b,
                                jaccard = inter / uni,
                                k = inter,
                                stringsAsFactors = FALSE)
    }
  }
  if (idx == 0L) return(data.frame(A = character(), B = character(),
                                    jaccard = numeric(), k = integer(),
                                    stringsAsFactors = FALSE))
  do.call(rbind, rows[seq_len(idx)])
}

# ---- Helper: apply sweep for one (tau, cons_j) combination ------------------
apply_thresholds <- function(stability_go, stability_wp,
                              all_pairs, n_candidates_fixed,
                              tau, cons_j, cons_k) {
  # 1. Stability -> selected
  go_selected <- names(stability_go)[stability_go >= tau]
  wp_selected <- names(stability_wp)[stability_wp >= tau]

  # 2. Filter pairs by Jaccard and k thresholds (same logic as step06)
  pairs_filt <- all_pairs[all_pairs$jaccard >= cons_j & all_pairs$k >= cons_k, ,
                           drop = FALSE]

  # 3. Consensus: any term appearing in a filtered pair
  go_cons_ids <- unique(pairs_filt$A)
  wp_cons_ids <- unique(pairs_filt$B)

  # 4. Tier counts
  n_go_stable   <- length(go_selected)
  n_wp_stable   <- length(wp_selected)
  n_tier1_go    <- sum(go_selected %in% go_cons_ids)
  n_tier1_wp    <- sum(wp_selected %in% wp_cons_ids)
  n_tier1       <- n_tier1_go + n_tier1_wp
  n_tier2_go    <- n_go_stable - n_tier1_go
  n_tier2_wp    <- n_wp_stable - n_tier1_wp
  n_tier2       <- n_tier2_go + n_tier2_wp
  n_pairs       <- nrow(pairs_filt)

  # 5. n_final via step98 logic
  n_final <- if (n_tier1 >= 5L) {
    n_tier1
  } else if (n_tier2 > 0L) {
    n_tier2
  } else {
    n_candidates_fixed  # fixed fallback: all reps from step07
  }

  # 6. final_source label (mirrors step98 console output)
  final_source <- if (n_tier1 >= 5L) {
    "tier1"
  } else if (n_tier2 > 0L) {
    "tier2"
  } else {
    "candidates"
  }

  list(
    n_go_stable      = n_go_stable,
    n_wp_stable      = n_wp_stable,
    n_tier1          = n_tier1,
    n_tier2          = n_tier2,
    n_consensus_pairs = n_pairs,
    n_final          = n_final,
    final_source     = final_source
  )
}

# ---- Main: loop over datasets ------------------------------------------------

cat("=== 01_sweep_pipelineA.R ===\n")
cat(sprintf("TAU_STABILITY values:    %s\n", paste(TAU_VALUES, collapse = ", ")))
cat(sprintf("CONS_JACCARD_MIN values: %s\n", paste(CONS_J_VALUES, collapse = ", ")))
cat(sprintf("CONS_K_MIN (fixed):      %d\n\n", CONS_K_MIN))

all_rows   <- list()
row_idx    <- 0L
skip_notes <- character(0)

for (entry in RUNS) {
  dataset <- entry$dataset
  dir     <- file.path(ROOT, entry$dir)

  cat(sprintf("--- %s ---\n", dataset))

  if (!dir.exists(dir)) {
    msg <- sprintf("SKIP: run directory not found: %s", entry$dir)
    cat(msg, "\n\n")
    skip_notes <- c(skip_notes, msg)
    next
  }

  # ---- Load bootstrap matrices -----------------------------------------------
  go_stab <- load_stability(file.path(dir, "06_bootstrap_consensus", "go_bp_boot_matrix.csv"))
  wp_stab <- load_stability(file.path(dir, "06_bootstrap_consensus", "wp_boot_matrix.csv"))

  if (is.null(go_stab) || is.null(wp_stab)) {
    msg <- sprintf("SKIP %s: boot matrix missing", dataset)
    cat(msg, "\n\n")
    skip_notes <- c(skip_notes, msg)
    next
  }
  cat(sprintf("  Boot matrix: %d GO reps, %d WP reps\n",
              length(go_stab), length(wp_stab)))

  # ---- Load universe ----------------------------------------------------------
  universe <- read_genes(file.path(dir, "01_universe", "universe_genes.txt"))
  if (length(universe) == 0) {
    msg <- sprintf("SKIP %s: universe_genes.txt missing or empty", dataset)
    cat(msg, "\n\n")
    skip_notes <- c(skip_notes, msg)
    next
  }

  # ---- Load representative IDs (needed to restrict gene set building) --------
  go_reps_path <- file.path(dir, "05_overlap", "go_bp", "representatives.csv")
  wp_reps_path <- file.path(dir, "05_overlap", "wp",    "representatives.csv")

  if (!file.exists(go_reps_path) || !file.exists(wp_reps_path)) {
    msg <- sprintf("SKIP %s: step05 representatives.csv missing", dataset)
    cat(msg, "\n\n")
    skip_notes <- c(skip_notes, msg)
    next
  }
  go_rep_ids <- as.character(read_csv(go_reps_path, show_col_types = FALSE)$ID)
  wp_rep_ids <- as.character(read_csv(wp_reps_path, show_col_types = FALSE)$ID)

  # Restrict stability vectors to the representatives (guard against stale matrices)
  go_stab <- go_stab[names(go_stab) %in% go_rep_ids]
  wp_stab <- wp_stab[names(wp_stab) %in% wp_rep_ids]

  # ---- Load Pathway2Gene path from config_used.yml ---------------------------
  config_path <- file.path(dir, "config_used.yml")
  wp_map_path <- NULL
  if (file.exists(config_path)) {
    cfg <- tryCatch(yaml::read_yaml(config_path), error = function(e) NULL)
    wp_map_path <- cfg$dataset$pathway2gene_file
    if (!is.null(wp_map_path) && !startsWith(wp_map_path, "/")) {
      wp_map_path <- file.path(ROOT, wp_map_path)
    }
  }
  if (is.null(wp_map_path) || !file.exists(wp_map_path)) {
    # Fallback to known shared path
    wp_map_path <- file.path(ROOT, "data", "dataset_0", "processed", "Pathway2Gene.csv")
  }

  # ---- Build gene sets (once per dataset, not per sweep iteration) -----------
  cat("  Building GO gene sets from org.Hs.eg.db ...\n")
  go_sets <- tryCatch(
    build_go_sets(names(go_stab), universe),
    error = function(e) {
      cat("  WARNING: GO set build failed:", conditionMessage(e), "\n")
      list()
    }
  )

  cat(sprintf("  GO gene sets built: %d / %d reps have sets\n",
              length(go_sets), length(go_stab)))

  cat("  Building WP gene sets from Pathway2Gene.csv ...\n")
  wp_sets <- tryCatch(
    build_wp_sets(names(wp_stab), universe, wp_map_path),
    error = function(e) {
      cat("  WARNING: WP set build failed:", conditionMessage(e), "\n")
      list()
    }
  )
  cat(sprintf("  WP gene sets built: %d / %d reps have sets\n",
              length(wp_sets), length(wp_stab)))

  # Drop reps without gene sets (mirrors step06 behaviour)
  go_stab <- go_stab[names(go_stab) %in% names(go_sets)]
  wp_stab <- wp_stab[names(wp_stab) %in% names(wp_sets)]

  # ---- Compute all pairwise Jaccard (once; no threshold filter yet) ----------
  cat("  Computing all pairwise Jaccard (GO x WP) ...\n")
  all_pairs <- compute_all_jaccard(go_sets, wp_sets)
  cat(sprintf("  Total GO×WP pairs computed: %d\n", nrow(all_pairs)))

  # ---- Fixed n_candidates (from existing step07 output; does not vary with TAU)
  n_cand_path <- file.path(dir, "07_tiers", "tiered_all_candidates.csv")
  n_candidates_fixed <- if (file.exists(n_cand_path)) {
    max(0L, length(readLines(n_cand_path, warn = FALSE)) - 1L)
  } else {
    NA_integer_
  }

  # ---- Sweep -----------------------------------------------------------------
  cat("  Sweeping thresholds ...\n")
  for (tau in TAU_VALUES) {
    for (cons_j in CONS_J_VALUES) {
      res <- apply_thresholds(
        go_stab, wp_stab, all_pairs, n_candidates_fixed,
        tau = tau, cons_j = cons_j, cons_k = CONS_K_MIN
      )

      row_idx <- row_idx + 1L
      all_rows[[row_idx]] <- data.frame(
        dataset           = dataset,
        run_dir           = entry$dir,
        tau_stability     = tau,
        cons_jaccard_min  = cons_j,
        cons_k_min        = CONS_K_MIN,
        n_go_reps         = length(go_stab),
        n_wp_reps         = length(wp_stab),
        n_go_stable       = res$n_go_stable,
        n_wp_stable       = res$n_wp_stable,
        n_tier1           = res$n_tier1,
        n_tier2           = res$n_tier2,
        n_consensus_pairs = res$n_consensus_pairs,
        n_final           = res$n_final,
        final_source      = res$final_source,
        stringsAsFactors  = FALSE
      )
    }
  }
  cat("\n")
}

if (row_idx == 0L) stop("No rows produced — check run directories and inputs.")

sweep_df <- do.call(rbind, all_rows)
rownames(sweep_df) <- NULL

# ---- Write output ------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(sweep_df, OUT_FILE, row.names = FALSE)
cat(sprintf("Wrote: %s\n\n", OUT_FILE))

# ---- Print summary -----------------------------------------------------------

cat("=== SUMMARY ===\n\n")

# Per-dataset: show Tier1 count across TAU (at baseline CONS_J=0.15)
for (dataset in unique(sweep_df$dataset)) {
  sub <- sweep_df[sweep_df$dataset == dataset & sweep_df$cons_jaccard_min == 0.15, ]
  cat(sprintf("%-12s  (CONS_J=0.15 fixed)\n", dataset))
  cat(sprintf("  %-6s  %10s  %10s  %10s  %14s  %12s\n",
              "TAU", "n_go_stab", "n_wp_stab", "n_tier1", "n_consensus_pr", "final_source"))
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    marker <- if (r$n_tier1 == 0) " [no Tier1]" else ""
    cat(sprintf("  %-6.2f  %10d  %10d  %10d  %14d  %12s%s\n",
                r$tau_stability,
                r$n_go_stable, r$n_wp_stable,
                r$n_tier1, r$n_consensus_pairs,
                r$final_source,
                marker))
  }
  cat("\n")
}

# Show how CONS_J affects Tier1 (at baseline TAU=0.70)
cat("Effect of CONS_JACCARD_MIN on Tier1 (TAU=0.70 fixed):\n")
cat(sprintf("  %-12s  %-8s  %14s  %8s  %12s\n",
            "dataset", "CONS_J", "n_consensus_pr", "n_tier1", "final_source"))
for (dataset in unique(sweep_df$dataset)) {
  sub <- sweep_df[sweep_df$dataset == dataset & sweep_df$tau_stability == 0.70, ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    cat(sprintf("  %-12s  %-8.2f  %14d  %8d  %12s\n",
                r$dataset, r$cons_jaccard_min,
                r$n_consensus_pairs, r$n_tier1, r$final_source))
  }
}
cat("\n")

# Flag datasets where Tier1 is always 0
dataset_names    <- unique(sweep_df$dataset)
always_zero_mask <- vapply(dataset_names, function(d) {
  all(sweep_df$n_tier1[sweep_df$dataset == d] == 0L)
}, logical(1))
always_zero <- dataset_names[always_zero_mask]
if (length(always_zero) > 0) {
  cat("DATASETS WHERE TIER1 IS ALWAYS EMPTY (across all threshold combinations):\n")
  for (d in always_zero) cat(" *", d, "\n")
  cat("\n")
}

# Flag settings where Tier1 first appears per dataset
cat("First threshold setting where Tier1 > 0 (per dataset):\n")
for (dataset in unique(sweep_df$dataset)) {
  sub <- sweep_df[sweep_df$dataset == dataset & sweep_df$n_tier1 > 0, ]
  if (nrow(sub) == 0L) {
    cat(sprintf("  %-12s: never — Tier1 remains empty at all tested thresholds\n", dataset))
  } else {
    # Pick the strictest TAU + loosest CONS_J combo that first shows Tier1
    # Sort by TAU desc (strictest first), CONS_J asc (loosest first for same TAU)
    sub_ord <- sub[order(-sub$tau_stability, sub$cons_jaccard_min), ]
    r <- sub_ord[1, ]
    cat(sprintf("  %-12s: TAU=%.2f, CONS_J=%.2f  -> Tier1=%d\n",
                dataset, r$tau_stability, r$cons_jaccard_min, r$n_tier1))
  }
}

if (length(skip_notes) > 0) {
  cat("\nSkipped entries:\n")
  for (s in skip_notes) cat(" -", s, "\n")
}

cat("\nDone.\n")
