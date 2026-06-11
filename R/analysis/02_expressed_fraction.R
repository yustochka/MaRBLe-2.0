# R/analysis/02_expressed_fraction.R
# QC: for every FINAL pathway (all 4 datasets × both pipelines), report what
# fraction of its member genes are present in the expression matrix and what
# fraction are low-expressed.
#
# Does NOT re-run any pipeline step. Reads from existing baseline run outputs.
#
# Output: results/analysis/expressed_fraction.csv
#
# Run from project root:
#   Rscript R/analysis/02_expressed_fraction.R

suppressPackageStartupMessages({
  library(readr)
  library(org.Hs.eg.db)   # load before anything that might mask select()
  library(AnnotationDbi)
})

ROOT <- getwd()  # must be run from project root

LOW_EXPR_QUANTILE <- 0.10   # bottom-10% of gene means → "low expression"

# ---- Dataset registry -------------------------------------------------------
# Update this list if new baseline runs are added or paths change.
DATASETS <- list(
  list(
    name          = "dataset_0",
    expr_file     = "data/dataset_0/processed/processedMatrix.csv",
    wp_file       = "data/dataset_0/processed/Pathway2Gene.csv",
    gene_id_type  = "SYMBOL",   # expression row IDs are already gene symbols
    pipelineA_dir = "results/pipelineA/DATASE_0_BASELINE_2026-03-01_18-00-42",
    pipelineB_dir = "results/pipelineB/DATASE_0_BASELINE_2026-03-01_18-04-29"
  ),
  list(
    name          = "dataset_1",
    expr_file     = "data/dataset_1/normCounts_GSE199939.csv",
    wp_file       = "data/dataset_0/processed/Pathway2Gene.csv",
    gene_id_type  = "ENTREZ",
    pipelineA_dir = "results/pipelineA/DATASE_1_BASELINE_2026-03-03_15-00-39",
    pipelineB_dir = "results/pipelineB/DATASE_1_BASELINE_2026-03-04_09-54-24"
  ),
  list(
    name          = "dataset_2",
    expr_file     = "data/dataset_2/normCounts_GSE243836.csv",
    wp_file       = "data/dataset_0/processed/Pathway2Gene.csv",
    gene_id_type  = "ENTREZ",
    pipelineA_dir = "results/pipelineA/DATASET_2_BASELINE_2026-03-04_11-29-07",
    pipelineB_dir = "results/pipelineB/DATASET_2_BASELINE_2026-03-04_14-43-36"
  ),
  list(
    name          = "dataset_3",
    expr_file     = "data/dataset_3/normCounts_GSE247345.csv",
    wp_file       = "data/dataset_0/processed/Pathway2Gene.csv",
    gene_id_type  = "ENTREZ",
    pipelineA_dir = "results/pipelineA/DATASET_3_BASELINE_2026-03-04_14-54-56",
    pipelineB_dir = "results/pipelineB/DATASET_3_BASELINE_2026-03-04_12-32-56"
  )
)

# ---- Helper: load expression matrix, return named gene-mean vector ----------
# For ENTREZ datasets, converts row IDs to SYMBOL via org.Hs.eg.db (base R,
# no dplyr, to avoid masking AnnotationDbi::select).
# Returns a numeric named vector: names = gene SYMBOL, values = row means.
load_gene_means <- function(expr_file, gene_id_type) {
  raw <- read_csv(file.path(ROOT, expr_file), show_col_types = FALSE)

  # First column holds gene IDs (unnamed in the CSV → read_csv calls it "...1")
  id_col  <- names(raw)[1]
  ids     <- as.character(raw[[id_col]])
  expr_mat <- as.matrix(raw[, -1, drop = FALSE])
  mode(expr_mat) <- "numeric"

  if (gene_id_type == "SYMBOL") {
    symbols <- ids
  } else if (gene_id_type == "ENTREZ") {
    # Map ENTREZ → first SYMBOL returned by org.Hs.eg.db (base R only)
    map_raw <- suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys    = unique(ids),
        keytype = "ENTREZID",
        columns = "SYMBOL"
      )
    )
    # Build lookup: ENTREZ -> first SYMBOL (drop unmapped)
    map_raw <- map_raw[!is.na(map_raw$SYMBOL), ]
    map_raw <- map_raw[!duplicated(map_raw$ENTREZID), ]
    lut <- setNames(map_raw$SYMBOL, as.character(map_raw$ENTREZID))

    symbols <- lut[ids]   # NA for unmapped IDs
  } else {
    stop("Unsupported gene_id_type: ", gene_id_type)
  }

  # Keep only rows with a valid symbol
  keep <- !is.na(symbols) & symbols != ""
  expr_mat <- expr_mat[keep, , drop = FALSE]
  symbols  <- symbols[keep]

  # For duplicate symbols, keep the row with the highest mean expression
  row_means_all <- rowMeans(expr_mat, na.rm = TRUE)
  ord  <- order(row_means_all, decreasing = TRUE)
  expr_mat <- expr_mat[ord, , drop = FALSE]
  symbols  <- symbols[ord]
  keep2 <- !duplicated(symbols)
  expr_mat <- expr_mat[keep2, , drop = FALSE]
  symbols  <- symbols[keep2]

  gene_means <- rowMeans(expr_mat, na.rm = TRUE)
  names(gene_means) <- symbols
  gene_means
}

# ---- Helper: build WP gene sets (base R, no dplyr) --------------------------
# Returns a named list of character vectors (ALL symbols in pathway, not filtered).
build_wp_sets <- function(wp_ids, wp_file) {
  if (length(wp_ids) == 0) return(list())
  raw <- read_csv(file.path(ROOT, wp_file), show_col_types = FALSE)

  id_cands  <- c("wpid","id","pathway_id","pathway","wp_id","wp","WPID","WP_ID")
  sym_cands <- c("SYMBOL","symbol","gene","gene_symbol","genesymbol","hgnc_symbol","hugo","GeneSymbol")

  id_col  <- names(raw)[tolower(names(raw)) %in% tolower(id_cands)][1]
  sym_col <- names(raw)[tolower(names(raw)) %in% tolower(sym_cands)][1]

  if (is.na(id_col) || is.na(sym_col))
    stop("Pathway2Gene.csv missing ID/SYMBOL columns. Found: ", paste(names(raw), collapse = ", "))

  tbl <- data.frame(
    ID     = as.character(raw[[id_col]]),
    SYMBOL = as.character(raw[[sym_col]]),
    stringsAsFactors = FALSE
  )
  tbl <- tbl[tbl$ID %in% wp_ids, ]   # NOT filtered to universe here
  tbl <- tbl[!duplicated(tbl), ]

  gene_list <- split(tbl$SYMBOL, tbl$ID)
  lapply(gene_list, unique)
}

# ---- Helper: build GO-BP gene sets (base R, no dplyr) -----------------------
# Returns a named list of character vectors (ALL symbols in term, not filtered).
build_go_sets <- function(go_ids) {
  if (length(go_ids) == 0) return(list())

  raw <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = unique(go_ids),
      keytype = "GO",
      columns = c("SYMBOL", "ONTOLOGY")
    )
  )

  id_col <- intersect(c("GO", "GOID", "GOALL"), names(raw))[1]
  if (is.na(id_col))
    stop("org.Hs.eg.db returned no GO ID column. Got: ", paste(names(raw), collapse = ", "))
  if (!"ONTOLOGY" %in% names(raw)) raw$ONTOLOGY <- NA_character_

  raw$ID <- as.character(raw[[id_col]])
  raw    <- raw[!is.na(raw$ID) & !is.na(raw$SYMBOL), ]
  raw    <- raw[is.na(raw$ONTOLOGY) | raw$ONTOLOGY == "BP", ]
  raw    <- raw[!duplicated(raw[, c("ID", "SYMBOL")]), ]
  # NOT filtered to universe — keep all annotated genes

  gene_list <- split(raw$SYMBOL, raw$ID)
  lapply(gene_list, unique)
}

# ---- Helper: compute QC metrics for one pathway term -----------------------
qc_term <- function(term_id, members_in_universe, gene_means, low_cut) {
  genes   <- members_in_universe
  n_total <- length(genes)

  in_mat  <- genes[genes %in% names(gene_means)]
  n_in    <- length(in_mat)

  if (n_in == 0) {
    return(list(
      n_members       = n_total,
      n_in_matrix     = 0L,
      frac_in_matrix  = 0,
      frac_low_expr   = NA_real_,
      mean_expr_median = NA_real_
    ))
  }

  vals <- gene_means[in_mat]
  list(
    n_members        = n_total,
    n_in_matrix      = n_in,
    frac_in_matrix   = n_in / max(1L, n_total),
    frac_low_expr    = mean(vals <= low_cut),
    mean_expr_median = median(vals, na.rm = TRUE)
  )
}

# ---- Helper: process one FINAL.csv table ------------------------------------
# Returns a data.frame with one row per term.
process_final <- function(final_csv, dataset_name, pipeline, gene_means, wp_sets, go_sets) {
  if (!file.exists(final_csv)) {
    message("  [skip] not found: ", final_csv)
    return(data.frame())
  }

  fin <- read_csv(final_csv, show_col_types = FALSE)
  if (nrow(fin) == 0) return(data.frame())

  # Standardise the ID / Description / collection / direction columns
  fin$ID          <- as.character(fin$ID)
  fin$Description <- as.character(fin[["Description"]])
  fin$collection  <- as.character(fin$collection)
  fin$direction   <- if ("direction" %in% names(fin)) as.character(fin$direction) else NA_character_

  low_cut  <- quantile(gene_means, probs = LOW_EXPR_QUANTILE, na.rm = TRUE)

  rows <- lapply(seq_len(nrow(fin)), function(i) {
    id  <- fin$ID[i]
    col <- fin$collection[i]

    members <- if (col == "WP") {
      if (!is.null(wp_sets[[id]])) wp_sets[[id]] else character(0)
    } else {
      if (!is.null(go_sets[[id]])) go_sets[[id]] else character(0)
    }

    qc <- qc_term(id, members, gene_means, low_cut)

    data.frame(
      dataset          = dataset_name,
      pipeline         = pipeline,
      ID               = id,
      Description      = fin$Description[i],
      collection       = col,
      direction        = fin$direction[i],
      n_members        = qc$n_members,
      n_in_matrix      = qc$n_in_matrix,
      frac_in_matrix   = round(qc$frac_in_matrix,  4),
      frac_low_expr    = round(qc$frac_low_expr,    4),
      mean_expr_median = round(qc$mean_expr_median, 4),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# ---- Main loop --------------------------------------------------------------
cat("=== 02_expressed_fraction.R ===\n")
cat("Low-expression cutoff: bottom", LOW_EXPR_QUANTILE * 100, "% of gene means\n\n")

all_rows <- list()

for (ds in DATASETS) {
  cat("--- Dataset:", ds$name, "---\n")

  # 1) Load expression matrix → gene means
  cat("  Loading expression matrix:", ds$expr_file, "\n")
  gene_means <- tryCatch(
    load_gene_means(ds$expr_file, ds$gene_id_type),
    error = function(e) { message("  ERROR loading expr: ", e$message); NULL }
  )
  if (is.null(gene_means)) next
  cat("  Genes in matrix (post-ID conversion):", length(gene_means), "\n")

  universe <- names(gene_means)

  # 2) Collect all pathway IDs needed from both pipelines
  a_final <- file.path(ROOT, ds$pipelineA_dir, "FINAL", "FINAL.csv")
  b_final <- file.path(ROOT, ds$pipelineB_dir, "FINAL", "FINAL.csv")

  all_ids <- character(0)
  for (fp in c(a_final, b_final)) {
    if (file.exists(fp)) {
      tmp <- read_csv(fp, show_col_types = FALSE)
      all_ids <- unique(c(all_ids, as.character(tmp$ID)))
    }
  }
  go_ids <- grep("^GO:", all_ids, value = TRUE)
  wp_ids <- grep("^WP",  all_ids, value = TRUE)

  # 3) Build gene sets filtered to this dataset's universe
  cat("  Building GO sets (", length(go_ids), " terms)...\n", sep = "")
  go_sets <- tryCatch(
    build_go_sets(go_ids),
    error = function(e) { message("  ERROR building GO sets: ", e$message); list() }
  )

  cat("  Building WP sets (", length(wp_ids), " terms)...\n", sep = "")
  wp_sets <- tryCatch(
    build_wp_sets(wp_ids, ds$wp_file),
    error = function(e) { message("  ERROR building WP sets: ", e$message); list() }
  )

  # 4) Process each pipeline's FINAL.csv
  for (pipe_info in list(
    list(label = "pipelineA", csv = a_final),
    list(label = "pipelineB", csv = b_final)
  )) {
    cat("  Processing", pipe_info$label, "...\n")
    rows <- process_final(pipe_info$csv, ds$name, pipe_info$label,
                          gene_means, wp_sets, go_sets)
    if (nrow(rows) > 0) {
      all_rows[[length(all_rows) + 1]] <- rows
      cat("    Rows added:", nrow(rows), "\n")
    }
  }

  cat("\n")
}

# ---- Write output -----------------------------------------------------------
if (length(all_rows) == 0) {
  stop("No rows collected — check that FINAL.csv files exist.")
}

result <- do.call(rbind, all_rows)
rownames(result) <- NULL

out_dir <- file.path(ROOT, "results", "analysis")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(out_dir, "expressed_fraction.csv")
write_csv(result, out_file)

# ---- Print summary ----------------------------------------------------------
cat("=== Summary ===\n")
cat(sprintf("%-12s  %-10s  %4s  %6s  %6s  %6s  %7s\n",
            "dataset", "pipeline", "N", "medFrac", "lo10%", "medMed", "problem"))
cat(strrep("-", 65), "\n", sep = "")

for (ds_name in sapply(DATASETS, `[[`, "name")) {
  for (pipe in c("pipelineA", "pipelineB")) {
    sub <- result[result$dataset == ds_name & result$pipeline == pipe, ]
    if (nrow(sub) == 0) next

    med_frac   <- median(sub$frac_in_matrix,   na.rm = TRUE)
    med_low    <- median(sub$frac_low_expr,    na.rm = TRUE)
    med_expr   <- median(sub$mean_expr_median, na.rm = TRUE)
    n_prob     <- sum(sub$frac_in_matrix < 0.5, na.rm = TRUE)

    cat(sprintf("%-12s  %-10s  %4d  %6.3f  %6.3f  %6.1f  %7d\n",
                ds_name, pipe, nrow(sub),
                med_frac, med_low, med_expr, n_prob))
  }
}

cat(strrep("-", 65), "\n", sep = "")
cat("problem = pathways with < 50% of member genes in the expression matrix\n")
cat("\nOutput written to:", out_file, "\n")
