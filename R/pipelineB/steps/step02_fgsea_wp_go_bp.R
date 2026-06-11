# R/pipelineB/steps/step02_fgsea_wp_go_bp.R
# Pipeline B — Step 02: Sensitivity analysis (ranked, cutoff-free) via fgsea
# Runs fgsea on:
#   (1) WikiPathways (from Pathway2Gene.csv)
#   (2) GO:BP        (from org.Hs.eg.db SYMBOL -> GO mapping)
#
# Reads:
#   cfg$dataset$expression_file
#   cfg$dataset$metadata_file
#   paths$step01/universe_genes.txt   (written by step01)
#   cfg$dataset$pathway2gene_file     (e.g. data/dataset_0/processed/Pathway2Gene.csv)
#
# Writes into paths$step02  (results/pipelineB/<run_id>/02_fgsea/):
#   fgsea_wikipathways.csv
#   fgsea_go_bp.csv
#   run_note.txt
#
# Returns (invisibly):
#   fgsea_wp_file   path to fgsea_wikipathways.csv
#   fgsea_go_file   path to fgsea_go_bp.csv

suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(org.Hs.eg.db)
})

if (!exists("convert_expr_to_symbol")) source("R/utils/gene_id_convert.R")

step02_fgsea_wp_go_bp <- function(cfg, paths, ctx) {

  # ---- Parameters (from cfg with script-matching defaults) ---------------
  EXPR_FILE     <- cfg$dataset$expression_file
  META_FILE     <- cfg$dataset$metadata_file
  UNIVERSE_FILE <- file.path(paths$step01, "universe_genes.txt")
  WP_MAP_FILE   <- cfg$dataset$pathway2gene_file
  TARGET_TISSUE <- if (!is.null(cfg$pipelineB$target_tissue)) cfg$pipelineB$target_tissue else "Temporal cortex"
  MIN_SET_SIZE  <- if (!is.null(cfg$pipelineB$min_set_size))  as.integer(cfg$pipelineB$min_set_size)  else 10L
  MAX_SET_SIZE  <- if (!is.null(cfg$pipelineB$max_set_size))  as.integer(cfg$pipelineB$max_set_size)  else 500L
  INCLUDE_AGE   <- if (!is.null(cfg$pipelineB$include_age))   isTRUE(cfg$pipelineB$include_age)       else FALSE
  SET_SEED      <- if (!is.null(cfg$pipelineB$gsea_seed))     as.integer(cfg$pipelineB$gsea_seed)     else 7L

  stopifnot(file.exists(EXPR_FILE), file.exists(META_FILE),
            file.exists(UNIVERSE_FILE), file.exists(WP_MAP_FILE))

  # ---- Helpers -----------------------------------------------------------
  read_gene_vec <- function(path) {
    x <- readr::read_lines(path)
    x <- stringr::str_trim(x)
    unique(x[x != ""])
  }

  getcol <- function(df, candidates) {
    nm <- names(df); nml <- tolower(nm)
    hit <- nm[match(tolower(candidates), nml, nomatch = 0)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }

  collapse_duplicate_genes <- function(expr_mat) {
    if (!anyDuplicated(rownames(expr_mat))) return(expr_mat)
    v        <- apply(expr_mat, 1, stats::var, na.rm = TRUE)
    keep_idx <- tapply(seq_along(v), rownames(expr_mat), function(ix) ix[which.max(v[ix])])
    keep_idx <- unname(unlist(keep_idx))
    expr_mat[keep_idx, , drop = FALSE]
  }

  make_pathways_list <- function(set_df, min_size, max_size) {
    tmp <- set_df |>
      dplyr::filter(!is.na(set_id), !is.na(gene), gene != "") |>
      dplyr::distinct(set_id, gene) |>
      dplyr::group_by(set_id) |>
      dplyr::summarise(genes = list(unique(gene)), .groups = "drop") |>
      dplyr::mutate(size = lengths(genes)) |>
      dplyr::filter(size >= min_size, size <= max_size)
    stats::setNames(tmp$genes, tmp$set_id)
  }

  leading_edge_to_string <- function(x) {
    if (is.null(x) || length(x) == 0) return("")
    paste(x, collapse = ";")
  }

  # ---- 1) Read expression + metadata -------------------------------------
  expr_raw   <- readr::read_csv(EXPR_FILE, show_col_types = FALSE)
  meta       <- readr::read_csv(META_FILE, show_col_types = FALSE)

  gene_col   <- names(expr_raw)[1]
  expr_genes <- as.character(expr_raw[[gene_col]])

  expr_mat <- expr_raw |>
    dplyr::select(-dplyr::all_of(gene_col)) |>
    as.matrix()

  mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- expr_genes
  expr_samples <- colnames(expr_mat)

  # ---- Gene ID conversion (ENTREZ/ENSEMBL → SYMBOL if configured) ---------
  expr_mat <- convert_expr_to_symbol(expr_mat, cfg)

  # ---- Column mapping (from config; defaults match dataset_0) ---------------
  SID_COL <- if (!is.null(cfg$dataset$sample_id_col) && nzchar(cfg$dataset$sample_id_col))
    cfg$dataset$sample_id_col else "SampleID"
  GRP_COL <- if (!is.null(cfg$dataset$group_col) && nzchar(cfg$dataset$group_col))
    cfg$dataset$group_col else "Genotype"

  if (!SID_COL %in% names(meta))
    stop("sample_id_col '", SID_COL, "' not found in metadata.\n",
         "  Available columns: ", paste(names(meta), collapse = ", "), "\n",
         "  Set cfg$dataset$sample_id_col in config.")
  if (!GRP_COL %in% names(meta))
    stop("group_col '", GRP_COL, "' not found in metadata.\n",
         "  Available columns: ", paste(names(meta), collapse = ", "), "\n",
         "  Set cfg$dataset$group_col in config.")

  # Trim whitespace on sample IDs in both metadata and expression matrix
  expr_samples_trimmed <- stringr::str_trim(expr_samples)
  colnames(expr_mat)   <- expr_samples_trimmed

  meta <- meta |> dplyr::mutate(dplyr::across(dplyr::all_of(SID_COL), stringr::str_trim))

  # Duplicate sample ID guard
  meta_ids <- as.character(meta[[SID_COL]])
  dupes    <- meta_ids[duplicated(meta_ids)]
  if (length(dupes) > 0)
    stop("Duplicate sample IDs in metadata: ", paste(unique(dupes), collapse = ", "))

  # Rename to stable internal names for model.matrix
  meta_use <- meta |>
    dplyr::rename(.sid = dplyr::all_of(SID_COL), .grp = dplyr::all_of(GRP_COL)) |>
    dplyr::filter(.sid %in% expr_samples_trimmed)

  if ("Tissue" %in% names(meta_use)) {
    meta_use <- meta_use |> dplyr::filter(Tissue == TARGET_TISSUE)
  }

  common_samples <- intersect(expr_samples_trimmed, meta_use$.sid)

  # Warn about mismatched samples
  extra_expr <- setdiff(expr_samples_trimmed, meta_use$.sid)
  if (length(extra_expr) > 0)
    warning(length(extra_expr), " expression column(s) not in metadata; dropping: ",
            paste(head(extra_expr, 5), collapse = ", "))

  extra_meta <- setdiff(meta_use$.sid, expr_samples_trimmed)
  if (length(extra_meta) > 0)
    message(length(extra_meta), " metadata sample(s) not in expression matrix; dropping: ",
            paste(head(extra_meta, 5), collapse = ", "))

  if (length(common_samples) < 4)
    stop("Too few samples after matching metadata to expression (", length(common_samples), "). ",
         "Check sample_id_col ('", SID_COL, "') and TARGET_TISSUE.")

  expr_mat <- expr_mat[, common_samples, drop = FALSE]
  meta_use <- meta_use |>
    dplyr::filter(.sid %in% common_samples) |>
    dplyr::mutate(.sid = factor(.sid, levels = common_samples)) |>
    dplyr::arrange(.sid)

  stopifnot(identical(as.character(meta_use$.sid), colnames(expr_mat)))

  meta_use <- meta_use |> dplyr::mutate(.grp = as.factor(.grp))

  if (length(levels(meta_use$.grp)) < 2)
    stop("Only one group level found in '", GRP_COL, "': ",
         paste(levels(meta_use$.grp), collapse = ", "),
         ". Check cfg$dataset$group_col.")

  cat(sprintf("Group '%s' levels: %s\n", GRP_COL,
              paste(levels(meta_use$.grp), collapse = " / ")))

  # ---- 2) Restrict to universe + clean duplicates ------------------------
  universe   <- read_gene_vec(UNIVERSE_FILE)
  keep_genes <- intersect(rownames(expr_mat), universe)
  expr_mat_u <- expr_mat[keep_genes, , drop = FALSE]

  if (nrow(expr_mat_u) < 1000) {
    warning("Very few genes left after universe restriction (", nrow(expr_mat_u), "). ",
            "Check that universe IDs match rownames (SYMBOL vs ENSEMBL).")
  }

  expr_mat_u <- collapse_duplicate_genes(expr_mat_u)

  # ---- 3) limma model -> ranked gene list (moderated t) ------------------
  if (INCLUDE_AGE) {
    AGE_COL <- if (!is.null(cfg$dataset$age_col) && nzchar(cfg$dataset$age_col))
      cfg$dataset$age_col else "Age"
    if (!AGE_COL %in% names(meta_use))
      stop("INCLUDE_AGE=TRUE but age_col '", AGE_COL, "' not found in metadata.")
    meta_use <- meta_use |> dplyr::rename(.age = dplyr::all_of(AGE_COL))
    design <- stats::model.matrix(~ .grp + .age, data = meta_use)
  } else {
    design <- stats::model.matrix(~ .grp, data = meta_use)
  }

  coef_names <- colnames(design)
  geno_coef  <- coef_names[stringr::str_detect(coef_names, "^\\.grp")]
  if (length(geno_coef) == 0) {
    stop("Could not find group coefficient in design matrix. Coefs: ",
         paste(coef_names, collapse = ", "))
  }
  geno_coef <- geno_coef[1]

  fit <- limma::lmFit(expr_mat_u, design)
  fit <- limma::eBayes(fit)

  tt <- limma::topTable(fit, coef = geno_coef, number = Inf, sort.by = "none")
  stopifnot(nrow(tt) == nrow(expr_mat_u))

  ranks <- tt$t
  names(ranks) <- rownames(expr_mat_u)
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)

  # ---- 4) fgsea on WikiPathways ------------------------------------------
  wp_map   <- readr::read_csv(WP_MAP_FILE, show_col_types = FALSE)
  wpid_col <- getcol(wp_map, c("wpid","id","pathway_id","pathway","WPID","WP_ID"))
  sym_col  <- getcol(wp_map, c("SYMBOL","symbol","gene","gene_symbol","GeneSymbol","hgnc_symbol"))
  name_col <- getcol(wp_map, c("name","pathway_name","PathwayName","title","Description"))

  if (is.na(wpid_col) || is.na(sym_col)) {
    stop("Pathway2Gene.csv missing required columns. Found: ", paste(names(wp_map), collapse = ", "))
  }

  wp_sets_df <- wp_map |>
    dplyr::transmute(
      set_id = as.character(.data[[wpid_col]]),
      gene   = as.character(.data[[sym_col]])
    ) |>
    dplyr::filter(gene %in% names(ranks))

  wp_pathways <- make_pathways_list(wp_sets_df, MIN_SET_SIZE, MAX_SET_SIZE)

  set.seed(SET_SEED)
  wp_fg <- fgsea::fgseaMultilevel(
    pathways = wp_pathways,
    stats    = ranks,
    minSize  = MIN_SET_SIZE,
    maxSize  = MAX_SET_SIZE
  )

  wp_out <- wp_fg |>
    tibble::as_tibble() |>
    dplyr::transmute(
      ID = pathway,
      NES, pval, padj, size,
      leadingEdge = purrr::map_chr(leadingEdge, leading_edge_to_string)
    ) |>
    dplyr::arrange(padj, dplyr::desc(abs(NES)))

  if (!is.na(name_col)) {
    wp_names <- wp_map |>
      dplyr::transmute(
        ID          = as.character(.data[[wpid_col]]),
        Description = as.character(.data[[name_col]])
      ) |>
      dplyr::distinct(ID, Description)

    wp_out <- wp_out |>
      dplyr::left_join(wp_names, by = "ID") |>
      dplyr::relocate(Description, .after = ID)
  }

  FG_WP_FILE <- file.path(paths$step02, "fgsea_wikipathways.csv")
  readr::write_csv(wp_out, FG_WP_FILE)

  # ---- 5) fgsea on GO:BP -------------------------------------------------
  go_map_raw <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = names(ranks),
      keytype = "SYMBOL",
      columns = c("GO", "ONTOLOGY")
    )
  )

  if (!all(c("GO", "SYMBOL") %in% names(go_map_raw))) {
    stop("org.Hs.eg.db::select did not return expected columns. Got: ",
         paste(names(go_map_raw), collapse = ", "))
  }
  if (!"ONTOLOGY" %in% names(go_map_raw)) go_map_raw$ONTOLOGY <- NA_character_

  go_sets_df <- go_map_raw |>
    tibble::as_tibble() |>
    dplyr::transmute(
      set_id = as.character(GO),
      gene   = as.character(SYMBOL),
      ont    = as.character(ONTOLOGY)
    ) |>
    dplyr::filter(!is.na(set_id), !is.na(gene)) |>
    dplyr::filter(is.na(ont) | ont == "BP") |>
    dplyr::filter(gene %in% names(ranks))

  go_pathways <- make_pathways_list(
    dplyr::select(go_sets_df, set_id, gene),
    MIN_SET_SIZE, MAX_SET_SIZE
  )

  set.seed(SET_SEED)
  go_fg <- fgsea::fgseaMultilevel(
    pathways = go_pathways,
    stats    = ranks,
    minSize  = MIN_SET_SIZE,
    maxSize  = MAX_SET_SIZE
  )

  go_out <- go_fg |>
    tibble::as_tibble() |>
    dplyr::transmute(
      ID = pathway,
      NES, pval, padj, size,
      leadingEdge = purrr::map_chr(leadingEdge, leading_edge_to_string)
    ) |>
    dplyr::arrange(padj, dplyr::desc(abs(NES)))

  FG_GO_FILE <- file.path(paths$step02, "fgsea_go_bp.csv")
  readr::write_csv(go_out, FG_GO_FILE)

  # ---- 6) Run note -------------------------------------------------------
  note <- glue::glue("
[Pipeline B \u2014 fgsea sensitivity run note]
Expression: {EXPR_FILE}
Metadata:   {META_FILE}
Universe:   {UNIVERSE_FILE}

Samples used: {ncol(expr_mat_u)}
Genes used (after universe restriction + duplicate collapse): {nrow(expr_mat_u)}

Design: {paste(colnames(design), collapse = ', ')}
Rank statistic: limma moderated t-statistics (coef = {geno_coef})
Interpretation note: NES sign follows the sign of the ranking statistic.

Gene set size filter: {MIN_SET_SIZE}\u2013{MAX_SET_SIZE}

WikiPathways sets: {length(wp_pathways)}
GO:BP sets: {length(go_pathways)}

Outputs:
- {FG_WP_FILE}
- {FG_GO_FILE}
")

  writeLines(note, file.path(paths$step02, "run_note.txt"))
  cat(note, "\nDone.\n")

  invisible(list(
    fgsea_wp_file = FG_WP_FILE,
    fgsea_go_file = FG_GO_FILE
  ))
}
