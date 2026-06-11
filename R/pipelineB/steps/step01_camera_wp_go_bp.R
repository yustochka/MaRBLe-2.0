# R/pipelineB/steps/step01_camera_wp_go_bp.R
# Pipeline B — Step 01: CAMERA enrichment (correlation-adjusted competitive test)
#
# Reads:
#   cfg$dataset$expression_file       (e.g. data/dataset_0/processed/processedMatrix.csv)
#   cfg$dataset$metadata_file         (e.g. data/dataset_0/processed/metaData.csv)
#   cfg$dataset$pathway2gene_file     (e.g. data/dataset_0/processed/Pathway2Gene.csv)
#
# Writes into paths$step01  (results/pipelineB/<run_id>/01_camera/):
#   universe_genes.txt                genes in expression matrix covered by any pathway
#   camera_wikipathways.csv
#   camera_go_bp.csv
#   run_note.txt
#
# Returns (invisibly):
#   universe_file    path to universe_genes.txt
#   camera_wp_file   path to camera_wikipathways.csv
#   camera_go_file   path to camera_go_bp.csv

suppressPackageStartupMessages({
  library(limma)
  library(org.Hs.eg.db)
})

if (!exists("convert_expr_to_symbol")) source("R/utils/gene_id_convert.R")

step01_camera_wp_go_bp <- function(cfg, paths, ctx) {

  # ---- Parameters (from cfg with script-matching defaults) ---------------
  EXPR_FILE     <- cfg$dataset$expression_file
  META_FILE     <- cfg$dataset$metadata_file
  WP_MAP_FILE   <- cfg$dataset$pathway2gene_file
  TARGET_TISSUE <- if (!is.null(cfg$pipelineB$target_tissue)) cfg$pipelineB$target_tissue else "Temporal cortex"
  MIN_SET_SIZE  <- if (!is.null(cfg$pipelineB$min_set_size))  as.integer(cfg$pipelineB$min_set_size)  else 10L
  MAX_SET_SIZE  <- if (!is.null(cfg$pipelineB$max_set_size))  as.integer(cfg$pipelineB$max_set_size)  else 500L
  INCLUDE_AGE   <- if (!is.null(cfg$pipelineB$include_age))   isTRUE(cfg$pipelineB$include_age)       else FALSE

  stopifnot(file.exists(EXPR_FILE), file.exists(META_FILE), file.exists(WP_MAP_FILE))

  # ---- Helpers -----------------------------------------------------------
  read_gene_vec <- function(path) {
    x <- readr::read_lines(path)
    x <- stringr::str_trim(x)
    unique(x[x != ""])
  }

  make_index_list <- function(set_df, genes_in_expr) {
    set_df |>
      dplyr::filter(!is.na(set_id), !is.na(gene)) |>
      dplyr::distinct(set_id, gene) |>
      dplyr::group_by(set_id) |>
      dplyr::summarise(genes = list(unique(gene)), .groups = "drop") |>
      dplyr::mutate(
        idx  = purrr::map(genes, ~ unname(genes_in_expr[.x])),
        idx  = purrr::map(idx,   ~ .x[!is.na(.x)]),
        size = purrr::map_int(idx, length)
      ) |>
      dplyr::filter(size >= MIN_SET_SIZE, size <= MAX_SET_SIZE) |>
      dplyr::select(set_id, idx, size)
  }

  getcol <- function(df, candidates) {
    nm <- names(df); nml <- tolower(nm)
    hit <- nm[match(tolower(candidates), nml, nomatch = 0)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }

  # ---- 1) Read expression + metadata -------------------------------------
  expr_raw <- readr::read_csv(EXPR_FILE, show_col_types = FALSE)
  meta     <- readr::read_csv(META_FILE, show_col_types = FALSE)

  gene_col   <- names(expr_raw)[1]
  expr_genes <- as.character(expr_raw[[gene_col]])

  expr_mat <- expr_raw |>
    dplyr::select(-dplyr::all_of(gene_col)) |>
    as.matrix()

  mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- expr_genes

  expr_samples <- colnames(expr_mat)

  # ---- Gene ID conversion (ENTREZ/ENSEMBL → SYMBOL if configured) ---------
  expr_mat <- convert_expr_to_symbol(expr_mat, cfg, out_dir = paths$step01)

  # ---- Universe generation -----------------------------------------------
  # Universe = gene symbols in expression matrix covered by ≥1 pathway.
  wp_map_u      <- readr::read_csv(WP_MAP_FILE, show_col_types = FALSE)
  sym_col_u     <- getcol(wp_map_u, c("SYMBOL", "symbol", "GeneSymbol", "gene_symbol",
                                       "Gene", "gene", "HGNC_symbol", "HUGO"))
  if (is.na(sym_col_u)) {
    stop("Pathway2Gene.csv: cannot detect gene symbol column. Found: ",
         paste(names(wp_map_u), collapse = ", "))
  }
  pathway_genes <- unique(as.character(wp_map_u[[sym_col_u]]))
  universe      <- intersect(rownames(expr_mat), pathway_genes)
  universe      <- universe[!is.na(universe) & universe != ""]
  UNIVERSE_FILE <- file.path(paths$step01, "universe_genes.txt")
  writeLines(universe, UNIVERSE_FILE)
  cat("Universe size:", length(universe), "genes ->", UNIVERSE_FILE, "\n")

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
    meta_use <- meta_use |>
      dplyr::filter(Tissue == TARGET_TISSUE)
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

  # ---- 2) Restrict to universe -------------------------------------------
  keep_genes <- intersect(rownames(expr_mat), universe)
  expr_mat_u <- expr_mat[keep_genes, , drop = FALSE]

  if (nrow(expr_mat_u) < 1000) {
    warning("Very few genes left after universe restriction (", nrow(expr_mat_u), "). ",
            "Check ID type (SYMBOL vs ENSEMBL).")
  }

  gene_to_row <- seq_len(nrow(expr_mat_u))
  names(gene_to_row) <- rownames(expr_mat_u)

  # ---- 3) Fit limma model ------------------------------------------------
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

  fit <- limma::lmFit(expr_mat_u, design)
  fit <- limma::eBayes(fit)

  # ---- 4) CAMERA on WikiPathways -----------------------------------------
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
    dplyr::filter(gene %in% names(gene_to_row))

  wp_index_tbl <- make_index_list(wp_sets_df, gene_to_row)
  wp_index     <- wp_index_tbl$idx
  names(wp_index) <- wp_index_tbl$set_id

  wp_cam <- limma::camera(expr_mat_u, index = wp_index, design = design, contrast = geno_coef[1])

  wp_out <- wp_cam |>
    tibble::rownames_to_column("ID") |>
    tibble::as_tibble()

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

  CAM_WP_FILE <- file.path(paths$step01, "camera_wikipathways.csv")
  readr::write_csv(wp_out, CAM_WP_FILE)

  # ---- 5) CAMERA on GO:BP ------------------------------------------------
  go_map_raw <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = rownames(expr_mat_u),
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
    dplyr::filter(gene %in% names(gene_to_row))

  go_index_tbl <- make_index_list(
    dplyr::select(go_sets_df, set_id, gene),
    gene_to_row
  )
  go_index <- go_index_tbl$idx
  names(go_index) <- go_index_tbl$set_id

  go_cam <- limma::camera(expr_mat_u, index = go_index, design = design, contrast = geno_coef[1])

  go_out <- go_cam |>
    tibble::rownames_to_column("ID") |>
    tibble::as_tibble()

  CAM_GO_FILE <- file.path(paths$step01, "camera_go_bp.csv")
  readr::write_csv(go_out, CAM_GO_FILE)

  # ---- 6) Run note -------------------------------------------------------
  note <- glue::glue("
[Pipeline B \u2014 CAMERA run note]
Expression: {EXPR_FILE}
Universe:   {UNIVERSE_FILE}
Samples used: {ncol(expr_mat_u)}
Genes used (after universe restriction): {nrow(expr_mat_u)}

Design: {paste(colnames(design), collapse = ', ')}
Contrast coefficient: {geno_coef[1]}
Include Age covariate: {INCLUDE_AGE}

WikiPathways sets kept (size {MIN_SET_SIZE}-{MAX_SET_SIZE}): {length(wp_index)}
GO:BP sets kept (size {MIN_SET_SIZE}-{MAX_SET_SIZE}): {length(go_index)}

Outputs:
- {CAM_WP_FILE}
- {CAM_GO_FILE}
")

  writeLines(note, file.path(paths$step01, "run_note.txt"))
  cat(note, "\nDone.\n")

  invisible(list(
    universe_file  = UNIVERSE_FILE,
    camera_wp_file = CAM_WP_FILE,
    camera_go_file = CAM_GO_FILE
  ))
}
