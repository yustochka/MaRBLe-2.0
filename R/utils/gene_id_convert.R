# R/utils/gene_id_convert.R
# Gene identifier conversion for Pipeline B.
#
# Converts expression matrix row names to HGNC gene symbols using org.Hs.eg.db.
# Controlled entirely by cfg$dataset$gene_id_type.
#
# Supported gene_id_type values (case-insensitive):
#   "SYMBOL"   — no conversion; row names are already gene symbols (default)
#   "ENTREZ"   — Entrez Gene IDs (integer-like strings, e.g. "653635")
#   "ENSEMBL"  — Ensembl gene IDs (e.g. "ENSG00000099998")
#
# Duplicate handling (both directions):
#   1-to-many (one input ID → multiple symbols): keep first symbol returned by db
#   many-to-1 (multiple input IDs → same symbol): keep row with highest mean expression
#
# Prints one log line:
#   "Gene ID conversion: ENTREZ->SYMBOL. Input: X, mapped: Y, dropped: D, duplicates collapsed: C, unique symbols kept: Z."
# If out_dir is provided, writes gene_id_mapping.csv (original_id, symbol, mapped) there.

convert_expr_to_symbol <- function(expr_mat, cfg, out_dir = NULL) {

  gene_id_type <- if (!is.null(cfg$dataset$gene_id_type) &&
                      nzchar(cfg$dataset$gene_id_type))
    toupper(trimws(cfg$dataset$gene_id_type)) else "SYMBOL"

  if (gene_id_type == "SYMBOL") return(expr_mat)   # no-op for dataset_0

  keytype <- switch(gene_id_type,
    "ENTREZ"   = "ENTREZID",
    "ENTREZID" = "ENTREZID",
    "ENSEMBL"  = "ENSEMBL",
    stop("Unsupported gene_id_type: '", gene_id_type, "'.",
         " Use SYMBOL, ENTREZ, or ENSEMBL.")
  )

  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
  })

  input_ids <- rownames(expr_mat)
  n_input   <- length(input_ids)

  map_raw <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = unique(input_ids),
      keytype = keytype,
      columns = "SYMBOL"
    )
  )

  map_raw <- map_raw[!is.na(map_raw$SYMBOL) & nzchar(map_raw$SYMBOL), ]

  # 1-to-many: one input ID maps to multiple symbols → keep first
  map_clean <- map_raw[!duplicated(map_raw[[keytype]]), ]
  id_to_sym <- stats::setNames(map_clean$SYMBOL, map_clean[[keytype]])

  new_names <- id_to_sym[input_ids]
  keep      <- !is.na(new_names)
  n_mapped  <- sum(keep)
  n_dropped <- n_input - n_mapped

  expr_conv <- expr_mat[keep, , drop = FALSE]
  rownames(expr_conv) <- new_names[keep]

  # many-to-1: multiple input IDs map to same symbol → keep highest-mean row
  n_dup_collapsed <- 0L
  if (anyDuplicated(rownames(expr_conv))) {
    row_means <- rowMeans(expr_conv, na.rm = TRUE)
    best_idx  <- tapply(seq_along(row_means), rownames(expr_conv),
                        function(ix) ix[which.max(row_means[ix])])
    n_dup_collapsed <- nrow(expr_conv) - length(best_idx)
    expr_conv <- expr_conv[unname(unlist(best_idx)), , drop = FALSE]
  }

  n_symbols <- nrow(expr_conv)
  cat(sprintf(
    "Gene ID conversion: %s->SYMBOL. Input: %d, mapped: %d, dropped: %d, duplicates collapsed: %d, unique symbols kept: %d.\n",
    gene_id_type, n_input, n_mapped, n_dropped, n_dup_collapsed, n_symbols
  ))

  # Save mapping CSV for transparency/debugging if out_dir is provided
  if (!is.null(out_dir) && nzchar(out_dir)) {
    mapping_df <- data.frame(
      original_id = input_ids,
      symbol      = as.character(new_names),
      mapped      = keep,
      stringsAsFactors = FALSE
    )
    mapping_file <- file.path(out_dir, "gene_id_mapping.csv")
    utils::write.csv(mapping_df, mapping_file, row.names = FALSE, quote = TRUE)
    cat(sprintf("Gene ID mapping saved: %s\n", mapping_file))
  }

  expr_conv
}
