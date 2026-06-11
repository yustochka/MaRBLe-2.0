# R/pipelineA/steps/step06_bootstrap_consensus.R
# Pipeline A — Step 06: Bootstrap stability + cross-collection consensus -> final shortlist
#
# Reads:
#   paths$step01/universe_genes.txt                   (written by step01)
#   paths$step02/de_genes_rawp05.txt                  (written by step02)
#   paths$step05/go_bp/representatives.csv            (written by step05)
#   paths$step05/wp/representatives.csv               (written by step05)
#   cfg$dataset$pathway2gene_file
#
# Writes into paths$step06:
#   go_bp_boot_matrix.csv
#   go_bp_stability.csv
#   wp_boot_matrix.csv
#   wp_stability.csv
#   cross_collection_pairs.csv
#   final_shortlist.csv
#   run_note.txt
#
# Returns:
#   go_stability_table    path to go_bp_stability.csv
#   wp_stability_table    path to wp_stability.csv
#   consensus_pairs_table path to cross_collection_pairs.csv
#   final_shortlist_table path to final_shortlist.csv

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
})

step06_bootstrap_consensus <- function(cfg, paths, ctx) {

  # ---- Parameters (from cfg with script-matching defaults) -------------------
  B             <- if (!is.null(cfg$bootstrap$B))           cfg$bootstrap$B           else 200
  SUBSAMPLE_P   <- if (!is.null(cfg$bootstrap$SUBSAMPLE_P)) cfg$bootstrap$SUBSAMPLE_P else 0.85
  DROP_DELTA    <- 0.10
  FDR_CUT       <- 0.05
  TAU_STABILITY    <- if (!is.null(cfg$thresholds$TAU_STABILITY))    cfg$thresholds$TAU_STABILITY    else 0.70
  CONS_JACCARD_MIN <- if (!is.null(cfg$thresholds$CONS_JACCARD_MIN)) cfg$thresholds$CONS_JACCARD_MIN else 0.10
  CONS_K_MIN       <- if (!is.null(cfg$thresholds$CONS_K_MIN))       cfg$thresholds$CONS_K_MIN       else 3

  # Reset seed for reproducibility within this step
  set.seed(ctx$seed)

  # ---- Inputs ----------------------------------------------------------------
  UNIVERSE_FILE <- file.path(paths$step01, "universe_genes.txt")
  DE_FILE       <- file.path(paths$step02, "de_genes_rawp05.txt")
  GO_REPS_FILE  <- file.path(paths$step05, "go_bp", "representatives.csv")
  WP_REPS_FILE  <- file.path(paths$step05, "wp",    "representatives.csv")
  WP_MAP_FILE   <- cfg$dataset$pathway2gene_file

  stopifnot(
    file.exists(UNIVERSE_FILE), file.exists(DE_FILE),
    file.exists(GO_REPS_FILE),  file.exists(WP_REPS_FILE),
    file.exists(WP_MAP_FILE)
  )

  # ---- Helpers ---------------------------------------------------------------
  read_gene_vec <- function(path) {
    x <- readr::read_lines(path)
    x <- stringr::str_trim(x)
    unique(x[x != ""])
  }

  ora_pval <- function(de_set, term_set, universe_set) {
    x <- length(intersect(de_set, term_set))
    m <- length(term_set)
    N <- length(universe_set)
    K <- length(de_set)
    if (m == 0L || N == 0L || K == 0L) return(1)
    stats::phyper(q = x - 1, m, N - m, K, lower.tail = FALSE)
  }

  boot_stability <- function(reps_tbl, sets_list, universe, de, label) {
    ids  <- reps_tbl$ID
    nset <- length(ids)

    sel_mat <- matrix(0L, nrow = B, ncol = nset,
                      dimnames = list(paste0("b", seq_len(B)), ids))

    for (b in seq_len(B)) {
      u_b  <- sample(universe, size = floor(length(universe) * SUBSAMPLE_P), replace = FALSE)
      de_b <- intersect(de, u_b)

      if (length(de_b) > 0L && DROP_DELTA > 0) {
        drop_n <- floor(length(de_b) * DROP_DELTA)
        if (drop_n > 0L) de_b <- setdiff(de_b, sample(de_b, drop_n))
      }

      pvals <- vapply(ids, function(id) {
        s_b <- intersect(sets_list[[id]], u_b)
        ora_pval(de_b, s_b, u_b)
      }, numeric(1))

      padj <- stats::p.adjust(pvals, method = "BH")
      sel_mat[b, ] <- as.integer(padj <= FDR_CUT)

      if (b %% 25 == 0) cat(glue::glue("[{label}] bootstrap {b}/{B}\n"))
    }

    stability <- colMeans(sel_mat)
    tib <- reps_tbl |>
      dplyr::mutate(stability = stability,
                    selected  = stability >= TAU_STABILITY)

    list(matrix = sel_mat, summary = tib)
  }

  pairwise_jaccard <- function(listA, listB, j_min = 0, k_min = 0) {
    if (!length(listA) || !length(listB)) {
      return(tibble::tibble(A = character(), B = character(),
                            jaccard = numeric(), k = integer()))
    }
    a_ids <- names(listA); b_ids <- names(listB)
    purrr::map_dfr(a_ids, function(a) {
      purrr::map_dfr(b_ids, function(b) {
        inter <- length(intersect(listA[[a]], listB[[b]]))
        uni   <- length(union(listA[[a]], listB[[b]]))
        if (uni == 0) return(NULL)
        j <- inter / uni
        if (inter >= k_min && j >= j_min)
          tibble::tibble(A = a, B = b, jaccard = j, k = inter)
        else NULL
      })
    })
  }

  # ---- 1) Read inputs --------------------------------------------------------
  universe <- read_gene_vec(UNIVERSE_FILE)
  de_genes <- intersect(read_gene_vec(DE_FILE), universe)

  go_reps <- readr::read_csv(GO_REPS_FILE, show_col_types = FALSE) |>
    dplyr::transmute(ID          = as.character(ID),
                     Description = as.character(Description),
                     p.adjust    = as.numeric(p.adjust),
                     k           = as.integer(k))

  wp_reps <- readr::read_csv(WP_REPS_FILE, show_col_types = FALSE) |>
    dplyr::transmute(ID          = as.character(ID),
                     Description = as.character(Description),
                     p.adjust    = as.numeric(p.adjust),
                     k           = as.integer(k))

  # ---- 2) Build fixed gene sets for the representative terms -----------------
  # --- GO-BP mapping ---
  go_map_raw <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = unique(go_reps$ID),
      keytype = "GO",
      columns = c("SYMBOL", "ONTOLOGY")
    )
  )

  id_col <- intersect(c("GO", "GOID", "GOALL"), names(go_map_raw))
  if (length(id_col) == 0) {
    stop("org.Hs.eg.db::select() returned no GO id column (GO/GOID/GOALL). Got: ",
         paste(names(go_map_raw), collapse = ", "))
  }
  if (!"ONTOLOGY" %in% names(go_map_raw)) go_map_raw$ONTOLOGY <- NA_character_

  go_map <- go_map_raw |>
    tibble::as_tibble() |>
    dplyr::rename(ID = !!id_col[1]) |>
    dplyr::mutate(
      ID       = as.character(ID),
      SYMBOL   = as.character(SYMBOL),
      ONTOLOGY = as.character(ONTOLOGY)
    ) |>
    dplyr::filter(!is.na(ID), !is.na(SYMBOL)) |>
    dplyr::filter(is.na(ONTOLOGY) | ONTOLOGY == "BP") |>
    dplyr::select(ID, SYMBOL) |>
    dplyr::filter(SYMBOL %in% universe) |>
    dplyr::distinct()

  go_sets_tbl <- go_map |>
    dplyr::group_by(ID) |>
    dplyr::summarise(genes = list(unique(SYMBOL)), .groups = "drop")

  go_sets <- go_sets_tbl$genes
  names(go_sets) <- go_sets_tbl$ID
  go_sets <- go_sets[names(go_sets) %in% go_reps$ID]

  # --- WikiPathways mapping ---
  wp_map_raw <- readr::read_csv(WP_MAP_FILE, show_col_types = FALSE)

  id_candidates  <- c("wpid", "ID", "id", "Pathway", "pathway", "pathway_id", "WPID", "WP_ID", "WP")
  sym_candidates <- c("SYMBOL", "symbol", "GeneSymbol", "gene_symbol", "Gene", "gene", "HGNC_symbol", "HUGO")

  id_col  <- names(wp_map_raw)[tolower(names(wp_map_raw)) %in% tolower(id_candidates)]
  sym_col <- names(wp_map_raw)[tolower(names(wp_map_raw)) %in% tolower(sym_candidates)]

  if (length(id_col) == 0 || length(sym_col) == 0) {
    stop(
      "WP_MAP_FILE is missing expected columns.\n",
      "Looked for ID in: ",     paste(id_candidates,  collapse = ", "),
      "\nLooked for SYMBOL in: ", paste(sym_candidates, collapse = ", "),
      "\nFound columns: ",        paste(names(wp_map_raw), collapse = ", ")
    )
  }

  wp_map <- wp_map_raw |>
    dplyr::transmute(
      ID     = as.character(.data[[ id_col[1]  ]]),
      SYMBOL = as.character(.data[[ sym_col[1] ]])
    ) |>
    dplyr::filter(SYMBOL %in% universe, ID %in% wp_reps$ID) |>
    dplyr::distinct()

  wp_sets_tbl <- wp_map |>
    dplyr::group_by(ID) |>
    dplyr::summarise(genes = list(unique(SYMBOL)), .groups = "drop")

  wp_sets <- wp_sets_tbl$genes
  names(wp_sets) <- wp_sets_tbl$ID

  # Drop reps with empty gene sets
  go_reps <- go_reps |> dplyr::filter(ID %in% names(go_sets))
  wp_reps <- wp_reps |> dplyr::filter(ID %in% names(wp_sets))

  # ---- 3) Bootstrap stability ------------------------------------------------
  cat("Running bootstrap stability…\n")
  go_boot <- boot_stability(go_reps, go_sets, universe, de_genes, label = "GO")
  wp_boot <- boot_stability(wp_reps, wp_sets, universe, de_genes, label = "WP")

  readr::write_csv(tibble::as_tibble(go_boot$matrix, rownames = "bootstrap"),
                   file.path(paths$step06, "go_bp_boot_matrix.csv"))
  readr::write_csv(go_boot$summary, file.path(paths$step06, "go_bp_stability.csv"))

  readr::write_csv(tibble::as_tibble(wp_boot$matrix, rownames = "bootstrap"),
                   file.path(paths$step06, "wp_boot_matrix.csv"))
  readr::write_csv(wp_boot$summary, file.path(paths$step06, "wp_stability.csv"))

  # ---- 4) Cross-collection consensus (GO <-> WP) -----------------------------
  cat("Computing GO<->WP consensus…\n")
  pairs <- pairwise_jaccard(go_sets, wp_sets, j_min = CONS_JACCARD_MIN, k_min = CONS_K_MIN)
  readr::write_csv(pairs, file.path(paths$step06, "cross_collection_pairs.csv"))
  go_cons_ids <- unique(pairs$A)
  wp_cons_ids <- unique(pairs$B)

  # ---- 5) Final shortlist ----------------------------------------------------
  go_final <- go_boot$summary |>
    dplyr::mutate(collection = "GO",
                  consensus  = ID %in% go_cons_ids,
                  keep_final = selected & consensus)

  wp_final <- wp_boot$summary |>
    dplyr::mutate(collection = "WP",
                  consensus  = ID %in% wp_cons_ids,
                  keep_final = selected & consensus)

  final_shortlist <- dplyr::bind_rows(go_final, wp_final) |>
    dplyr::arrange(dplyr::desc(keep_final), dplyr::desc(stability), p.adjust, collection)

  readr::write_csv(final_shortlist, file.path(paths$step06, "final_shortlist.csv"))

  # ---- 6) Run note -----------------------------------------------------------
  note <- glue::glue("
[Bootstrap + Consensus]
B = {B}, subsample_p = {SUBSAMPLE_P}, drop_delta = {DROP_DELTA}, FDR_cut = {FDR_CUT}
stability_cut = {TAU_STABILITY}
consensus: Jaccard >= {CONS_JACCARD_MIN}, shared genes >= {CONS_K_MIN}

GO reps: {nrow(go_reps)}   WP reps: {nrow(wp_reps)}
GO kept (stable): {sum(go_boot$summary$selected)}  WP kept (stable): {sum(wp_boot$summary$selected)}
GO consensus-supported: {length(go_cons_ids)}      WP consensus-supported: {length(wp_cons_ids)}
Final keepers (stable & consensus): {sum(final_shortlist$keep_final)}
")
  writeLines(note, file.path(paths$step06, "run_note.txt"))
  cat(note, "\nDone -> results in:", paths$step06, "\n")

  invisible(list(
    go_stability_table    = file.path(paths$step06, "go_bp_stability.csv"),
    wp_stability_table    = file.path(paths$step06, "wp_stability.csv"),
    consensus_pairs_table = file.path(paths$step06, "cross_collection_pairs.csv"),
    final_shortlist_table = file.path(paths$step06, "final_shortlist.csv")
  ))
}
