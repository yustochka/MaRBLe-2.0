# R/pipelineA/steps/step05_overlap_clustering.R
# Pipeline A — Step 05: Jaccard overlap clustering -> representatives
#
# Reads:
#   paths$step04/semantic/reps.csv              (GO semantic collapse reps, step04a)
#   paths$step03_wp/ora_wikipathways_filtered.csv (WP ORA filtered, step03a)
#   cfg$dataset$pathway2gene_file
#   paths$step01/universe_genes.txt
#
# Writes into paths$step05/go_bp/ and paths$step05/wp/:
#   representatives.csv    top-ranked term per cluster
#   clusters.csv           all terms with cluster assignments
#   edges.csv              all kept Jaccard edges
#   run_note.txt
#
# Returns:
#   go_reps_file     path to GO representatives.csv
#   wp_reps_file     path to WP representatives.csv
#   go_clusters_file path to GO clusters.csv
#   wp_clusters_file path to WP clusters.csv

suppressPackageStartupMessages({
  library(igraph)
  library(org.Hs.eg.db)
})

step05_overlap_clustering <- function(cfg, paths, ctx) {

  # ---- Knobs ------------------------------------------------------------------
  # These thresholds are intentionally hardcoded at Pipeline A's original values
  # and NOT wired to cfg, to preserve strict pipeline separation.
  # Pipeline B uses cfg$thresholds$JACCARD_MIN (default 0.20) for a looser
  # pre-consensus overlap pass; Pipeline A uses 0.30 here for a tighter single
  # pass that directly feeds bootstrap consensus.  Do not unify these without
  # careful validation that results are unchanged.
  JACCARD_MIN      <- 0.30
  MIN_INTERSECTION <- 3
  MIN_SET_SIZE     <- 10
  MAX_SET_SIZE     <- 500

  # ---- Inputs ----------------------------------------------------------------
  GO_SEM_REPS   <- file.path(paths$step04, "semantic", "reps.csv")
  WP_ORA_FILE   <- file.path(paths$step03_wp, "ora_wikipathways_filtered.csv")
  WP_MAP_FILE   <- cfg$dataset$pathway2gene_file
  UNIVERSE_FILE <- file.path(paths$step01, "universe_genes.txt")

  # Hierarchy tagging not yet refactored; disabled safely
  GO_HIER_REPS <- NULL

  stopifnot(
    file.exists(GO_SEM_REPS),
    file.exists(WP_ORA_FILE),
    file.exists(WP_MAP_FILE),
    file.exists(UNIVERSE_FILE)
  )

  # ---- Output dirs -----------------------------------------------------------
  OUT_GO_DIR <- file.path(paths$step05, "go_bp")
  OUT_WP_DIR <- file.path(paths$step05, "wp")
  dir.create(OUT_GO_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(OUT_WP_DIR, recursive = TRUE, showWarnings = FALSE)

  # ---- Helpers ---------------------------------------------------------------
  read_universe <- function(path) {
    x <- readr::read_lines(path)
    x <- stringr::str_trim(x)
    unique(x[x != ""])
  }

  pairwise_jaccard <- function(sets, min_j, min_k) {
    empty <- tibble::tibble(from = character(), to = character(),
                            jaccard = double(), k = integer())
    ids <- names(sets)
    if (length(ids) < 2L) return(empty)

    comb <- t(combn(ids, 2))
    res <- purrr::map(seq_len(nrow(comb)), function(i) {
      a <- comb[i, 1]; b <- comb[i, 2]
      A <- sets[[a]];   B <- sets[[b]]
      inter <- length(intersect(A, B))
      if (inter < min_k) return(NULL)
      j <- inter / length(union(A, B))
      if (j < min_j) return(NULL)
      tibble::tibble(from = a, to = b, jaccard = j, k = inter)
    }) |> dplyr::bind_rows()

    if (nrow(res) == 0) empty else res
  }

  pick_representatives <- function(terms_tbl, edges, out_dir, label) {
    nm <- names(terms_tbl)

    id_col <- intersect(c("ID", "id", "go_id", "GOID", "goid", "term", "term_id"), nm)
    if (!length(id_col)) stop("terms_tbl: no ID column; has: ", paste(nm, collapse = ", "))
    terms_tbl <- terms_tbl |> dplyr::rename(ID = !!id_col[1])

    desc_col <- intersect(c("Description", "description", "term_name", "name", "title"), names(terms_tbl))
    if (!length(desc_col)) stop("terms_tbl: no description-like column (Description/name).")
    terms_tbl <- terms_tbl |> dplyr::rename(Description = !!desc_col[1])

    if (!"p.adjust" %in% names(terms_tbl)) {
      pa_col <- intersect(c("padj", "p_adjust", "p_adj", "qvalue", "q_value"), names(terms_tbl))
      if (!length(pa_col)) stop("terms_tbl: no FDR-like column (p.adjust/padj/qvalue).")
      terms_tbl <- dplyr::rename(terms_tbl, p.adjust = !!pa_col[1])
    }
    if (!"k" %in% names(terms_tbl)) {
      k_col <- intersect(c("K", "k1", "count", "overlap", "n_in_table_1", "n_de"), names(terms_tbl))
      if (length(k_col)) terms_tbl <- dplyr::rename(terms_tbl, k = !!k_col[1]) else terms_tbl$k <- NA_integer_
    }

    if (is.null(edges) || !all(c("from", "to") %in% names(edges))) {
      edges <- tibble::tibble(from = character(), to = character(),
                              jaccard = double(), k = integer())
    }

    # ensure K exists (pathway size)
    if (!"K" %in% names(terms_tbl)) {
      if ("n" %in% names(terms_tbl)) {
        terms_tbl <- dplyr::rename(terms_tbl, K = n)
      } else {
        terms_tbl$K <- NA_integer_
      }
    }

    vs <- terms_tbl |> dplyr::distinct(ID) |> dplyr::rename(name = ID)
    g  <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vs)

    cl   <- if (igraph::ecount(g) > 0) igraph::cluster_louvain(g) else igraph::components(g)
    memb <- if (inherits(cl, "communities")) igraph::membership(cl) else cl$membership
    cluster_df <- tibble::tibble(ID = names(memb), cluster = as.integer(memb))

    ranked <- terms_tbl |>
      dplyr::left_join(cluster_df, by = "ID") |>
      dplyr::mutate(
        cluster  = ifelse(is.na(cluster), dplyr::row_number() + 10^6, cluster),
        hit_frac = ifelse(is.na(K), NA_real_, k / pmax(K, 1))
      ) |>
      dplyr::group_by(cluster) |>
      dplyr::arrange(
        p.adjust,
        dplyr::desc(hit_frac),
        K,
        dplyr::desc(k),
        nchar(Description),
        .by_group = TRUE
      ) |>
      dplyr::mutate(
        rank_in_cluster = dplyr::row_number(),
        cluster_size    = dplyr::n()
      ) |>
      dplyr::ungroup()

    reps <- ranked |> dplyr::filter(rank_in_cluster == 1L)

    readr::write_csv(reps, file.path(out_dir, "representatives.csv"))
    readr::write_csv(
      ranked |> dplyr::select(ID, Description, p.adjust, k,
                               cluster, rank_in_cluster, cluster_size),
      file.path(out_dir, "clusters.csv")
    )
    readr::write_csv(edges, file.path(out_dir, "edges.csv"))

    note <- glue::glue(
      "[{label}] Jaccard clustering\n",
      "terms_in   = {nrow(terms_tbl)}\n",
      "edges_kept = {nrow(edges)}  (jaccard>={JACCARD_MIN}, k>={MIN_INTERSECTION})\n",
      "clusters   = {length(unique(ranked$cluster))}\n",
      "reps_kept  = {nrow(reps)}\n",
      "filters    = size[{MIN_SET_SIZE},{MAX_SET_SIZE}]"
    )
    writeLines(note, file.path(out_dir, "run_note.txt"))

    reps
  }

  # ---- 0) Universe -----------------------------------------------------------
  universe <- read_universe(UNIVERSE_FILE)

  # ---- 1) GO-BP (semantic collapse reps) ------------------------------------
  go_sem <- readr::read_csv(GO_SEM_REPS, show_col_types = FALSE) |>
    dplyr::mutate(
      ID       = as.character(ID),
      p.adjust = as.numeric(p.adjust),
      k        = as.integer(k)
    )

  if (!is.null(GO_HIER_REPS) && file.exists(GO_HIER_REPS)) {
    go_hier <- readr::read_csv(GO_HIER_REPS, show_col_types = FALSE) |>
      dplyr::mutate(ID = as.character(ID))
    go_sem <- go_sem |>
      dplyr::mutate(supported_by_hierarchy = ID %in% go_hier$ID)
  } else {
    go_sem <- go_sem |> dplyr::mutate(supported_by_hierarchy = NA)
  }

  go_ids <- unique(go_sem$ID)

  go_map_raw <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = go_ids,
      keytype = "GO",
      columns = c("SYMBOL", "ONTOLOGY")
    )
  )

  id_col <- intersect(c("GO", "GOID", "GOALL"), names(go_map_raw))
  if (length(id_col) == 0) {
    stop("Could not find a GO ID column (GO/GOID/GOALL) in org.Hs.eg.db::select() output.")
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
    dplyr::distinct() |>
    dplyr::filter(is.na(ONTOLOGY) | ONTOLOGY == "BP")

  go_sets <- go_map |>
    dplyr::mutate(SYMBOL = ifelse(SYMBOL %in% universe, SYMBOL, NA_character_)) |>
    dplyr::filter(!is.na(SYMBOL)) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(genes = list(unique(SYMBOL)), .groups = "drop") |>
    dplyr::mutate(n = lengths(genes)) |>
    dplyr::filter(ID %in% go_sem$ID, n >= MIN_SET_SIZE, n <= MAX_SET_SIZE)

  go_list  <- stats::setNames(go_sets$genes, go_sets$ID)
  go_edges <- pairwise_jaccard(go_list, JACCARD_MIN, MIN_INTERSECTION)

  # Include all go_sem rows; terms not in go_sets (size < MIN_SET_SIZE or > MAX_SET_SIZE)
  # have no Jaccard edges and become singleton clusters — still valid representatives.
  go_terms_tbl <- go_sem |>
    dplyr::left_join(go_sets |> dplyr::transmute(ID, K = n), by = "ID") |>
    dplyr::select(dplyr::any_of(c("ID", "Description", "p.adjust", "k", "K",
                                  "supported_by_hierarchy")))

  n_go_singletons <- sum(!go_sem$ID %in% go_sets$ID)
  if (n_go_singletons > 0L)
    cat(sprintf("  Keeping %d singleton GO representative(s) (n < %d or n > %d in universe)\n",
                n_go_singletons, MIN_SET_SIZE, MAX_SET_SIZE))

  go_reps <- pick_representatives(go_terms_tbl, go_edges, OUT_GO_DIR, "GO-BP")

  # ---- 2) WikiPathways -------------------------------------------------------
  getcol <- function(df, candidates) {
    nm  <- names(df); nml <- tolower(nm)
    hit <- nm[match(tolower(candidates), nml, nomatch = 0)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }

  wp_ora_raw <- readr::read_csv(WP_ORA_FILE, show_col_types = FALSE)

  idc  <- getcol(wp_ora_raw, c("ID", "wpid", "pathway", "pathway_id", "term_id", "wp_id"))
  desc <- getcol(wp_ora_raw, c("Description", "Term", "Name", "PathwayName", "pathway_name"))
  padj <- getcol(wp_ora_raw, c("p.adjust", "padj", "qvalue", "q_value"))
  kcol <- getcol(wp_ora_raw, c("k", "overlap", "count", "n_de", "n_overlap"))

  if (is.na(idc) || is.na(desc) || is.na(padj)) {
    stop("WP ORA file is missing required columns. Looked for:\n",
         "  ID: ID/wpID/pathway/pathway_id/term_id/wp_id\n",
         "  Description: Description/Term/Name/PathwayName\n",
         "  FDR: p.adjust/padj/qvalue")
  }

  k_vec <- if (!is.na(kcol)) as.integer(wp_ora_raw[[kcol]]) else rep(NA_integer_, nrow(wp_ora_raw))

  wp_ora <- wp_ora_raw |>
    dplyr::transmute(
      ID          = as.character(.data[[idc]]),
      Description = as.character(.data[[desc]]),
      p.adjust    = as.numeric(.data[[padj]]),
      k           = k_vec
    )

  wp_map_raw <- readr::read_csv(WP_MAP_FILE, show_col_types = FALSE)

  idc2 <- getcol(wp_map_raw, c("ID", "wpid", "pathway", "pathway_id", "wp_id"))
  symc <- getcol(wp_map_raw, c("SYMBOL", "gene", "GeneSymbol", "gene_symbol", "hgnc_symbol"))

  if (is.na(idc2) || is.na(symc)) {
    stop("Pathway2Gene.csv is missing required columns. Looked for:\n",
         "  Pathway ID: ID/wpID/pathway/pathway_id\n",
         "  Gene symbol: SYMBOL/gene/GeneSymbol/gene_symbol")
  }

  wp_map <- wp_map_raw |>
    dplyr::transmute(
      ID     = as.character(.data[[idc2]]),
      SYMBOL = as.character(.data[[symc]])
    )

  wp_sets <- wp_map |>
    dplyr::filter(SYMBOL %in% universe) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(genes = list(unique(SYMBOL)), .groups = "drop") |>
    dplyr::mutate(n = lengths(genes)) |>
    dplyr::filter(ID %in% wp_ora$ID, n >= MIN_SET_SIZE, n <= MAX_SET_SIZE)

  wp_list  <- stats::setNames(wp_sets$genes, wp_sets$ID)
  wp_edges <- pairwise_jaccard(wp_list, JACCARD_MIN, MIN_INTERSECTION)

  # Include all wp_ora rows; same singleton logic as GO above.
  wp_terms_tbl <- wp_ora |>
    dplyr::left_join(wp_sets |> dplyr::transmute(ID, K = n), by = "ID") |>
    dplyr::select(dplyr::any_of(c("ID", "Description", "p.adjust", "k", "K")))

  n_wp_singletons <- sum(!wp_ora$ID %in% wp_sets$ID)
  if (n_wp_singletons > 0L)
    cat(sprintf("  Keeping %d singleton WP representative(s) (n < %d or n > %d in universe)\n",
                n_wp_singletons, MIN_SET_SIZE, MAX_SET_SIZE))

  wp_reps <- pick_representatives(wp_terms_tbl, wp_edges, OUT_WP_DIR, "WikiPathways")

  cat("\n== Step 05 done ==\n",
      "GO reps kept: ", nrow(go_reps), "  (", OUT_GO_DIR, ")\n",
      "WP reps kept: ", nrow(wp_reps), "  (", OUT_WP_DIR, ")\n", sep = "")

  invisible(list(
    go_reps_file     = file.path(OUT_GO_DIR, "representatives.csv"),
    wp_reps_file     = file.path(OUT_WP_DIR, "representatives.csv"),
    go_clusters_file = file.path(OUT_GO_DIR, "clusters.csv"),
    wp_clusters_file = file.path(OUT_WP_DIR, "clusters.csv")
  ))
}
