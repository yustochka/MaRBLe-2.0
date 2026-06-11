# R/pipelineB/steps/step05_overlap_cluster_B.R
# Pipeline B — Step 05: Overlap-based clustering (Jaccard graph) on agreement-filtered
# pathway lists; pick 1 representative per cluster with a "biologist-safe" ranking rule.
#
# Reads (from current run):
#   paths$step03/wp_input_agreement.csv   (WP: CAMERA ∩ fgsea agreement list)
#   paths$step04/reps.csv                 (GO:BP: semantically collapsed representatives)
#   paths$step01/universe_genes.txt
#   cfg$dataset$pathway2gene_file     (e.g. data/dataset_0/processed/Pathway2Gene.csv)
#
# NOTE: GO terms come from step04 (not step03 directly) so that semantic
# collapse is applied before overlap clustering. The ablation runner writes a
# step04 passthrough when use_semantic_collapse = false, which copies the raw
# step03 GO list to step04/reps.csv — so this path is always valid.
#
# Writes into paths$step05/  (results/pipelineB/<run_id>/05_overlap/):
#   edges_up.csv
#   edges_down.csv
#   clusters_up.csv
#   clusters_down.csv
#   representatives_up.csv
#   representatives_down.csv
#   run_note.txt
#
# Returns (invisibly):
#   out_dir               path to the step05 output folder
#   representatives_up    path to representatives_up.csv
#   representatives_down  path to representatives_down.csv

suppressPackageStartupMessages({
  library(igraph)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

step05_overlap_cluster_B <- function(cfg, paths, ctx) {

  # ---- Parameters -----------------------------------------------------------
  WP_IN_FILE    <- file.path(paths$step03, "wp_input_agreement.csv")
  GO_IN_FILE    <- file.path(paths$step04, "reps.csv")   # semantically collapsed GO reps
  UNIVERSE_FILE <- file.path(paths$step01, "universe_genes.txt")
  WP_MAP_FILE   <- cfg$dataset$pathway2gene_file

  JACCARD_MIN      <- if (!is.null(cfg$thresholds$JACCARD_MIN))      as.numeric(cfg$thresholds$JACCARD_MIN)      else 0.20
  MIN_INTERSECTION <- if (!is.null(cfg$thresholds$MIN_INTERSECTION)) as.integer(cfg$thresholds$MIN_INTERSECTION) else 3L
  MIN_SET_SIZE     <- if (!is.null(cfg$pipelineB$min_set_size))     as.integer(cfg$pipelineB$min_set_size)     else 10L
  MAX_SET_SIZE     <- if (!is.null(cfg$pipelineB$max_set_size))     as.integer(cfg$pipelineB$max_set_size)     else 500L

  OUT_DIR <- paths$step05

  stopifnot(
    file.exists(WP_IN_FILE), file.exists(GO_IN_FILE),
    file.exists(UNIVERSE_FILE), file.exists(WP_MAP_FILE)
  )

  # ---- Helpers --------------------------------------------------------------
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

  make_sets <- function(df, id_col = "ID", gene_col = "SYMBOL") {
    tmp <- df |>
      dplyr::filter(!is.na(.data[[id_col]]), !is.na(.data[[gene_col]]), .data[[gene_col]] != "") |>
      dplyr::distinct(.data[[id_col]], .data[[gene_col]]) |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::summarise(genes = list(unique(.data[[gene_col]])), .groups = "drop") |>
      dplyr::mutate(size = lengths(genes)) |>
      dplyr::filter(size >= MIN_SET_SIZE, size <= MAX_SET_SIZE)

    if (nrow(tmp) == 0) return(stats::setNames(list(), character()))
    stats::setNames(tmp$genes, tmp[[id_col]])
  }

  pairwise_edges <- function(sets, min_j, min_k) {
    empty <- tibble::tibble(from = character(), to = character(), jaccard = double(), k = integer())
    ids <- names(sets)
    if (length(ids) < 2L) return(empty)

    comb <- t(combn(ids, 2))
    res <- purrr::map_dfr(seq_len(nrow(comb)), function(i) {
      a <- comb[i, 1]; b <- comb[i, 2]
      A <- sets[[a]];   B <- sets[[b]]
      inter <- length(intersect(A, B))
      if (inter < min_k) return(NULL)
      uni <- length(union(A, B))
      if (uni == 0) return(NULL)
      j <- inter / uni
      if (j < min_j) return(NULL)
      tibble::tibble(from = a, to = b, jaccard = j, k = inter)
    })

    if (nrow(res) == 0) empty else res
  }

  # Representative ranking rule:
  # 1) minimize agreement_q = max(CAMERA_FDR, fgsea_padj)
  # 2) maximize |NES|  (if available; otherwise it becomes NA and is ignored)
  # 3) prefer smaller set_size (avoid overly broad generic term)
  rank_and_pick_rep <- function(terms_tbl) {
    terms_tbl |>
      dplyr::mutate(
        agreement_q = pmax(camera_FDR, fgsea_padj, na.rm = TRUE),
        absNES = if ("fgsea_NES" %in% names(terms_tbl)) abs(fgsea_NES) else NA_real_
      ) |>
      dplyr::arrange(agreement_q, dplyr::desc(absNES), set_size, nchar(Description)) |>
      dplyr::slice(1)
  }

  # ---- 1) Read inputs -------------------------------------------------------
  universe <- read_gene_vec(UNIVERSE_FILE)

  wp_in <- readr::read_csv(WP_IN_FILE, show_col_types = FALSE)
  go_in <- readr::read_csv(GO_IN_FILE, show_col_types = FALSE)

  # ---- 1a) Normalize WP input columns --------------------------------------
  dir_col_wp  <- getcol(wp_in, c("direction","Direction"))
  cam_col_wp  <- getcol(wp_in, c("camera_FDR","cam_FDR","cam_fdr","FDR","camera_fdr"))
  pad_col_wp  <- getcol(wp_in, c("fgsea_padj","fg_FDR","fg_fdr","padj","fgseaPadj"))
  nes_col_wp  <- getcol(wp_in, c("fgsea_NES","NES","fgseaNES"))
  desc_col_wp <- getcol(wp_in, c("Description","description","name","PathwayName","title"))

  if (is.na(dir_col_wp) || is.na(cam_col_wp) || is.na(pad_col_wp) || is.na(desc_col_wp)) {
    stop("WP agreement input missing columns. Found: ", paste(names(wp_in), collapse = ", "))
  }

  wp_terms <- wp_in |>
    dplyr::transmute(
      ID         = as.character(ID),
      Description = as.character(.data[[desc_col_wp]]),
      direction  = stringr::str_to_title(as.character(.data[[dir_col_wp]])),
      camera_FDR = as.numeric(.data[[cam_col_wp]]),
      fgsea_padj = as.numeric(.data[[pad_col_wp]]),
      fgsea_NES  = if (!is.na(nes_col_wp)) as.numeric(.data[[nes_col_wp]]) else NA_real_,
      collection = "WP"
    ) |>
    dplyr::filter(!is.na(ID), !is.na(direction))

  # ---- 1b) Normalize GO input columns --------------------------------------
  dir_col_go  <- getcol(go_in, c("direction","Direction"))
  cam_col_go  <- getcol(go_in, c("camera_FDR","cam_FDR","cam_fdr","FDR","camera_fdr"))
  pad_col_go  <- getcol(go_in, c("fgsea_padj","fg_FDR","fg_fdr","padj","fgseaPadj"))
  nes_col_go  <- getcol(go_in, c("fgsea_NES","NES","fgseaNES","fg_NES"))  # fg_NES = step04 column name
  desc_col_go <- getcol(go_in, c("Description","description","term","name"))

  if (is.na(dir_col_go) || is.na(cam_col_go) || is.na(pad_col_go) || is.na(desc_col_go)) {
    stop("GO agreement input missing columns. Found: ", paste(names(go_in), collapse = ", "))
  }

  go_terms <- go_in |>
    dplyr::transmute(
      ID         = as.character(ID),
      Description = as.character(.data[[desc_col_go]]),
      direction  = stringr::str_to_title(as.character(.data[[dir_col_go]])),
      camera_FDR = as.numeric(.data[[cam_col_go]]),
      fgsea_padj = as.numeric(.data[[pad_col_go]]),
      fgsea_NES  = if (!is.na(nes_col_go)) as.numeric(.data[[nes_col_go]]) else NA_real_,
      collection = "GO"
    ) |>
    dplyr::filter(!is.na(ID), !is.na(direction))

  # ---- Early exit: 0 terms from step03 ------------------------------------
  if (nrow(wp_terms) == 0 && nrow(go_terms) == 0) {
    cat("Step B5: 0 terms from step03 — writing canonical empty files and returning.\n")
    empty_reps  <- tibble::tibble(
      ID=character(), Description=character(), direction=character(),
      camera_FDR=double(), fgsea_padj=double(), fgsea_NES=double(),
      collection=character(), set_size=integer(),
      agreement_q=double(), absNES=double(), cluster=integer()
    )
    empty_edges <- tibble::tibble(from=character(), to=character(), jaccard=double(), k=integer())
    REPS_UP_FILE   <- file.path(OUT_DIR, "representatives_up.csv")
    REPS_DOWN_FILE <- file.path(OUT_DIR, "representatives_down.csv")
    for (fname in c("edges_up.csv", "edges_down.csv"))
      readr::write_csv(empty_edges, file.path(OUT_DIR, fname))
    for (fname in c("clusters_up.csv", "clusters_down.csv"))
      readr::write_csv(empty_reps, file.path(OUT_DIR, fname))
    readr::write_csv(empty_reps, REPS_UP_FILE)
    readr::write_csv(empty_reps, REPS_DOWN_FILE)
    writeLines("[Step B5: 0 terms in — 0 kept]", file.path(OUT_DIR, "run_note.txt"))
    return(invisible(list(
      out_dir              = OUT_DIR,
      representatives_up   = REPS_UP_FILE,
      representatives_down = REPS_DOWN_FILE
    )))
  }

  # ---- 2) Build gene sets (restricted to universe) --------------------------
  wp_map   <- readr::read_csv(WP_MAP_FILE, show_col_types = FALSE)
  wpid_col <- getcol(wp_map, c("wpid","id","pathway_id","pathway","WPID","WP_ID"))
  sym_col  <- getcol(wp_map, c("SYMBOL","symbol","gene","gene_symbol","GeneSymbol","hgnc_symbol"))
  if (is.na(wpid_col) || is.na(sym_col)) {
    stop("Pathway2Gene.csv missing ID/SYMBOL columns. Found: ", paste(names(wp_map), collapse = ", "))
  }

  wp_sets_df <- wp_map |>
    dplyr::transmute(
      ID     = as.character(.data[[wpid_col]]),
      SYMBOL = as.character(.data[[sym_col]])
    ) |>
    dplyr::filter(SYMBOL %in% universe, ID %in% wp_terms$ID)

  wp_sets <- make_sets(wp_sets_df, id_col = "ID", gene_col = "SYMBOL")

  go_ids_to_query <- unique(go_terms$ID)
  go_map_raw <- if (length(go_ids_to_query) == 0) {
    tibble::tibble(GO = character(), SYMBOL = character(), ONTOLOGY = character())
  } else {
    suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys    = go_ids_to_query,
        keytype = "GO",
        columns = c("SYMBOL", "ONTOLOGY", "GO")
      )
    )
  }

  id_col_go_map <- intersect(c("GO","GOID","GOALL"), names(go_map_raw))
  if (length(id_col_go_map) == 0) {
    stop("GO mapping did not return GO id column. Got: ", paste(names(go_map_raw), collapse = ", "))
  }
  if (!"ONTOLOGY" %in% names(go_map_raw)) go_map_raw$ONTOLOGY <- NA_character_

  go_sets_df <- go_map_raw |>
    tibble::as_tibble() |>
    dplyr::transmute(
      ID       = as.character(.data[[id_col_go_map[1]]]),
      SYMBOL   = as.character(SYMBOL),
      ONTOLOGY = as.character(ONTOLOGY)
    ) |>
    dplyr::filter(!is.na(ID), !is.na(SYMBOL)) |>
    dplyr::filter(is.na(ONTOLOGY) | ONTOLOGY == "BP") |>
    dplyr::filter(SYMBOL %in% universe, ID %in% go_terms$ID) |>
    dplyr::distinct(ID, SYMBOL)

  go_sets <- make_sets(go_sets_df, id_col = "ID", gene_col = "SYMBOL")

  wp_terms <- wp_terms |>
    dplyr::mutate(set_size = purrr::map_int(ID, ~ length(wp_sets[[.x]]))) |>
    dplyr::filter(set_size >= MIN_SET_SIZE, set_size <= MAX_SET_SIZE)

  go_terms <- go_terms |>
    dplyr::mutate(set_size = purrr::map_int(ID, ~ length(go_sets[[.x]]))) |>
    dplyr::filter(set_size >= MIN_SET_SIZE, set_size <= MAX_SET_SIZE)

  wp_sets <- wp_sets[names(wp_sets) %in% wp_terms$ID]
  go_sets <- go_sets[names(go_sets) %in% go_terms$ID]

  # ---- 3) Direction-aware clustering (Up and Down separately) ---------------
  cluster_one_direction <- function(direction_label) {
    terms_dir <- dplyr::bind_rows(
      wp_terms |> dplyr::filter(direction == direction_label),
      go_terms |> dplyr::filter(direction == direction_label)
    )

    if (nrow(terms_dir) == 0) {
      empty_reps <- terms_dir |>
        dplyr::mutate(cluster = integer(), agreement_q = double(), absNES = double())
      return(list(
        edges      = tibble::tibble(from=character(), to=character(), jaccard=double(), k=integer()),
        clusters   = dplyr::mutate(terms_dir, cluster = integer()),
        reps       = empty_reps,
        n_clusters = 0L
      ))
    }

    ids      <- terms_dir$ID
    sets_dir <- c(wp_sets[ids], go_sets[ids])

    edges <- pairwise_edges(sets_dir, JACCARD_MIN, MIN_INTERSECTION)

    vs <- tibble::tibble(name = terms_dir$ID) |> dplyr::distinct()
    g  <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vs)

    cl   <- if (igraph::ecount(g) > 0) igraph::cluster_louvain(g) else igraph::components(g)
    memb <- if (inherits(cl, "communities")) igraph::membership(cl) else cl$membership

    cluster_df <- tibble::tibble(ID = names(memb), cluster = as.integer(memb))

    clusters_full <- terms_dir |>
      dplyr::left_join(cluster_df, by = "ID") |>
      dplyr::mutate(cluster = ifelse(is.na(cluster), dplyr::row_number() + 10^6, cluster)) |>
      dplyr::arrange(cluster, camera_FDR, fgsea_padj)

    reps <- clusters_full |>
      dplyr::group_by(cluster) |>
      dplyr::group_modify(~ rank_and_pick_rep(.x)) |>
      dplyr::ungroup()

    list(
      edges      = edges,
      clusters   = clusters_full,
      reps       = reps,
      n_clusters = dplyr::n_distinct(clusters_full$cluster)
    )
  }

  up_res   <- cluster_one_direction("Up")
  down_res <- cluster_one_direction("Down")

  # ---- 4) Write outputs -----------------------------------------------------
  readr::write_csv(up_res$edges,   file.path(OUT_DIR, "edges_up.csv"))
  readr::write_csv(down_res$edges, file.path(OUT_DIR, "edges_down.csv"))

  readr::write_csv(up_res$clusters,   file.path(OUT_DIR, "clusters_up.csv"))
  readr::write_csv(down_res$clusters, file.path(OUT_DIR, "clusters_down.csv"))

  REPS_UP_FILE   <- file.path(OUT_DIR, "representatives_up.csv")
  REPS_DOWN_FILE <- file.path(OUT_DIR, "representatives_down.csv")

  # Fill any missing GO descriptions before writing representatives
  if (!exists("add_go_term_name", mode = "function"))
    source("R/utils/go_term_names.R")
  reps_up_out   <- add_go_term_name(up_res$reps,   id_col = "ID", desc_col = "Description")
  reps_down_out <- add_go_term_name(down_res$reps, id_col = "ID", desc_col = "Description")

  readr::write_csv(reps_up_out,   REPS_UP_FILE)
  readr::write_csv(reps_down_out, REPS_DOWN_FILE)

  note <- glue::glue("
[Pipeline B \u2014 Step B5: Overlap clustering]
Direction-aware clustering (Up and Down separately)

Gene set filters:
- size {MIN_SET_SIZE}\u2013{MAX_SET_SIZE}
Overlap edge rule:
- Jaccard \u2265 {JACCARD_MIN}
- shared genes \u2265 {MIN_INTERSECTION}

Representative selection rule (biologist-safe):
1) agreement_q = max(CAMERA_FDR, fgsea_padj) (smaller is better)
2) maximize |NES| (if NES exists in the input; otherwise ignored)
3) prefer smaller gene set size (avoid overly broad generic terms)

Counts:
- WP terms (after size filter): {nrow(wp_terms)}
- GO terms (after size filter): {nrow(go_terms)}

UP:
- clusters: {up_res$n_clusters}
- reps: {nrow(up_res$reps)}
DOWN:
- clusters: {down_res$n_clusters}
- reps: {nrow(down_res$reps)}

Outputs in:
- {OUT_DIR}
")

  writeLines(note, file.path(OUT_DIR, "run_note.txt"))
  cat(note, "\nDone.\n")

  invisible(list(
    out_dir              = OUT_DIR,
    representatives_up   = REPS_UP_FILE,
    representatives_down = REPS_DOWN_FILE
  ))
}
