# R/pipelineB/steps/step06_consensus_and_themes_A1A2.R
# Pipeline B — Step 06: Consensus + Themes (single run, cfg-driven).
# Keeps terms if GO<->WP consensus OR pass a solo high-confidence rule.
#
# Reads (from current run):
#   paths$step05/representatives_up.csv
#   paths$step05/representatives_down.csv
#   paths$step01/universe_genes.txt
#   cfg$dataset$pathway2gene_file     (e.g. data/dataset_0/processed/Pathway2Gene.csv)
#
# Writes into paths$step06  (results/pipelineB/<run_id>/06_consensus/):
#   final_up.csv
#   final_down.csv
#
# Returns (invisibly):
#   final_up_file    path to final_up.csv
#   final_down_file  path to final_down.csv

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

step06_consensus_and_themes_A1A2 <- function(cfg, paths, ctx) {

  # ---- Parameters -----------------------------------------------------------
  REPS_UP       <- file.path(paths$step05, "representatives_up.csv")
  REPS_DOWN     <- file.path(paths$step05, "representatives_down.csv")
  UNIVERSE_FILE <- file.path(paths$step01, "universe_genes.txt")
  WP_MAP_FILE   <- cfg$dataset$pathway2gene_file
  OUT_DIR       <- paths$step06

  # All thresholds from cfg — stop with clear error if any are missing
  if (is.null(cfg$pipelineB$min_set_size))      stop("Missing config key: pipelineB$min_set_size")
  if (is.null(cfg$pipelineB$max_set_size))      stop("Missing config key: pipelineB$max_set_size")
  if (is.null(cfg$thresholds$CONS_JACCARD_MIN)) stop("Missing config key: thresholds$CONS_JACCARD_MIN")
  if (is.null(cfg$thresholds$CONS_K_MIN))       stop("Missing config key: thresholds$CONS_K_MIN")
  if (is.null(cfg$pipelineB$solo_q_max))        stop("Missing config key: pipelineB$solo_q_max")
  if (is.null(cfg$pipelineB$solo_abs_nes_min))  stop("Missing config key: pipelineB$solo_abs_nes_min")

  MIN_SET_SIZE     <- as.integer(cfg$pipelineB$min_set_size)
  MAX_SET_SIZE     <- as.integer(cfg$pipelineB$max_set_size)
  CONS_JACCARD_MIN <- as.numeric(cfg$thresholds$CONS_JACCARD_MIN)
  CONS_K_MIN       <- as.integer(cfg$thresholds$CONS_K_MIN)
  SOLO_Q_MAX       <- as.numeric(cfg$pipelineB$solo_q_max)
  SOLO_ABS_NES_MIN <- as.numeric(cfg$pipelineB$solo_abs_nes_min)

  stopifnot(
    file.exists(REPS_UP), file.exists(REPS_DOWN),
    file.exists(UNIVERSE_FILE), file.exists(WP_MAP_FILE)
  )

  LAST_RESORT_N <- if (!is.null(cfg$pipelineB$final_top_n_per_direction))
    as.integer(cfg$pipelineB$final_top_n_per_direction) else 5L

  # ---- Option A: fallback-mode candidate expansion --------------------------
  OPTION_A_ENABLE         <- if (!is.null(cfg$pipelineB$optionA_enable))
    isTRUE(cfg$pipelineB$optionA_enable) else TRUE
  OPTION_A_Q_FLOOR        <- if (!is.null(cfg$pipelineB$optionA_quality_floor_q))
    as.numeric(cfg$pipelineB$optionA_quality_floor_q) else 0.20
  OPTION_A_ALLOW_LOW_CONF <- if (!is.null(cfg$pipelineB$optionA_allow_low_conf))
    isTRUE(cfg$pipelineB$optionA_allow_low_conf) else FALSE

  # Detect whether Step03 used fallback (read agreement_source from step03 CSVs)
  detect_agree_src <- function(fpath) {
    if (!file.exists(fpath)) return("strict_intersection")
    tmp <- tryCatch(readr::read_csv(fpath, show_col_types = FALSE),
                   error = function(e) tibble::tibble())
    if ("agreement_source" %in% names(tmp) && nrow(tmp) > 0)
      as.character(tmp$agreement_source[1])
    else "strict_intersection"
  }
  agree_src_wp <- detect_agree_src(file.path(paths$step03, "wp_input_agreement.csv"))
  agree_src_go <- detect_agree_src(file.path(paths$step03, "go_bp_input_agreement.csv"))
  is_fallback  <- any(c(agree_src_wp, agree_src_go) %in% c("relaxed_fdr", "topn"))

  cat(sprintf("[Step06] Agreement sources \u2014 WP: %s, GO: %s | Option A active: %s\n",
              agree_src_wp, agree_src_go,
              if (is_fallback && OPTION_A_ENABLE) "YES" else "NO"))

  # ---- Helpers --------------------------------------------------------------
  read_gene_vec <- function(path) {
    x <- readr::read_lines(path)
    x <- stringr::str_trim(x)
    unique(x[x != ""])
  }

  normalize_direction <- function(x) {
    x2 <- toupper(trimws(as.character(x)))
    dplyr::case_when(
      x2 %in% c("UP", "UPREGULATED", "UPREGULATION")       ~ "UP",
      x2 %in% c("DOWN", "DOWNREGULATED", "DOWNREGULATION") ~ "DOWN",
      TRUE ~ x2
    )
  }

  getcol_ci <- function(df, candidates) {
    nm <- names(df); nml <- tolower(nm)
    hit <- nm[match(tolower(candidates), nml, nomatch = 0)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }

  get_wp_sets <- function(wp_ids, universe) {
    wp_map_raw <- readr::read_csv(WP_MAP_FILE, show_col_types = FALSE)
    id_col  <- getcol_ci(wp_map_raw, c("wpid","id","pathway_id","pathway","wp_id","wp"))
    sym_col <- getcol_ci(wp_map_raw, c("symbol","gene","gene_symbol","genesymbol","hgnc_symbol","hugo"))

    if (is.na(id_col) || is.na(sym_col)) {
      stop("Pathway2Gene.csv missing required columns. Found: ", paste(names(wp_map_raw), collapse = ", "))
    }

    wp_map <- wp_map_raw |>
      dplyr::transmute(
        ID     = as.character(.data[[id_col]]),
        SYMBOL = as.character(.data[[sym_col]])
      ) |>
      dplyr::filter(ID %in% wp_ids, SYMBOL %in% universe) |>
      dplyr::distinct()

    sets_tbl <- wp_map |>
      dplyr::group_by(ID) |>
      dplyr::summarise(genes = list(unique(SYMBOL)), .groups = "drop") |>
      dplyr::mutate(size = lengths(genes)) |>
      dplyr::filter(size >= MIN_SET_SIZE, size <= MAX_SET_SIZE)

    sets <- sets_tbl$genes
    names(sets) <- sets_tbl$ID
    sets
  }

  get_go_sets <- function(go_ids, universe) {
    go_map_raw <- suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys    = unique(go_ids),
        keytype = "GO",
        columns = c("SYMBOL", "ONTOLOGY", "GO")
      )
    )

    if (!"GO" %in% names(go_map_raw)) {
      stop("org.Hs.eg.db::select() returned no 'GO' column. Got: ", paste(names(go_map_raw), collapse = ", "))
    }
    if (!"ONTOLOGY" %in% names(go_map_raw)) go_map_raw$ONTOLOGY <- NA_character_

    go_map <- tibble::as_tibble(go_map_raw) |>
      dplyr::rename(ID = GO) |>
      dplyr::mutate(
        ID       = as.character(ID),
        SYMBOL   = as.character(SYMBOL),
        ONTOLOGY = as.character(ONTOLOGY)
      ) |>
      dplyr::filter(!is.na(ID), !is.na(SYMBOL)) |>
      dplyr::filter(is.na(ONTOLOGY) | ONTOLOGY == "BP") |>
      dplyr::filter(SYMBOL %in% universe) |>
      dplyr::distinct(ID, SYMBOL)

    sets_tbl <- go_map |>
      dplyr::group_by(ID) |>
      dplyr::summarise(genes = list(unique(SYMBOL)), .groups = "drop") |>
      dplyr::mutate(size = lengths(genes)) |>
      dplyr::filter(size >= MIN_SET_SIZE, size <= MAX_SET_SIZE)

    sets <- sets_tbl$genes
    names(sets) <- sets_tbl$ID
    sets
  }

  pairwise_go_wp <- function(go_sets, wp_sets, j_min, k_min) {
    empty <- tibble::tibble(GO = character(), WP = character(), jaccard = numeric(), k = integer())
    if (!length(go_sets) || !length(wp_sets)) return(empty)
    go_ids <- names(go_sets); wp_ids <- names(wp_sets)
    result <- purrr::map_dfr(go_ids, function(g) {
      purrr::map_dfr(wp_ids, function(w) {
        inter <- length(intersect(go_sets[[g]], wp_sets[[w]]))
        if (inter < k_min) return(NULL)
        uni <- length(union(go_sets[[g]], wp_sets[[w]]))
        if (uni == 0) return(NULL)
        j <- inter / uni
        if (j < j_min) return(NULL)
        tibble::tibble(GO = g, WP = w, jaccard = j, k = inter)
      })
    })
    if (nrow(result) == 0) empty else result
  }

  prep_reps <- function(df) {
    out <- df |>
      dplyr::mutate(
        direction   = normalize_direction(direction),
        collection  = as.character(collection),
        ID          = as.character(ID),
        Description = as.character(Description),
        set_size    = suppressWarnings(as.integer(set_size)),
        camera_FDR  = suppressWarnings(as.numeric(camera_FDR)),
        fgsea_padj  = suppressWarnings(as.numeric(fgsea_padj)),
        fgsea_NES   = suppressWarnings(as.numeric(fgsea_NES))
      ) |>
      dplyr::mutate(
        agreement_q = pmax(camera_FDR, fgsea_padj, na.rm = TRUE),
        absNES      = ifelse(!is.na(fgsea_NES), abs(fgsea_NES), NA_real_)
      )
    if ("cluster" %in% names(out)) out$cluster <- as.character(out$cluster)
    out
  }

  process_one_direction <- function(reps_all, dir_label, universe) {
    reps_dir <- reps_all |>
      dplyr::filter(direction == dir_label) |>
      dplyr::filter(set_size >= MIN_SET_SIZE, set_size <= MAX_SET_SIZE)

    if (nrow(reps_dir) == 0) {
      return(list(
        reps_in = 0L,
        final   = reps_dir |>
          dplyr::mutate(
            in_consensus   = logical(0),
            solo_high_conf = logical(0),
            keep_final     = logical(0),
            final_reason   = character(0)
          )
      ))
    }

    go_reps <- reps_dir |> dplyr::filter(collection == "GO")
    wp_reps <- reps_dir |> dplyr::filter(collection == "WP")

    go_sets <- get_go_sets(go_reps$ID, universe)
    wp_sets <- get_wp_sets(wp_reps$ID, universe)

    pairs <- pairwise_go_wp(go_sets, wp_sets, j_min = CONS_JACCARD_MIN, k_min = CONS_K_MIN)

    go_cons <- unique(pairs$GO)
    wp_cons <- unique(pairs$WP)

    final <- reps_dir |>
      dplyr::mutate(
        in_consensus = dplyr::case_when(
          collection == "GO" ~ ID %in% go_cons,
          collection == "WP" ~ ID %in% wp_cons,
          TRUE ~ FALSE
        ),
        solo_high_conf = (!is.na(agreement_q) & agreement_q <= SOLO_Q_MAX) &
                         (SOLO_ABS_NES_MIN < 0 |             # negative = no NES requirement (test sentinel)
                          (!is.na(absNES) & absNES >= SOLO_ABS_NES_MIN)),
        final_reason = dplyr::case_when(
          in_consensus                                                            ~ "consensus",
          solo_high_conf                                                          ~ "solo_high_conf",
          is_fallback & OPTION_A_ENABLE &
            !is.na(agreement_q) & agreement_q <= OPTION_A_Q_FLOOR               ~ "fallback_ranked",
          is_fallback & OPTION_A_ENABLE & OPTION_A_ALLOW_LOW_CONF               ~ "fallback_ranked_low_conf",
          TRUE                                                                    ~ NA_character_
        ),
        keep_final = !is.na(final_reason)
      )

    # ---- Last-resort fallback: no pathways kept but reps exist + Step03 used fallback
    if (sum(final$keep_final, na.rm = TRUE) == 0 && nrow(final) > 0 && is_fallback) {
      cat(sprintf(
        "Step06: no pathways passed strict rules; applying last-resort fallback_topn_last_resort (n = %d) for %s.\n",
        LAST_RESORT_N, dir_label
      ))
      top_ids <- final |>
        dplyr::arrange(agreement_q, dplyr::desc(absNES), set_size) |>
        dplyr::slice_head(n = LAST_RESORT_N) |>
        dplyr::pull(ID)
      final <- final |>
        dplyr::mutate(
          final_reason = dplyr::if_else(ID %in% top_ids, "fallback_topn_last_resort", final_reason),
          keep_final   = !is.na(final_reason)
        )
    }

    list(reps_in = nrow(reps_dir), final = final)
  }

  # ---- Run ------------------------------------------------------------------
  universe  <- read_gene_vec(UNIVERSE_FILE)

  reps_up   <- readr::read_csv(REPS_UP,   show_col_types = FALSE) |> prep_reps()
  reps_down <- readr::read_csv(REPS_DOWN, show_col_types = FALSE) |> prep_reps()

  need_cols    <- c("ID","Description","direction","collection","set_size","camera_FDR","fgsea_padj","agreement_q","absNES")
  missing_up   <- setdiff(need_cols, names(reps_up))
  missing_down <- setdiff(need_cols, names(reps_down))
  if (length(missing_up))   stop("Missing columns in UP reps: ",   paste(missing_up,   collapse = ", "))
  if (length(missing_down)) stop("Missing columns in DOWN reps: ", paste(missing_down, collapse = ", "))

  cat(glue::glue(
    "[Pipeline B \u2014 Step 06: Consensus + Themes]\n",
    "Inputs:\n",
    "- {REPS_UP}\n",
    "- {REPS_DOWN}\n",
    "Universe:\n",
    "- {UNIVERSE_FILE}\n",
    "WP map:\n",
    "- {WP_MAP_FILE}\n\n",
    "Thresholds (from cfg):\n",
    "- consensus Jaccard >= {CONS_JACCARD_MIN}  [thresholds$CONS_JACCARD_MIN]\n",
    "- consensus shared  >= {CONS_K_MIN}         [thresholds$CONS_K_MIN]\n",
    "- solo_q_max        <= {SOLO_Q_MAX}          [pipelineB$solo_q_max]\n",
    "- solo_abs_nes_min  >= {SOLO_ABS_NES_MIN}    [pipelineB$solo_abs_nes_min]\n"
  ), "\n")

  reps_all <- dplyr::bind_rows(reps_up, reps_down)

  message("Direction distribution in input:")
  print(reps_all |> dplyr::count(direction, sort = TRUE))

  up_res   <- process_one_direction(reps_all, "UP",   universe)
  down_res <- process_one_direction(reps_all, "DOWN", universe)

  UP_FILE   <- file.path(OUT_DIR, "final_up.csv")
  DOWN_FILE <- file.path(OUT_DIR, "final_down.csv")

  # Fill any missing GO descriptions before writing
  if (!exists("add_go_term_name", mode = "function"))
    source("R/utils/go_term_names.R")
  final_up_out   <- add_go_term_name(up_res$final,   id_col = "ID", desc_col = "Description")
  final_down_out <- add_go_term_name(down_res$final, id_col = "ID", desc_col = "Description")

  readr::write_csv(final_up_out,   UP_FILE)
  readr::write_csv(final_down_out, DOWN_FILE)

  cat(glue::glue(
    "\n[Step 06 complete]\n",
    "UP:   {up_res$reps_in} reps in -> {sum(up_res$final$keep_final, na.rm=TRUE)} kept\n",
    "DOWN: {down_res$reps_in} reps in -> {sum(down_res$final$keep_final, na.rm=TRUE)} kept\n",
    "Outputs in: {OUT_DIR}\n"
  ), "\n")

  invisible(list(
    final_up_file   = UP_FILE,
    final_down_file = DOWN_FILE
  ))
}
