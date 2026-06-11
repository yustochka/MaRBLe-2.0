# R/utils/go_term_names.R
# Helper: fill human-readable Description for GO term rows using GO.db.
#
# add_go_term_name(df, id_col = "ID", desc_col = "Description")
#   - Only fills rows where id starts with "GO:" AND desc is NA or empty.
#   - For WP rows (id starts with "WP") the column is left unchanged.
#   - Uses GO.db::GOTERM[[go_id]] and AnnotationDbi::Term() for direct lookup.
#   - If GO.db is unavailable, prints a warning and leaves NAs unchanged.
#   - Never removes or reorders rows.
#   - Prints one log line: "GO term name lookup: filled X / Y GO descriptions."

add_go_term_name <- function(df, id_col = "ID", desc_col = "Description") {
  if (nrow(df) == 0) return(df)

  # Ensure desc_col exists
  if (!desc_col %in% names(df)) df[[desc_col]] <- NA_character_

  ids     <- as.character(df[[id_col]])
  current <- as.character(df[[desc_col]])
  needs   <- grepl("^GO:", ids) & (is.na(current) | trimws(current) == "")

  if (!any(needs)) return(df)

  go_ids <- unique(ids[needs])

  term_map <- tryCatch({
    suppressPackageStartupMessages(library(GO.db))
    # Look up each GO ID via GOTERM environment; Term() extracts the name string
    terms <- vapply(go_ids, function(gid) {
      obj <- tryCatch(GO.db::GOTERM[[gid]], error = function(e) NULL)
      if (!is.null(obj)) as.character(AnnotationDbi::Term(obj)) else NA_character_
    }, character(1L))
    stats::setNames(terms, go_ids)
  }, error = function(e) {
    message("go_term_names: GO.db lookup failed (", conditionMessage(e),
            "). GO descriptions will remain NA.")
    NULL
  })

  if (!is.null(term_map)) {
    filled     <- term_map[ids]
    df[[desc_col]] <- dplyr::if_else(needs & !is.na(filled), filled, current)
    n_filled   <- sum(needs & !is.na(term_map[ids]), na.rm = TRUE)
    cat(sprintf("GO term name lookup: filled %d / %d GO descriptions.\n",
                n_filled, sum(needs)))
  }

  df
}
