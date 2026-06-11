# shiny/test_dashboard_loaders.R
# Smoke test for the dashboard data loader functions.
# Run from the project root: source("shiny/test_dashboard_loaders.R")
# No external test framework required — uses stopifnot() for basic assertions.

message("=== MaRBLe Dashboard — data loader smoke test ===\n")

# Source the loader
dl_path <- if (file.exists("shiny/R/data_loader.R")) "shiny/R/data_loader.R" else
           if (file.exists("R/data_loader.R"))  "R/data_loader.R" else
           stop("Cannot find data_loader.R")
source(dl_path, local = FALSE)

ok <- function(label) message("  [OK] ", label)
fail <- function(label, err) message("  [FAIL] ", label, " — ", err)

# ── 1. detect_datasets ────────────────────────────────────────────────────────
tryCatch({
  ds <- detect_datasets()
  stopifnot(is.character(ds), length(ds) >= 1L)
  ok(paste0("detect_datasets() → ", paste(ds, collapse=", ")))
}, error = function(e) fail("detect_datasets()", conditionMessage(e)))

# ── 2. load_shortlist ─────────────────────────────────────────────────────────
for (cfg in list(
  list("dataset_0", "A", "Default",     "ds0 A Default"),
  list("dataset_0", "A", "Recommended", "ds0 A Recommended"),
  list("dataset_0", "B", "Default",     "ds0 B Default"),
  list("dataset_3", "A", "Default",     "ds3 A Default (22 terms)")
)) {
  tryCatch({
    r <- load_shortlist(cfg[[1]], cfg[[2]], cfg[[3]])
    stopifnot(is.list(r), "found" %in% names(r))
    if (r$found) {
      stopifnot(is.data.frame(r$data), nrow(r$data) > 0L)
      ok(sprintf("load_shortlist(%s) → found=TRUE, %d rows", cfg[[4]], nrow(r$data)))
    } else {
      ok(sprintf("load_shortlist(%s) → not found (expected if run missing)", cfg[[4]]))
    }
  }, error = function(e) fail(paste0("load_shortlist(", cfg[[4]], ")"), conditionMessage(e)))
}

# An unknown preset (e.g. "Relaxed") falls back to Default rather than erroring
tryCatch({
  r <- load_shortlist("dataset_0", "A", "Relaxed")
  stopifnot(isTRUE(r$found))
  ok("Unknown preset → falls back to Default (found=TRUE)")
}, error = function(e) fail("Unknown preset fallback", conditionMessage(e)))

# ── 3. load_summary_cards ─────────────────────────────────────────────────────
tryCatch({
  sc <- load_summary_cards("dataset_0", "A", "Default")
  stopifnot(is.list(sc), "dataset" %in% names(sc), "final_terms" %in% names(sc))
  ok(sprintf("load_summary_cards() → final_terms=%s, policy=%s", sc$final_terms, sc$policy_class))
}, error = function(e) fail("load_summary_cards()", conditionMessage(e)))

# ── 4. load_funnel_counts ─────────────────────────────────────────────────────
for (cfg in list(
  list("dataset_0", "A", "ds0 A"),
  list("dataset_0", "B", "ds0 B"),
  list("dataset_3", "A", "ds3 A")
)) {
  tryCatch({
    fd <- load_funnel_counts(cfg[[1]], cfg[[2]])
    stopifnot(is.data.frame(fd), nrow(fd) >= 1L)
    ok(sprintf("load_funnel_counts(%s) → %d stages: %s",
               cfg[[3]], nrow(fd), paste(fd$stage, collapse=" → ")))
  }, error = function(e) fail(paste0("load_funnel_counts(", cfg[[3]], ")"), conditionMessage(e)))
}

# ── 5. find_pathway2gene_file ─────────────────────────────────────────────────
tryCatch({
  p <- find_pathway2gene_file("dataset_0")
  stopifnot(!is.na(p), file.exists(p))
  ok(paste0("find_pathway2gene_file(dataset_0) → ", basename(p)))
}, error = function(e) fail("find_pathway2gene_file()", conditionMessage(e)))

message("\n=== Smoke test complete ===")
