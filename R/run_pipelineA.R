# run_pipelineA.R
# Entry point for Pipeline A.

library(yaml)

source("R/utils/config_loader.R")
source("R/utils/run_context.R")
source("R/utils/paths_builder.R")
source("R/utils/validate_inputs.R")
source("R/pipelineA/steps/step01_make_universe.R")
source("R/pipelineA/steps/step02_make_de.R")
source("R/pipelineA/steps/step03a_ora_wp.R")
source("R/pipelineA/steps/step03b_ora_go_bp.R")
source("R/pipelineA/steps/step04a_go_collapse_semantic.R")
source("R/pipelineA/steps/step05_overlap_clustering.R")
source("R/pipelineA/steps/step06_bootstrap_consensus.R")
source("R/pipelineA/steps/step07_make_tiers.R")
source("R/pipelineA/steps/step98_export_final.R")

cfg   <- load_config("config/default.yml")
ctx   <- create_run_context("pipelineA", cfg)
paths <- build_paths(ctx)

message("Pipeline A run initialized at ", ctx$run_dir)
message("Dataset: ", cfg$dataset$name)

validate_config_and_inputs(cfg, "pipelineA")

# ---------------------------------------------------------------------------
# Step-toggle helpers
# ---------------------------------------------------------------------------

enabled_steps <- if (!is.null(cfg$pipelineA$run_steps)) {
  as.integer(cfg$pipelineA$run_steps)
} else {
  c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 98L)
}

step_enabled <- function(n) as.integer(n) %in% enabled_steps

message("Steps enabled: ", paste(enabled_steps, collapse = ", "))

# ---------------------------------------------------------------------------
# Manifest helper — writes manifest.yml into a step folder after it runs
# ---------------------------------------------------------------------------

write_manifest <- function(step_dir, step_num, inputs, outputs) {
  yaml::write_yaml(
    list(
      step      = step_num,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      inputs    = as.list(inputs),
      outputs   = as.list(outputs)
    ),
    file.path(step_dir, "manifest.yml")
  )
}

ran_steps         <- integer(0)
go_bp_ora_empty   <- FALSE
early_stop_reason <- NA_character_

# ---------------------------------------------------------------------------
# Step 01: Gene universe
# ---------------------------------------------------------------------------

if (step_enabled(1)) {
  result01 <- step01_make_universe(cfg, paths, ctx)
  message("Universe written to: ", result01$universe_file)
  write_manifest(paths$step01, 1,
    inputs  = list(expression_file   = cfg$dataset$expression_file,
                   pathway2gene_file = cfg$dataset$pathway2gene_file),
    outputs = list(universe_file   = result01$universe_file)
  )
  message("--- Pipeline A Step01 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 1L)
}

# ---------------------------------------------------------------------------
# Step 02: DE gene list
# ---------------------------------------------------------------------------

if (step_enabled(2)) {
  result02 <- step02_make_de(cfg, paths, ctx)
  message("DE genes written to:  ", result02$de_genes_file)
  message("DE policy written to: ", result02$de_policy_file)
  write_manifest(paths$step02, 2,
    inputs  = list(step01_dir      = paths$step01,
                   statistics_file = cfg$dataset$statistics_file),
    outputs = list(de_genes_file   = result02$de_genes_file,
                   de_policy_file  = result02$de_policy_file)
  )
  message("--- Pipeline A Step02 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 2L)
}

# ---------------------------------------------------------------------------
# Step 03: ORA (WikiPathways + GO:BP)
# ---------------------------------------------------------------------------

if (step_enabled(3)) {
  result03a <- step03a_ora_wp(cfg, paths, ctx)
  message("WP ORA filtered written to:    ", result03a$ora_filtered_file)

  result03b <- step03b_ora_go_bp(cfg, paths, ctx)
  message("GO-BP ORA filtered written to: ", result03b$ora_filtered_file)

  write_manifest(paths$step03, 3,
    inputs  = list(step01_dir       = paths$step01,
                   step02_dir       = paths$step02),
    outputs = list(wp_filtered_file = result03a$ora_filtered_file,
                   go_filtered_file = result03b$ora_filtered_file)
  )
  message("--- Pipeline A Step03 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 3L)
}

# ---------------------------------------------------------------------------
# Early-stop check: if GO-BP filtered table is empty, skip steps 04 onward
# ---------------------------------------------------------------------------

if (3L %in% ran_steps) {
  go_bp_filt_nrow <- tryCatch(
    nrow(readr::read_csv(result03b$ora_filtered_file, show_col_types = FALSE)),
    error = function(e) 0L
  )
  if (go_bp_filt_nrow == 0L) {
    go_bp_ora_empty   <- TRUE
    early_stop_reason <- "GO-BP ORA filtered table is empty (0 pathways after filters + fallback)"
    message("NOTE: ", early_stop_reason)
    for (s in intersect(enabled_steps, c(4L, 5L, 6L, 7L, 98L)))
      message("  Step", formatC(s, width = 2, flag = "0"),
              " skipped: empty GO-BP ORA filtered table")
    enabled_steps <- setdiff(enabled_steps, c(4L, 5L, 6L, 7L, 98L))
  }
}

# ---------------------------------------------------------------------------
# Step 04: GO semantic collapse
# ---------------------------------------------------------------------------

if (step_enabled(4)) {
  result04a <- step04a_go_collapse_semantic(cfg, paths, ctx)
  message("GO collapse reps written to: ", result04a$reps_file)
  write_manifest(paths$step04, 4,
    inputs  = list(step03_go_bp_dir = paths$step03_go_bp),
    outputs = list(reps_file        = result04a$reps_file,
                   mapping_file     = result04a$mapping_file)
  )
  message("--- Pipeline A Step04 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 4L)
}

# ---------------------------------------------------------------------------
# Step 05: Overlap clustering
# ---------------------------------------------------------------------------

if (step_enabled(5)) {
  result05 <- step05_overlap_clustering(cfg, paths, ctx)
  message("GO overlap reps written to:  ", result05$go_reps_file)
  message("WP overlap reps written to:  ", result05$wp_reps_file)
  write_manifest(paths$step05, 5,
    inputs  = list(step04_dir    = paths$step04,
                   step03_wp_dir = paths$step03_wp,
                   step01_dir    = paths$step01),
    outputs = list(go_reps_file  = result05$go_reps_file,
                   wp_reps_file  = result05$wp_reps_file)
  )
  message("--- Pipeline A Step05 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 5L)
}

# ---------------------------------------------------------------------------
# Step 06: Bootstrap stability + consensus
# ---------------------------------------------------------------------------

if (step_enabled(6)) {
  result06 <- step06_bootstrap_consensus(cfg, paths, ctx)
  message("GO stability table written to:    ", result06$go_stability_table)
  message("WP stability table written to:    ", result06$wp_stability_table)
  message("Consensus pairs written to:       ", result06$consensus_pairs_table)
  message("Final shortlist written to:       ", result06$final_shortlist_table)
  write_manifest(paths$step06, 6,
    inputs  = list(step01_dir = paths$step01,
                   step05_dir = paths$step05),
    outputs = list(go_stability_table   = result06$go_stability_table,
                   wp_stability_table   = result06$wp_stability_table,
                   consensus_pairs      = result06$consensus_pairs_table,
                   final_shortlist      = result06$final_shortlist_table)
  )
  message("--- Pipeline A Step06 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 6L)
}

# ---------------------------------------------------------------------------
# Step 07: Tiers
# ---------------------------------------------------------------------------

if (step_enabled(7)) {
  result07 <- step07_make_tiers(cfg, paths, ctx)
  message("Tier 1 written to:                ", result07$tier1_file)
  message("Tier 2 written to:                ", result07$tier2_file)
  message("All candidates (tiered) to:       ", result07$final_table_file)
  write_manifest(paths$step07, 7,
    inputs  = list(step06_dir       = paths$step06),
    outputs = list(tier1_file       = result07$tier1_file,
                   tier2_file       = result07$tier2_file,
                   final_table_file = result07$final_table_file)
  )
  message("--- Pipeline A Step07 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 7L)
}

# ---------------------------------------------------------------------------
# Step 98: Export canonical final list to FINAL/
# ---------------------------------------------------------------------------

if (step_enabled(98)) {
  result98 <- step98_export_final(cfg, paths, ctx)
  message("Final pathways exported to: ", result98$final_file)
  message("Final summary written to:   ", result98$summary_file)
  message("--- Pipeline A Step98 complete [dataset: ", cfg$dataset$name, "] ---")
  ran_steps <- c(ran_steps, 98L)
}

# ---------------------------------------------------------------------------
# Run metadata
# ---------------------------------------------------------------------------

git_hash <- tryCatch(
  trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE)),
  error = function(e) NA_character_
)

run_meta <- list(
  pipeline          = "pipelineA",
  run_id            = ctx$run_id,
  timestamp         = ctx$timestamp,
  dataset           = cfg$dataset$name,
  ran_steps         = as.list(ran_steps),
  git_commit        = if (length(git_hash) == 1 && nchar(git_hash) > 0) git_hash else NA_character_,
  run_dir           = ctx$run_dir,
  early_stop_reason = early_stop_reason
)

yaml::write_yaml(run_meta, file.path(ctx$run_dir, "run_meta.yml"))

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------

message("")
message("=== Pipeline A run complete ===")
message("Run dir:   ", ctx$run_dir)
message("Dataset:   ", cfg$dataset$name)
message("Steps run: ",
        if (length(ran_steps) > 0) paste(sort(ran_steps), collapse = ", ") else "(none)")
if (!is.na(early_stop_reason))
  message("Early stop: ", early_stop_reason)
message("Metadata:  ", file.path(ctx$run_dir, "run_meta.yml"))
