# R/utils/run_context.R
# Create a timestamped run folder, save a config snapshot, and set the RNG seed.

create_run_context <- function(pipeline_name, cfg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  run_dir   <- file.path("results", pipeline_name, timestamp)

  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  yaml::write_yaml(cfg, file.path(run_dir, "config_used.yml"))

  seed <- cfg$run$seed
  set.seed(seed)

  list(
    run_id    = paste0(pipeline_name, "_", timestamp),
    run_dir   = run_dir,
    timestamp = timestamp,
    seed      = seed
  )
}