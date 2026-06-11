# R/utils/paths_builder.R
# Build and create per-step output directories inside a run folder.
# build_paths()   -> Pipeline A directory map
# build_paths_b() -> Pipeline B directory map

build_paths <- function(ctx) {
  make_dir <- function(name) {
    path <- file.path(ctx$run_dir, name)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
  }

  list(
    step01       = make_dir("01_universe"),
    step02       = make_dir("02_de"),
    step03       = make_dir("03_ora"),
    step03_wp    = make_dir("03_ora/wp"),
    step03_go_bp = make_dir("03_ora/go_bp"),
    step04       = make_dir("04_go_collapse"),
    step05 = make_dir("05_overlap"),
    step06 = make_dir("06_bootstrap_consensus"),
    step07 = make_dir("07_tiers")
  )
}

build_paths_b <- function(ctx) {
  make_dir <- function(name) {
    path <- file.path(ctx$run_dir, name)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
  }

  list(
    step01 = make_dir("01_camera"),
    step02 = make_dir("02_fgsea"),
    step03 = make_dir("03_prepare_inputs"),
    step04 = make_dir("04_go_collapse"),
    step05 = make_dir("05_overlap"),
    step06 = make_dir("06_consensus")
  )
}