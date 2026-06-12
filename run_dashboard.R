#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# One-command launcher for the MaRBLe dashboard.
#
#   Rscript run_dashboard.R
#
# Installs the three viewer packages (shiny, DT, ggplot2) if they are missing,
# then starts the app in your browser. This does NOT need renv or the
# Bioconductor pipeline stack — those are only required to re-run the
# pipelines themselves (see README).
# ---------------------------------------------------------------------------

pkgs    <- c("shiny", "DT", "ggplot2")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  message("Installing dashboard packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

shiny::runApp("shiny", launch.browser = TRUE)
