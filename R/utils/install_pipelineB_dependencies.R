# R/utils/install_pipelineB_dependencies.R
# Check and install missing dependencies for Pipeline B.
# Run once before executing any Pipeline B steps:
#   Rscript R/utils/install_pipelineB_dependencies.R

BIOC_PKGS <- c(
  "limma",
  "fgsea",
  "edgeR",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "BiocParallel"
)

already_present <- c()
to_install      <- c()

for (pkg in BIOC_PKGS) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    already_present <- c(already_present, pkg)
  } else {
    to_install <- c(to_install, pkg)
  }
}

cat("Already installed:", paste(already_present, collapse = ", "), "\n")
cat("Missing:          ", paste(if (length(to_install)) to_install else "(none)", collapse = ", "), "\n\n")

if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  cat("Installing via BiocManager:", paste(to_install, collapse = ", "), "\n")
  BiocManager::install(to_install, ask = FALSE, update = FALSE)
  cat("Done.\n")
} else {
  cat("All Pipeline B dependencies are present. Nothing to install.\n")
}
