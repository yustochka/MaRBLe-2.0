# R/ablation/utils/ablation_config.R
# Merge base config (config/default.yml) with an ablation variant YAML.
#
# The variant YAML only sets pipelineB.ablation — all dataset settings,
# thresholds, and other pipelineB parameters come from the base config.
# This guarantees that only one variable changes per ablation run.

merge_ablation_config <- function(base_cfg, ablation_yml_path) {
  if (!file.exists(ablation_yml_path))
    stop("Ablation variant config not found: ", ablation_yml_path)

  override <- yaml::read_yaml(ablation_yml_path)

  if (is.null(override$pipelineB) || is.null(override$pipelineB$ablation))
    stop("Ablation YAML must contain a pipelineB.ablation block.\n",
         "File: ", ablation_yml_path)

  abl <- override$pipelineB$ablation

  # Validate required keys
  required_keys <- c("variant_name", "use_fgsea", "use_semantic_collapse",
                     "use_overlap_clustering", "use_consensus")
  missing_keys  <- setdiff(required_keys, names(abl))
  if (length(missing_keys) > 0)
    stop("Ablation YAML missing required keys: ",
         paste(missing_keys, collapse = ", "), "\nFile: ", ablation_yml_path)

  # Graft ablation block onto base config (all other settings unchanged)
  base_cfg$pipelineB$ablation <- abl

  # Allow test variants to override other pipelineB keys (e.g. solo_abs_nes_min).
  # Only keys explicitly present in the override YAML (besides 'ablation') are applied.
  extra_keys <- setdiff(names(override$pipelineB), "ablation")
  for (key in extra_keys)
    base_cfg$pipelineB[[key]] <- override$pipelineB[[key]]

  base_cfg
}
