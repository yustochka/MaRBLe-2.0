# R/utils/config_loader.R
# Load pipeline configuration from a YAML file.

load_config <- function(path = "config/default.yml") {
  if (!file.exists(path)) {
    stop("Config file not found: ", path)
  }
  yaml::read_yaml(path)
}