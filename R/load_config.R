# ------------------------------------------------------------
# load_config.R
# Utilities for reading and merging YAML configuration files
# ------------------------------------------------------------

#' Load and merge common + model-specific configuration
#'
#' @param model Character string identifying the model config
#'              (e.g. "tree_monomolecular")
#' @return Named list of parameters
#'
load_config <- function(model) {
  
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required but not installed.")
  }
  
  common_path <- "config/common.yml"
  model_path  <- file.path("config", paste0(model, ".yml"))
  
  if (!file.exists(common_path)) {
    stop("Missing config file: ", common_path)
  }
  if (!file.exists(model_path)) {
    stop("Missing config file: ", model_path)
  }
  
  common_cfg <- yaml::read_yaml(common_path)
  model_cfg  <- yaml::read_yaml(model_path)
  
  # Model-specific values override common ones
  cfg <- modifyList(common_cfg, model_cfg)
  
  return(cfg)
}