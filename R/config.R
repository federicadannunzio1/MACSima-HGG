# Reads config/pipeline.yaml and checks that the keys the code needs are there.
read_pipeline_config <- function(path = file.path(getwd(), "config", "pipeline.yaml")) {
  cfg <- yaml::read_yaml(path)
  need_sig <- c("exhaustion", "m1", "m2", "hla", "checkpoint_ligand", "ki67", "tumor_antigens")
  stopifnot(is.numeric(cfg$pfs_cutoff_months), length(cfg$pfs_cutoff_months) == 1,
            all(need_sig %in% names(cfg$signatures)),
            all(vapply(cfg$signatures[need_sig], function(x) is.character(x) && length(x) > 0, logical(1))))
  cfg$signatures <- lapply(cfg$signatures, as.character)
  cfg
}
