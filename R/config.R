# Reads config/pipeline.yaml and checks that the keys the code needs are there.
read_pipeline_config <- function(path = file.path(getwd(), "config", "pipeline.yaml")) {
  cfg <- yaml::read_yaml(path)
  need_sig <- c("exhaustion", "m1", "m2", "hla", "checkpoint_ligand", "ki67", "tumor_antigens")
  stopifnot(is.numeric(cfg$pfs_cutoff_months), length(cfg$pfs_cutoff_months) == 1,
            all(need_sig %in% names(cfg$signatures)),
            all(vapply(cfg$signatures[need_sig], function(x) is.character(x) && length(x) > 0, logical(1))))
  cfg$signatures <- lapply(cfg$signatures, as.character)
  # params of steps 02-06: numbers as double (as when they were written in R), resolutions as a vector
  need_par <- c("cofactor_quantile", "min_total_mfi_pct", "n_pcs", "harmony_group", "umap_neighbors",
                "umap_mindist", "louvain_resolutions", "point_size", "dpi")
  miss <- setdiff(need_par, names(cfg$params))
  if (length(miss) > 0) stop("config/pipeline.yaml: params missing: ", paste(miss, collapse = ", "))
  cfg$params <- lapply(cfg$params[need_par], function(x) if (is.numeric(unlist(x))) as.numeric(unlist(x)) else as.character(x))
  stopifnot(length(cfg$params$harmony_group) == 1, all(lengths(cfg$params[setdiff(need_par, "louvain_resolutions")]) == 1))
  cfg
}
