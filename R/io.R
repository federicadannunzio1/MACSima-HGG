# Names of the intermediate RDS files, in one place. A step asks for a file by key, so a typo
# stops with an error instead of reading or writing a wrong file. Needs RDS_DIR (code/utils/constants.R).
RDS_FILES <- c(
  harmonization          = "00_harmonization.rds",
  canonical_to_seurat    = "01_canonical_to_seurat_names.rds",
  seurat_raw             = "01_seurat_raw.rds",
  seurat_qc              = "02_seurat_qc.rds",
  seurat_normalized      = "03_seurat_normalized.rds",
  seurat_dimred          = "04_seurat_dimred.rds",
  seurat_clustered       = "05_seurat_clustered.rds",
  seurat_annotated       = "05b_seurat_annotated.rds",
  seurat_gated           = "06_seurat_gated.rds",
  mfi_gated              = "06_mfi_gated.rds",
  patient_features       = "07_patient_features.rds"
)

# Full path of the RDS file with this key.
rds_path <- function(key) {
  if (length(key) != 1 || !key %in% names(RDS_FILES))
    stop("unknown RDS key: ", paste(key, collapse = ","), " (known: ", paste(names(RDS_FILES), collapse = ", "), ")")
  file.path(RDS_DIR, RDS_FILES[[key]])
}
