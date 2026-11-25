# Marker names. Seurat turns "_" into "-" in feature names; the canonical names (CSV headers,
# QC_FLAGS, gating YAML flags) use "_". These two functions are the only place of the conversion.

# canonical -> Seurat feature name
to_seurat_name <- function(x) gsub("_", "-", x, fixed = TRUE)

# Seurat feature name -> canonical
to_canonical_name <- function(x) gsub("-", "_", x, fixed = TRUE)
