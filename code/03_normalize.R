# =============================================================================
# 03_normalize.R
# Arcsinh normalization with empirical cofactor.
# cofactor[m] = median( quantile(MFI[m, cells of patient p], P50) ) over all p
# MFI_norm = arcsinh(MFI / cofactor)
# Output: Seurat with normalized values in "data" layer.
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))
source(file.path(getwd(), "code/utils/plot_utils.R"))

library(Seurat)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

message("=== Step 3: Arcsinh normalization ===")

seurat_obj <- readRDS(file.path(RDS_DIR, "02_seurat_qc.rds"))

mfi_raw    <- GetAssayData(seurat_obj, assay = "MFI", layer = "counts")
markers    <- rownames(mfi_raw)
patients   <- unique(seurat_obj$patient_id)
cell_meta  <- seurat_obj@meta.data

# ── 1. Compute per-patient P50 per marker ─────────────────────────────────────
p50_matrix <- matrix(NA, nrow = length(markers), ncol = length(patients),
                     dimnames = list(markers, patients))

for (pat in patients) {
  cells_pat <- rownames(cell_meta)[cell_meta$patient_id == pat]
  mfi_pat   <- mfi_raw[, cells_pat, drop = FALSE]
  p50_matrix[, pat] <- apply(mfi_pat, 1, quantile,
                              probs = PARAMS$cofactor_quantile, na.rm = TRUE)
}

# ── 2. Global cofactor = median of per-patient P50s ───────────────────────────
cofactors <- apply(p50_matrix, 1, median, na.rm = TRUE)

zero_cof <- names(cofactors)[cofactors == 0]
if (length(zero_cof) > 0) {
  message("WARNING: cofactor = 0 for markers (setting to 1): ",
          paste(zero_cof, collapse = ", "))
  cofactors[zero_cof] <- 1
}

message(sprintf("Cofactors computed for %d markers", length(cofactors)))
message("Cofactor range: ", round(min(cofactors), 3), " - ", round(max(cofactors), 3))

# ── 3. Apply arcsinh transformation ───────────────────────────────────────────
cofactor_mat <- matrix(cofactors, nrow = length(markers), ncol = ncol(mfi_raw),
                       byrow = FALSE, dimnames = dimnames(mfi_raw))
mfi_norm <- asinh(mfi_raw / cofactor_mat)

# ── 4. Store normalized values in Seurat "data" layer ─────────────────────────
seurat_obj <- SetAssayData(seurat_obj, assay = "MFI", layer = "data",
                            new.data = mfi_norm)

# ── 5. Save cofactor table ────────────────────────────────────────────────────
cof_df <- data.frame(
  marker    = names(cofactors),
  cofactor  = cofactors,
  stringsAsFactors = FALSE
) %>% arrange(desc(cofactor))

write_tsv(cof_df, file.path(TAB_DIR, "03_arcsinh_cofactors.tsv"))

# ── 6. Diagnostic plots ───────────────────────────────────────────────────────
p_cof <- cof_df %>% head(30) %>%
  ggplot(aes(x = reorder(marker, cofactor), y = cofactor)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = NULL, y = "Arcsinh cofactor (P50 median)", title = "Top 30 markers by cofactor") +
  theme_hgg()

key_markers <- c("CD3", "CD45", "GFAP", "CD68", "CD274", "H3.3K27M")
key_markers <- intersect(key_markers, markers)

plot_before_after <- function(marker) {
  raw  <- mfi_raw[marker, ]
  norm <- mfi_norm[marker, ]
  pat  <- cell_meta$patient_id

  df <- data.frame(
    raw  = as.numeric(raw),
    norm = as.numeric(norm),
    patient = pat
  ) %>% pivot_longer(c(raw, norm), names_to = "type", values_to = "value") %>%
    mutate(type = factor(type, levels = c("raw", "norm"),
                         labels = c("Raw MFI", "arcsinh(MFI/cofactor)")))

  ggplot(df, aes(x = value, color = patient)) +
    geom_density(linewidth = 0.7) +
    facet_wrap(~type, scales = "free_x") +
    scale_color_manual(values = PATIENT_COLORS) +
    labs(title = marker, x = NULL, y = "Density", color = "Patient") +
    theme_hgg() +
    theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"))
}

norm_plots <- lapply(key_markers, plot_before_after)
norm_panel <- wrap_plots(norm_plots, ncol = 2)

save_plot(p_cof,      file.path(FIG_DIR, "qc", "03_cofactors.pdf"),           width = 8,  height = 7)
save_plot(norm_panel,  file.path(FIG_DIR, "qc", "03_normalization_check.pdf"), width = 14, height = 3.5 * length(key_markers))

# ── 7. ScaleData (centering for PCA) ─────────────────────────────────────────
seurat_obj <- ScaleData(seurat_obj, assay = "MFI", features = markers, verbose = FALSE)
message("ScaleData complete.")

# ── 8. Save ──────────────────────────────────────────────────────────────────
saveRDS(seurat_obj, file.path(RDS_DIR, "03_seurat_normalized.rds"))
message("\nStep 3 complete. Saved: output/data/03_seurat_normalized.rds")
