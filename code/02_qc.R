# =============================================================================
# 02_qc.R
# Cell-level QC: remove debris, compute QC metrics, produce QC figures.
# Output: Seurat object post-QC + QC plots.
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))
source(file.path(getwd(), "code/utils/plot_utils.R"))

library(Seurat)
library(readr)
library(ggplot2)
library(dplyr)
library(patchwork)

message("=== Step 2: Cell QC ===")

seurat_obj <- readRDS(file.path(RDS_DIR, "01_seurat_raw.rds"))

# ── 1. Compute QC metrics ─────────────────────────────────────────────────────
mfi_mat <- GetAssayData(seurat_obj, assay = "MFI", layer = "counts")
seurat_obj$total_MFI  <- colSums(mfi_mat)
seurat_obj$n_pos_markers <- colSums(mfi_mat > 0)

# ── 2. Debris filter: remove cells below 1st percentile of total MFI ─────────
threshold_total_mfi <- quantile(seurat_obj$total_MFI, PARAMS$min_total_mfi_pct)
message(sprintf("Total MFI debris threshold (p%.0f): %.2f",
                PARAMS$min_total_mfi_pct * 100, threshold_total_mfi))

keep <- seurat_obj$total_MFI >= threshold_total_mfi
message(sprintf("Cells before QC: %d", ncol(seurat_obj)))
message(sprintf("Cells removed (debris): %d", sum(!keep)))
message(sprintf("Cells after QC: %d", sum(keep)))

seurat_obj <- seurat_obj[, keep]

# ── 3. Per-ROI cell count after QC ──────────────────────────────────────────
counts_post_qc <- seurat_obj@meta.data %>%
  group_by(patient_id, sample_id, roi_id, pfs_group) %>%
  summarise(n_cells_postQC = n(), .groups = "drop")
message("\n--- Cell counts post-QC ---")
print(as.data.frame(counts_post_qc))
write_tsv(counts_post_qc, file.path(TAB_DIR, "02_cell_counts_postQC.tsv"))

# ── 4. QC plots ───────────────────────────────────────────────────────────────
p1 <- qc_violin(seurat_obj, "total_MFI", group_by = "sample_id",
                y_lab = "Total MFI (sum over all markers)") +
  geom_hline(yintercept = threshold_total_mfi, linetype = "dashed", color = "red") +
  ggtitle("Total MFI per cell (post-QC)")

p2 <- qc_violin(seurat_obj, "n_pos_markers", group_by = "sample_id",
                y_lab = "N markers with MFI > 0", fill_col = "#4DAF4A") +
  ggtitle("Marker positivity per cell")

df_meta <- seurat_obj@meta.data
p3 <- ggplot(df_meta, aes(x = total_MFI, color = patient_id)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(values = PATIENT_COLORS) +
  scale_x_log10() +
  labs(x = "Total MFI (log10)", y = "Density", title = "MFI distribution per patient",
       color = "Patient") +
  theme_hgg()

# Per-patient cell counts colored by PFS group
counts_patient <- counts_post_qc %>%
  group_by(patient_id, pfs_group) %>%
  summarise(n_cells_postQC = sum(n_cells_postQC), .groups = "drop")

p4 <- ggplot(counts_patient, aes(x = reorder(patient_id, -n_cells_postQC),
                                  y = n_cells_postQC, fill = pfs_group)) +
  geom_col() +
  geom_text(aes(label = n_cells_postQC), vjust = -0.3, size = 3) +
  scale_fill_manual(values = PFS_GROUP_COLORS) +
  labs(x = NULL, y = "N cells", title = "Cells per patient (post-QC)", fill = "PFS group") +
  theme_hgg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

qc_panel <- (p1 + p2) / (p3 + p4)
save_plot(qc_panel,
          file.path(FIG_DIR, "qc", "02_qc_panel.pdf"),
          width = 16, height = 10)

# ── 5. Save ──────────────────────────────────────────────────────────────────
saveRDS(seurat_obj, file.path(RDS_DIR, "02_seurat_qc.rds"))
message("\nStep 2 complete. Saved: output/data/02_seurat_qc.rds")
