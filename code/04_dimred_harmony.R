# =============================================================================
# 04_dimred_harmony.R
# PCA → Harmony batch correction (by patient) → UMAP.
# Output: Seurat with "pca", "harmony", "umap" reductions.
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))
source(file.path(getwd(), "code/utils/plot_utils.R"))

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

message("=== Step 4: PCA + Harmony + UMAP ===")

seurat_obj <- readRDS(file.path(RDS_DIR, "03_seurat_normalized.rds"))
markers    <- rownames(seurat_obj)

# ── 1. PCA ───────────────────────────────────────────────────────────────────
message("Running PCA...")
seurat_obj <- RunPCA(
  seurat_obj,
  assay    = "MFI",
  features = markers,
  npcs     = PARAMS$n_pcs,
  reduction.name = "pca",
  verbose  = FALSE
)

p_elbow <- ElbowPlot(seurat_obj, ndims = PARAMS$n_pcs, reduction = "pca") +
  geom_vline(xintercept = 20, linetype = "dashed", color = "red") +
  labs(title = "PCA elbow plot") +
  theme_hgg()

p_pca_patient <- DimPlot(seurat_obj, reduction = "pca", group.by = "patient_id",
                          pt.size = PARAMS$point_size) +
  scale_color_manual(values = PATIENT_COLORS) +
  labs(title = "PCA — by patient (pre-Harmony)", color = "Patient") +
  theme_hgg()

# ── 2. Harmony batch correction ───────────────────────────────────────────────
message("Running Harmony (batch = patient_id)...")
seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars  = PARAMS$harmony_group,
  reduction      = "pca",
  dims.use       = 1:PARAMS$n_pcs,
  reduction.save = "harmony",
  verbose        = TRUE,
  assay.use      = "MFI"
)

p_harmony_patient <- DimPlot(seurat_obj, reduction = "harmony", group.by = "patient_id",
                              dims = c(1, 2), pt.size = PARAMS$point_size) +
  scale_color_manual(values = PATIENT_COLORS) +
  labs(title = "Harmony embedding — by patient", color = "Patient") +
  theme_hgg()

# ── 3. UMAP on Harmony embedding ─────────────────────────────────────────────
message("Running UMAP on Harmony embedding...")
seurat_obj <- RunUMAP(
  seurat_obj,
  reduction      = "harmony",
  dims           = 1:PARAMS$n_pcs,
  n.neighbors    = PARAMS$umap_neighbors,
  min.dist       = PARAMS$umap_mindist,
  reduction.name = "umap",
  verbose        = FALSE
)

# ── 4. UMAP plots ───────────────────────────────────────────────────────────
p_umap_patient <- DimPlot(seurat_obj, reduction = "umap", group.by = "patient_id",
                           pt.size = PARAMS$point_size) +
  scale_color_manual(values = PATIENT_COLORS) +
  labs(title = "UMAP — by patient", color = "Patient") +
  theme_hgg()

p_umap_pfs <- DimPlot(seurat_obj, reduction = "umap", group.by = "pfs_group",
                       pt.size = PARAMS$point_size) +
  scale_color_manual(values = PFS_GROUP_COLORS) +
  labs(title = "UMAP — by PFS group", color = "PFS group") +
  theme_hgg()

p_umap_sample <- DimPlot(seurat_obj, reduction = "umap", group.by = "sample_id",
                          pt.size = PARAMS$point_size) +
  labs(title = "UMAP — by ROI", color = "ROI") +
  theme_hgg() + NoLegend()

# ── 5. Feature plots ─────────────────────────────────────────────────────────
key_markers <- c("CD3", "CD45", "GFAP", "CD68", "CD56", "CD274", "FoxP3", "CD8a")
key_markers <- intersect(key_markers, markers)

feat_plots <- FeaturePlot(
  seurat_obj, features = key_markers, reduction = "umap",
  pt.size = PARAMS$point_size, ncol = 4, order = TRUE
) & theme_hgg() & theme(legend.position = "right")

# ── 6. Save plots ────────────────────────────────────────────────────────────
dimred_dir <- file.path(FIG_DIR, "dimred")

save_plot(p_elbow,           file.path(dimred_dir, "04_elbow.pdf"),             width = 6,  height = 4)
save_plot(p_pca_patient,     file.path(dimred_dir, "04_pca_pre_harmony.pdf"),   width = 8,  height = 6)
save_plot(p_harmony_patient, file.path(dimred_dir, "04_harmony_embedding.pdf"), width = 8,  height = 6)
save_plot(p_umap_patient + p_umap_pfs,
          file.path(dimred_dir, "04_umap_patient_pfs.pdf"), width = 16, height = 7)
save_plot(p_umap_sample,     file.path(dimred_dir, "04_umap_roi.pdf"),          width = 8,  height = 7)
save_plot(feat_plots,        file.path(dimred_dir, "04_umap_features.pdf"),     width = 16, height = 4 * ceiling(length(key_markers) / 4))

# ── 7. Save ──────────────────────────────────────────────────────────────────
saveRDS(seurat_obj, file.path(RDS_DIR, "04_seurat_dimred.rds"))
message("\nStep 4 complete. Saved: output/data/04_seurat_dimred.rds")
