# =============================================================================
# 05b_annotation.R
# Manual cell type annotation from heatmap + corrected feature plots.
# Based on 05_cluster_heatmap.pdf (avg normalized MFI, z-scored across 20 clusters)
# NOTE: Step 06 (gating) is independent — reads from 05, not 05b.
#       This step is for visualization and biological interpretation of clusters.
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))
source(file.path(getwd(), "code/utils/plot_utils.R"))

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(readr)

message("=== Step 5b: Manual annotation (20 clusters) ===")

seurat_obj <- readRDS(file.path(RDS_DIR, "05_seurat_clustered.rds"))
markers    <- rownames(seurat_obj)

# ── 1. Manual annotation from heatmap ────────────────────────────────────────
# Based on 05_cluster_heatmap.pdf (res=0.15, 20 clusters, 262k cells):
#
#   g0  → GFAP+, Vimentin+, CD44+ (largest cluster)                → Tumor
#   g1  → moderate mixed expression                                 → Tumor
#   g2  → CD163++, CD68++, TMEM119+, CD73+                         → TAM_M2
#   g3  → moderate expression, no strong immune signature           → Tumor
#   g4  → moderate, CD44+                                           → Tumor
#   g5  → GD2++++, p53+++, H3-3G34V+++, EZHIP+++, Vimentin++      → Tumor_G34
#   g6  → CD8a++++, CD57++, CD45RO++, Ki-67++, CCR7+               → T_CD8
#   g7  → moderate tumor markers                                    → Tumor
#   g8  → CD3+++, CD45+++, PD-1+++, TIM3+++, FoxP3+, BATF+, CD4+  → T_exhausted
#   g9  → NKG2D+, B7H3+, low immune                                → Tumor
#   g10 → CD105++ (Endoglin)                                        → Endothelial
#   g11 → very low expression                                       → Unclassified
#   g12 → moderate, near myeloid branch                             → Tumor
#   g13 → low expression                                            → Tumor
#   g14 → BCAN+, CD80+                                              → Tumor_inflamed
#   g15 → CD80++, HLA-DR++, CD45RA+, GATA3+                        → DC_activated
#   g16 → low expression                                            → Tumor
#   g17 → CD206+++, CD163++, CD68++, TMEM119++, CD73+              → Microglia
#   g18 → Podoplanin+, low immune                                   → Stromal
#   g19 → very low expression                                       → Unclassified

MANUAL_ANNOT <- c(
  "0"  = "Tumor",
  "1"  = "Tumor",
  "2"  = "TAM_M2",
  "3"  = "Tumor",
  "4"  = "Tumor",
  "5"  = "Tumor_G34",
  "6"  = "T_CD8",
  "7"  = "Tumor",
  "8"  = "T_exhausted",
  "9"  = "Tumor",
  "10" = "Endothelial",
  "11" = "Unclassified",
  "12" = "Tumor",
  "13" = "Tumor",
  "14" = "Tumor_inflamed",
  "15" = "DC_activated",
  "16" = "Tumor",
  "17" = "Microglia",
  "18" = "Stromal",
  "19" = "Unclassified"
)

# Assign to Seurat object
celltype_vec <- MANUAL_ANNOT[as.character(seurat_obj$cluster)]
names(celltype_vec) <- colnames(seurat_obj)
seurat_obj <- AddMetaData(seurat_obj, metadata = celltype_vec, col.name = "celltype")

# Ordered levels for consistent plotting
CT_LEVELS <- c(
  "Tumor", "Tumor_G34", "Tumor_inflamed",
  "T_CD8", "T_exhausted",
  "TAM_M2", "Microglia", "DC_activated",
  "Endothelial", "Stromal", "Unclassified"
)
seurat_obj$celltype <- factor(seurat_obj$celltype, levels = CT_LEVELS)

CELLTYPE_COLORS_MANUAL <- c(
  Tumor            = "#E41A1C",
  Tumor_G34        = "#FF6666",
  Tumor_inflamed   = "#FDAE61",
  T_CD8            = "#377EB8",
  T_exhausted      = "#1B3A6B",
  TAM_M2           = "#4DAF4A",
  Microglia        = "#984EA3",
  DC_activated     = "#FF7F00",
  Endothelial      = "#A65628",
  Stromal          = "#B3B3B3",
  Unclassified     = "#999999"
)

message("\n--- Manual annotation ---")
ct_counts <- table(seurat_obj$celltype)
for (ct in names(ct_counts)) {
  pct <- round(100 * ct_counts[[ct]] / sum(ct_counts), 1)
  message(sprintf("  %-20s : %6d cells (%4.1f%%)", ct, ct_counts[[ct]], pct))
}

# ── 2. UMAP plots with manual annotation ─────────────────────────────────────
p_umap_ct <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype",
                     label = TRUE, label.size = 3.5, pt.size = PARAMS$point_size,
                     repel = TRUE) +
  scale_color_manual(values = CELLTYPE_COLORS_MANUAL) +
  labs(title = "UMAP — Cell type (manual annotation)", color = "Cell type") +
  theme_hgg()

p_umap_patient <- DimPlot(seurat_obj, reduction = "umap", group.by = "patient_id",
                           pt.size = PARAMS$point_size) +
  scale_color_manual(values = PATIENT_COLORS) +
  labs(title = "UMAP — Patient", color = "Patient") +
  theme_hgg()

save_plot(p_umap_ct + p_umap_patient,
          file.path(FIG_DIR, "clustering", "05b_umap_annotated.pdf"),
          width = 16, height = 7)

# ── 3. Patient composition per cluster (key diagnostic plot) ─────────────────
# Who contributes to each cluster? Reveals patient-specific vs shared clusters.
comp_patient_cluster <- seurat_obj@meta.data %>%
  group_by(cluster, patient_id) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(frac = n / sum(n),
         celltype = MANUAL_ANNOT[as.character(cluster)]) %>%
  ungroup() %>%
  mutate(celltype = factor(celltype, levels = CT_LEVELS))

p_patient_per_cluster <- ggplot(comp_patient_cluster,
                                aes(x = reorder(paste0(cluster, "\n", celltype), as.numeric(cluster)),
                                    y = frac, fill = patient_id)) +
  geom_col(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = PATIENT_COLORS) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Cluster (cell type)", y = "Fraction", fill = "Patient",
       title = "Patient composition per cluster — which patient drives each cluster?") +
  theme_hgg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

save_plot(p_patient_per_cluster,
          file.path(FIG_DIR, "clustering", "05b_patient_per_cluster.pdf"),
          width = 12, height = 6)

# ── 4. Compositional analysis with correct annotation ─────────────────────────
comp_df <- seurat_obj@meta.data %>%
  group_by(patient_id, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

p_comp <- ggplot(comp_df, aes(x = patient_id, y = fraction, fill = celltype)) +
  geom_col(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = CELLTYPE_COLORS_MANUAL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Patient", y = "Fraction of cells", fill = "Cell type",
       title = "TME composition per patient (blind to outcome)") +
  theme_hgg()

save_plot(p_comp,
          file.path(FIG_DIR, "clustering", "05b_composition_annotated.pdf"),
          width = 10, height = 6)

# Heatmap: patient × cell type fraction
comp_wide <- comp_df %>%
  select(patient_id, celltype, fraction) %>%
  mutate(celltype = as.character(celltype)) %>%
  pivot_wider(names_from = celltype, values_from = fraction, values_fill = 0) %>%
  tibble::column_to_rownames("patient_id")

# Reorder columns
comp_wide <- comp_wide[, intersect(CT_LEVELS, colnames(comp_wide)), drop = FALSE]

pheatmap(as.matrix(comp_wide),
         color            = colorRampPalette(c("white", "#2166AC"))(100),
         display_numbers  = TRUE,
         number_format    = "%.2f",
         fontsize_row     = 11, fontsize_col = 9,
         cluster_rows     = TRUE,
         cluster_cols     = FALSE,
         clustering_method = "ward.D2",
         main   = "TME composition — patient x cell type fractions (blind to outcome)",
         filename = file.path(FIG_DIR, "clustering", "05b_composition_heatmap.pdf"),
         width = 12, height = 5)

write_tsv(comp_df, file.path(TAB_DIR, "05b_composition_annotated.tsv"))

# ── 5. Feature plots — corrected with proper scaling ─────────────────────────
# Fix: use the arcsinh-normalized values from the "data" layer explicitly
# Use min.cutoff/max.cutoff to clip extremes and improve color dynamic range

key_markers_canonical <- c(
  "CD3", "CD8a", "CD4", "CD45", "CD279",
  "TIM3", "FoxP3", "GFAP", "GD2", "H3-3G34V",
  "CD68", "CD163", "CD206", "TMEM119", "CD80",
  "HLA-DR", "CD105", "Podoplanin", "CD274", "B7H3"
)
key_markers_seurat <- intersect(key_markers_canonical, markers)

feat_plots <- FeaturePlot(
  seurat_obj,
  features   = key_markers_seurat,
  reduction  = "umap",
  pt.size    = PARAMS$point_size,
  ncol       = 5,
  order      = TRUE,
  min.cutoff = "q5",    # clip bottom 5% → improves color dynamic range
  max.cutoff = "q95"    # clip top 5%
) &
  scale_color_gradientn(colours = c("grey90", "#4292C6", "#08306B")) &
  theme_hgg() &
  theme(legend.position = "right", plot.title = element_text(size = 9))

save_plot(feat_plots,
          file.path(FIG_DIR, "dimred", "05b_umap_features_corrected.pdf"),
          width = 18, height = 3.5 * ceiling(length(key_markers_seurat) / 5))

# ── 6. Heatmap with manual annotation labels ──────────────────────────────────
avg_expr <- AverageExpression(seurat_obj, assay = "MFI", layer = "data",
                               group.by = "celltype", verbose = FALSE)$MFI

marker_var  <- apply(avg_expr, 1, var)
top_markers <- names(sort(marker_var, decreasing = TRUE))[1:min(45, nrow(avg_expr))]
top_markers <- top_markers[!grepl("DAPI", top_markers)]

pheatmap(avg_expr[top_markers, ],
         scale             = "row",
         clustering_method = "ward.D2",
         color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         fontsize_row      = 8,
         fontsize_col      = 9,
         main              = "Avg normalized MFI per cell type (z-scored)",
         filename          = file.path(FIG_DIR, "clustering", "05b_celltype_heatmap.pdf"),
         width = 10, height = 11)

# ── 7. Save annotated Seurat object ──────────────────────────────────────────
saveRDS(seurat_obj, file.path(RDS_DIR, "05b_seurat_annotated.rds"))
message("\nStep 5b complete. Saved: output/data/05b_seurat_annotated.rds")
