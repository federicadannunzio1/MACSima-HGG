# =============================================================================
# 05_clustering.R
# k-NN graph → Louvain clustering at multiple resolutions → clustree validation
# → automatic cell type annotation → compositional analysis per patient.
# Output: Seurat with cluster labels, cell type annotations, composition plots.
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))
source(file.path(getwd(), "code/utils/plot_utils.R"))

library(Seurat)
library(clustree)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)

message("=== Step 5: Clustering + annotation ===")

seurat_obj <- readRDS(file.path(RDS_DIR, "04_seurat_dimred.rds"))
markers    <- rownames(seurat_obj)

# Update CELLTYPE_RULES to use Seurat's feature names (_ → -)
normalize_marker <- function(x) gsub("_", "-", x)
CELLTYPE_RULES_SEURAT <- lapply(CELLTYPE_RULES, function(rule) {
  list(
    primary   = normalize_marker(rule$primary),
    secondary = normalize_marker(rule$secondary)
  )
})
# Patch specific name differences
CELLTYPE_RULES_SEURAT$B_cell$primary   <- c("CD20-Cytoplasmic")

# ── 1. k-NN graph ─────────────────────────────────────────────────────────────
message("Building k-NN graph on Harmony embedding...")
seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims      = 1:PARAMS$n_pcs,
  verbose   = FALSE
)

# ── 2. Louvain at multiple resolutions ───────────────────────────────────────
message("Running Louvain clustering...")
for (res in PARAMS$louvain_resolutions) {
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution   = res,
    algorithm    = 1,
    random.seed  = 42,
    verbose      = FALSE
  )
  col_to <- sprintf("louvain_res%s", gsub("\\.", "", sprintf("%.2f", res)))
  seurat_obj@meta.data[[col_to]] <- seurat_obj@meta.data[["seurat_clusters"]]
  n_clust <- length(unique(seurat_obj@meta.data[[col_to]]))
  message(sprintf("  res=%.2f -> %d clusters  [col: %s]", res, n_clust, col_to))
}

# ── 3. Clustree ──────────────────────────────────────────────────────────────
res_cols <- sprintf("louvain_res%s", gsub("\\.", "", sprintf("%.2f", PARAMS$louvain_resolutions)))
for (rc in res_cols) {
  seurat_obj@meta.data[[rc]] <- factor(seurat_obj@meta.data[[rc]])
}

p_clustree <- clustree(seurat_obj@meta.data, prefix = "louvain_res") +
  labs(title = "Clustree - stability across Louvain resolutions") +
  theme(plot.title = element_text(face = "bold"))

save_plot(p_clustree,
          file.path(FIG_DIR, "clustering", "05_clustree.pdf"),
          width = 10, height = 8)

# ── 4. Choose working resolution ─────────────────────────────────────────────
# Default: res=0.15 (good balance for 8 patients, ~260k cells)
# Inspect clustree to decide; can be changed in 05b
WORKING_RES    <- "louvain_res015"
seurat_obj$cluster <- seurat_obj@meta.data[[WORKING_RES]]
Idents(seurat_obj) <- "cluster"

n_clusters <- nlevels(seurat_obj$cluster)
message(sprintf("Working resolution: 0.15 -> %d clusters", n_clusters))

# UMAP colored by cluster
p_umap_clust <- DimPlot(seurat_obj, reduction = "umap", group.by = "cluster",
                         label = TRUE, label.size = 4, pt.size = PARAMS$point_size) +
  labs(title = sprintf("UMAP - Louvain clusters (res=0.15, n=%d)", n_clusters), color = "Cluster") +
  theme_hgg()

# All resolutions side by side
umap_res_plots <- lapply(res_cols, function(rc) {
  DimPlot(seurat_obj, reduction = "umap", group.by = rc,
          pt.size = PARAMS$point_size * 0.5, label = TRUE, label.size = 3) +
    labs(title = paste0("res=", sub("louvain_res0*", "0.", rc))) +
    theme_hgg() + NoLegend()
})
save_plot(wrap_plots(umap_res_plots, ncol = 3),
          file.path(FIG_DIR, "clustering", "05_umap_all_resolutions.pdf"),
          width = 18, height = 12)

# ── 5. Marker expression heatmap per cluster ──────────────────────────────────
avg_expr <- AverageExpression(seurat_obj, assay = "MFI", layer = "data",
                               group.by = "cluster", verbose = FALSE)$MFI

cluster_levels <- levels(seurat_obj$cluster)
avg_colnames   <- colnames(avg_expr)
cluster_to_col <- setNames(avg_colnames, cluster_levels)

marker_var  <- apply(avg_expr, 1, var)
top_markers <- names(sort(marker_var, decreasing = TRUE))[1:min(50, nrow(avg_expr))]
top_markers <- top_markers[!grepl("DAPI", top_markers)]

pheatmap(avg_expr[top_markers, ],
         scale        = "row",
         clustering_method = "ward.D2",
         color        = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         fontsize_row = 7,
         fontsize_col = 9,
         main         = "Avg normalized MFI per cluster (z-scored across clusters)",
         filename     = file.path(FIG_DIR, "clustering", "05_cluster_heatmap.pdf"),
         width = 12, height = 14)

# ── 6. Automatic cell type annotation ─────────────────────────────────────────
EXPR_THRESHOLD <- 0.5

annotate_cluster <- function(cluster_id) {
  col_id <- cluster_to_col[[as.character(cluster_id)]]
  scores <- sapply(names(CELLTYPE_RULES_SEURAT), function(ct) {
    rule       <- CELLTYPE_RULES_SEURAT[[ct]]
    prim_avail <- intersect(rule$primary, markers)
    if (length(prim_avail) == 0) return(0)
    prim_vals  <- avg_expr[prim_avail, col_id, drop = TRUE]
    mean(prim_vals > EXPR_THRESHOLD)
  })
  best_ct <- names(which.max(scores))
  if (max(scores) >= 0.5) best_ct else "Unknown"
}

cluster_ids   <- cluster_levels
cluster_annot <- setNames(sapply(cluster_ids, annotate_cluster), cluster_ids)

message("\n--- Automatic cluster annotations ---")
for (cl in cluster_ids) {
  message(sprintf("  Cluster %s -> %s", cl, cluster_annot[cl]))
}

celltype_vec <- cluster_annot[as.character(seurat_obj$cluster)]
names(celltype_vec) <- colnames(seurat_obj)
seurat_obj <- AddMetaData(seurat_obj, metadata = celltype_vec, col.name = "celltype_auto")

p_umap_ct <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype_auto",
                      label = TRUE, label.size = 3.5, pt.size = PARAMS$point_size,
                      repel = TRUE) +
  scale_color_manual(values = CELLTYPE_COLORS, na.value = "grey70") +
  labs(title = "UMAP - Cell type (automatic)", color = "Cell type") +
  theme_hgg()

save_plot(p_umap_clust + p_umap_ct,
          file.path(FIG_DIR, "clustering", "05_umap_clusters_celltypes.pdf"),
          width = 16, height = 7)

# ── 7. Compositional analysis ─────────────────────────────────────────────────
comp_df <- seurat_obj@meta.data %>%
  group_by(patient_id, celltype_auto) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

p_comp <- ggplot(comp_df, aes(x = patient_id, y = fraction, fill = celltype_auto)) +
  geom_col(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = CELLTYPE_COLORS, na.value = "grey70") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Patient", y = "Fraction of cells", fill = "Cell type",
       title = "TME composition per patient (automatic annotation)") +
  theme_hgg() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_comp,
          file.path(FIG_DIR, "clustering", "05_composition_per_patient.pdf"),
          width = 10, height = 6)

write_tsv(comp_df, file.path(TAB_DIR, "05_composition_per_patient.tsv"))

# ── 8. Save ──────────────────────────────────────────────────────────────────
saveRDS(seurat_obj, file.path(RDS_DIR, "05_seurat_clustered.rds"))
message("\nStep 5 complete. Saved: output/data/05_seurat_clustered.rds")
message("NOTE: Inspect 05_cluster_heatmap.pdf and update 05b_annotation.R with manual labels.")
