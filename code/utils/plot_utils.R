# =============================================================================
# plot_utils.R — Shared plotting helpers
# =============================================================================

library(ggplot2)
library(patchwork)

# Clean ggplot2 theme for all plots
theme_hgg <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey92", color = NA),
      plot.title        = element_text(size = base_size + 1, face = "bold"),
      legend.key.size   = unit(0.4, "cm"),
      legend.text       = element_text(size = base_size - 2)
    )
}

# Save a ggplot or patchwork object
save_plot <- function(plot, filename, width = 8, height = 6, dpi = PARAMS$dpi) {
  ggsave(filename, plot = plot, width = width, height = height,
         dpi = dpi, bg = "white")
  message("Saved: ", filename)
}

# DimPlot wrapper with consistent style
dim_plot <- function(seurat_obj, group_by, colors = NULL, title = NULL,
                     pt_size = PARAMS$point_size, reduction = "umap") {
  p <- DimPlot(seurat_obj, reduction = reduction, group.by = group_by,
               pt.size = pt_size, label = FALSE) +
    theme_hgg() +
    labs(title = title %||% group_by, color = group_by)
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  p
}

# Violin plot helper for QC metrics
qc_violin <- function(seurat_obj, feature, group_by = "sample_id",
                      y_lab = feature, fill_col = "steelblue") {
  df <- seurat_obj@meta.data
  df$feature_val <- df[[feature]]
  ggplot(df, aes(x = .data[[group_by]], y = feature_val)) +
    geom_violin(fill = fill_col, alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
    labs(x = NULL, y = y_lab) +
    theme_hgg() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
