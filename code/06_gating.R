# =============================================================================
# 06_gating.R — Rule-based cell type gating for MACSima imaging data
#
# Gating tree: config/gating_tree.yaml (single source; reference diagram
# info/gating strategy for the analysis.png; adaptations commented there;
# diagram of the tree: docs/gating_tree.md). Per-patient GMM thresholds for
# each gating marker; the tree assigns one cell type per cell.
#
# Output: cell type per cell + fractions per patient
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))
source(file.path(getwd(), "code/utils/plot_utils.R"))
source(file.path(getwd(), "code/utils/gating_tree.R"))

library(Seurat)
library(mclust)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)

message("=== Step 6: Rule-based gating ===")

GATING <- read_gating_tree(file.path(getwd(), "config/gating_tree.yaml"))
THR    <- GATING$thresholds

seurat_obj <- readRDS(file.path(RDS_DIR, "05_seurat_clustered.rds"))

# Extract normalized data as data.frame (cells x markers)
mfi_norm <- t(as.matrix(GetAssayData(seurat_obj, assay = "MFI", layer = "data")))
mfi_df   <- as.data.frame(mfi_norm)
mfi_df$cell_id    <- rownames(mfi_df)
mfi_df$patient_id <- seurat_obj$patient_id
mfi_df$sample_id  <- seurat_obj$sample_id
mfi_df$roi_id     <- seurat_obj$roi_id

markers_seurat <- colnames(mfi_norm)

# ── 1. GMM THRESHOLDS PER PATIENT PER MARKER ─────────────────────────────────

find_gmm_threshold <- function(values, marker_name) {
  values <- values[is.finite(values)]
  # Mclust starts from a random subset of cells: a fixed seed per fit makes every threshold
  # reproducible and independent of the order of the fits (the global random state is restored).
  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  }
  set.seed(THR$seed)
  tryCatch({
    fit <- mclust::Mclust(values, G = 2, verbose = FALSE)
    if (is.null(fit)) stop("GMM failed")

    means  <- fit$parameters$mean
    sds    <- sqrt(fit$parameters$variance$sigmasq)
    props  <- fit$parameters$pro

    ord    <- order(means)
    mu_neg <- means[ord[1]];  sd_neg <- sds[min(ord[1], length(sds))]
    mu_pos <- means[ord[2]];  sd_pos <- sds[min(ord[2], length(sds))]
    p_neg  <- props[ord[1]];  p_pos  <- props[ord[2]]

    x_seq  <- seq(mu_neg, mu_pos, length.out = 1000)
    d_neg  <- p_neg * dnorm(x_seq, mu_neg, sd_neg)
    d_pos  <- p_pos * dnorm(x_seq, mu_pos, sd_pos)
    thresh <- x_seq[which.min(abs(d_neg - d_pos))]

    if (thresh <= mu_neg || thresh >= mu_pos) {
      thresh <- (mu_neg + mu_pos) / 2
    }
    list(threshold = thresh, method = "GMM")
  }, error = function(e) {
    list(threshold = quantile(values, THR$on_fit_error), method = "p90")
  })
}

# Gating markers: logical name -> Seurat feature name (config/gating_tree.yaml)
gating_name_map <- unlist(GATING$markers)

# Every gating marker must be in the Seurat object: a missing one would silently
# empty every population gated on it.
missing_markers <- setdiff(gating_name_map, markers_seurat)
if (length(missing_markers) > 0) {
  stop("Gating markers not in the Seurat object: ", paste(missing_markers, collapse = ", "))
}
message(sprintf("Gating markers available: %d / %d", length(gating_name_map), length(gating_name_map)))

patients <- unique(mfi_df$patient_id)

message("Computing per-patient GMM thresholds...")
thresholds_per_patient <- list()
for (pat in patients) {
  thresholds_per_patient[[pat]] <- list()
  cells_pat <- mfi_df$patient_id == pat
  for (logic_name in names(gating_name_map)) {
    seurat_name <- gating_name_map[[logic_name]]
    vals        <- mfi_df[[seurat_name]][cells_pat]
    result      <- find_gmm_threshold(vals, seurat_name)

    pct_pos <- mean(vals > result$threshold, na.rm = TRUE)
    if (pct_pos > THR$max_positive) {
      result$threshold <- quantile(vals, THR$fallback_quantile, na.rm = TRUE)
      result$method    <- "p85_fallback"
    }
    thresholds_per_patient[[pat]][[logic_name]] <- result
  }
}

# Report thresholds
message("\n  Marker       | ", paste(sprintf("%-10s", patients), collapse=" | "), " | method")
for (logic_name in names(gating_name_map)) {
  thresh_vals <- sapply(patients, function(p)
    sprintf("%.2f", thresholds_per_patient[[p]][[logic_name]]$threshold))
  methods <- sapply(patients, function(p)
    thresholds_per_patient[[p]][[logic_name]]$method)
  method_str <- ifelse(all(methods == "GMM"), "GMM", paste(unique(methods), collapse="/"))
  message(sprintf("  %-12s | %s | %s",
                  logic_name,
                  paste(sprintf("%-10s", thresh_vals), collapse=" | "),
                  method_str))
}

# ── 2. THRESHOLD VALIDATION PLOTS ────────────────────────────────────────────
plot_threshold <- function(logic_name) {
  seurat_name <- gating_name_map[[logic_name]]
  df_plot <- data.frame(value = mfi_df[[seurat_name]], patient_id = mfi_df$patient_id)

  thresh_df <- data.frame(
    patient_id = patients,
    threshold  = sapply(patients, function(p)
      thresholds_per_patient[[p]][[logic_name]]$threshold)
  )

  ggplot(df_plot, aes(x = value, color = patient_id)) +
    geom_density(linewidth = 0.7, alpha = 0.8) +
    geom_vline(data = thresh_df,
               aes(xintercept = threshold, color = patient_id),
               linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = PATIENT_COLORS) +
    labs(title = sprintf("%s (%s)", logic_name, seurat_name),
         x = "arcsinh(MFI/cofactor)", y = "Density", color = "Patient") +
    theme_hgg() +
    theme(legend.position = "none", plot.title = element_text(size = 9))
}

thresh_plots <- lapply(names(gating_name_map), plot_threshold)
thresh_panel <- wrap_plots(thresh_plots, ncol = 4)
save_plot(thresh_panel,
          file.path(FIG_DIR, "qc", "06_gating_thresholds.pdf"),
          width = 16, height = 4 * ceiling(length(gating_name_map) / 4))

# ── 3. APPLY QC FLAGS ────────────────────────────────────────────────────────
flag_cols <- grep("^flag_", colnames(seurat_obj@meta.data), value = TRUE)
flag_mat  <- as.data.frame(seurat_obj@meta.data[, flag_cols, drop = FALSE])

# ── 4. GATING TREE ───────────────────────────────────────────────────────────
# Positivity per gating marker (per-patient thresholds); QC-flagged = NA.
compute_positivity <- function() {
  n   <- nrow(mfi_df)
  pos <- list()
  for (lname in names(gating_name_map)) {
    sname  <- gating_name_map[[lname]]
    vals   <- mfi_df[[sname]]
    result <- rep(NA, n)

    for (pat in patients) {
      idx    <- mfi_df$patient_id == pat
      thresh <- thresholds_per_patient[[pat]][[lname]]$threshold
      result[idx] <- vals[idx] > thresh
    }

    # Mask QC-flagged marker/sample combinations
    flag_col <- sub("{canonical}", gsub("-", "_", sname), GATING$flags$column, fixed = TRUE)
    if (flag_col %in% colnames(flag_mat)) {
      is_flagged         <- !is.na(flag_mat[[flag_col]]) & flag_mat[[flag_col]]
      result[is_flagged] <- NA
    }
    pos[[lname]] <- result
  }
  pos
}

message("Applying gating tree...")
mfi_df$celltype_gating <- apply_gating_tree(GATING, compute_positivity(), nrow(mfi_df))

# ── 4b. T CELL FUNCTIONAL ANNOTATIONS (metadata, not exclusive types) ────────
# config/gating_tree.yaml `annotations`: HLA-DR+ = activated, Ki-67+ = proliferating.
# These overlay on differentiation states and are stored as separate columns.
message("Adding T cell functional annotations...")

for (ann in GATING$annotations) {
  is_target <- grepl(ann$labels_regex, mfi_df$celltype_gating)
  col       <- rep(NA, nrow(mfi_df))
  for (pat in patients) {
    idx_pat  <- mfi_df$patient_id == pat & is_target
    thresh   <- thresholds_per_patient[[pat]][[ann$marker]]$threshold
    col[idx_pat] <- mfi_df[[gating_name_map[[ann$marker]]]][idx_pat] > thresh
  }
  mfi_df[[ann$column]] <- col
}

is_any_T <- grepl(GATING$annotations[[1]]$labels_regex, mfi_df$celltype_gating)
n_T <- sum(is_any_T)
n_act <- sum(mfi_df$T_activated, na.rm = TRUE)
n_prol <- sum(mfi_df$T_proliferating, na.rm = TRUE)
message(sprintf("  T cells: %d total, %d activated (HLA-DR+), %d proliferating (Ki-67+)",
                n_T, n_act, n_prol))

# ── 5. SUMMARY ────────────────────────────────────────────────────────────────
ct_table <- table(mfi_df$celltype_gating)
ct_df    <- as.data.frame(ct_table) %>%
  rename(celltype = Var1, n_cells = Freq) %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 1)) %>%
  arrange(desc(n_cells))

message("\n--- Gating results (all patients combined) ---")
print(ct_df)
write_tsv(ct_df, file.path(TAB_DIR, "06_gating_summary.tsv"))

# ── 6. COMPOSITION PER PATIENT ────────────────────────────────────────────────
CT_ORDER_GATING  <- names(GATING$labels)
CT_COLORS_GATING <- unlist(GATING$labels)

comp_gating <- mfi_df %>%
  mutate(celltype_gating = factor(celltype_gating, levels = CT_ORDER_GATING)) %>%
  group_by(patient_id, celltype_gating) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

write_tsv(comp_gating, file.path(TAB_DIR, "06_composition_gating.tsv"))

# Add PFS group for plotting
comp_gating <- comp_gating %>%
  left_join(CLINICAL[, c("patient_id", "pfs_group", "pfs_months")], by = "patient_id")

# Stacked bar per patient, ordered by PFS
p_comp <- ggplot(comp_gating,
                 aes(x = reorder(patient_id, -pfs_months),
                     y = fraction, fill = celltype_gating)) +
  geom_col(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = CT_COLORS_GATING, breaks = CT_ORDER_GATING, drop = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Patient (ordered by PFS)", y = "Fraction of cells", fill = "Cell type",
       title = "TME composition per patient - rule-based gating") +
  theme_hgg() +
  guides(fill = guide_legend(ncol = 2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_comp,
          file.path(FIG_DIR, "clustering", "06_composition_gating.pdf"),
          width = 14, height = 8)

# Heatmap patient x cell type
comp_wide <- comp_gating %>%
  select(patient_id, celltype_gating, fraction) %>%
  mutate(celltype_gating = as.character(celltype_gating)) %>%
  pivot_wider(names_from = celltype_gating, values_from = fraction, values_fill = 0) %>%
  tibble::column_to_rownames("patient_id")
comp_wide <- comp_wide[, intersect(CT_ORDER_GATING, colnames(comp_wide)), drop = FALSE]

# Annotation: PFS group
annot_row <- CLINICAL[, c("patient_id", "pfs_group", "pfs_months")] %>%
  tibble::column_to_rownames("patient_id")
annot_row <- annot_row[rownames(comp_wide), , drop = FALSE]
annot_colors <- list(pfs_group = PFS_GROUP_COLORS)

pheatmap(as.matrix(comp_wide),
         color = colorRampPalette(c("white", "#2166AC"))(100),
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_row = 11, fontsize_col = 7,
         cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_method = "ward.D2",
         annotation_row = annot_row,
         annotation_colors = annot_colors,
         main = "Patient x cell type fractions - gating",
         filename = file.path(FIG_DIR, "clustering", "06_composition_heatmap_gating.pdf"),
         width = 18, height = 7)

# ── 7. ADD GATING LABELS TO SEURAT ───────────────────────────────────────────
seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = setNames(mfi_df$celltype_gating, mfi_df$cell_id),
  col.name = "celltype_gating"
)

# Add T cell functional annotations
seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = setNames(mfi_df$T_activated, mfi_df$cell_id),
  col.name = "T_activated"
)
seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = setNames(mfi_df$T_proliferating, mfi_df$cell_id),
  col.name = "T_proliferating"
)

# UMAP gating vs clusters
p_umap_gating <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "celltype_gating",
  pt.size = PARAMS$point_size, label = FALSE
) +
  scale_color_manual(values = CT_COLORS_GATING, na.value = "grey80") +
  labs(title = "UMAP - cell type (rule-based gating)", color = "Cell type") +
  theme_hgg()

p_umap_cluster <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "cluster",
  pt.size = PARAMS$point_size, label = TRUE, label.size = 3
) +
  labs(title = "UMAP - Louvain clusters", color = "Cluster") +
  theme_hgg()

save_plot(p_umap_gating + p_umap_cluster,
          file.path(FIG_DIR, "clustering", "06_umap_gating_vs_clustering.pdf"),
          width = 20, height = 7)

# ── 8. SAVE ──────────────────────────────────────────────────────────────────
saveRDS(seurat_obj, file.path(RDS_DIR, "06_seurat_gated.rds"))
saveRDS(mfi_df,     file.path(RDS_DIR, "06_mfi_gated.rds"))

# ── 9. GATING DEVIATIONS TABLE ──────────────────────────────────────────────
# Where the implementation differs from the reference diagram (config/gating_tree.yaml)
deviations <- do.call(rbind, lapply(GATING$deviations, as.data.frame, stringsAsFactors = FALSE))
write_tsv(deviations, file.path(TAB_DIR, "06_gating_deviations.tsv"))

message("\nStep 6 complete.")
