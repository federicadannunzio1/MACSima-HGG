# =============================================================================
# 06_gating.R — Rule-based cell type gating for MACSima imaging data
#
# Per ogni marker chiave, soglia di positivita' via GMM a 2 componenti
# (per-paziente). Poi cascata di decisioni gerarchiche (gating tree).
#
# Output: cell type per cellula + frazioni per paziente
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))
source(file.path(getwd(), "code/utils/plot_utils.R"))

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
    list(threshold = quantile(values, 0.90), method = "p90")
  })
}

# Gating markers
gating_name_map <- c(
  CD45 = "CD45", CD3 = "CD3", CD8a = "CD8a", CD4 = "CD4",
  FoxP3 = "FoxP3", TIGIT = "TIGIT", PD1 = "CD279",
  CD56 = "CD56", CD68 = "CD68", P2Y12 = "P2Y12",
  CD163 = "CD163", CD80 = "CD80",
  GFAP = "GFAP", BCAN = "BCAN", CD31 = "CD31", Vimentin = "Vimentin"
)
gating_name_map <- gating_name_map[gating_name_map %in% markers_seurat]

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
    if (pct_pos > 0.60) {
      result$threshold <- quantile(vals, 0.85, na.rm = TRUE)
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
gate_cell_vectorized <- function() {
  n <- nrow(mfi_df)

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

    flag_col <- paste0("flag_", lname)
    if (flag_col %in% colnames(flag_mat)) {
      is_flagged         <- !is.na(flag_mat[[flag_col]]) & flag_mat[[flag_col]]
      result[is_flagged] <- NA
    }
    pos[[lname]] <- result
  }

  p <- function(marker) { ifelse(is.na(pos[[marker]]), FALSE, pos[[marker]]) }

  celltype <- rep("Unclassified", n)

  # IMMUNE branch (CD45+)
  is_immune <- if ("CD45" %in% names(pos)) p("CD45") else rep(FALSE, n)

  n_cd45_na <- sum(is.na(pos[["CD45"]]))
  if (n_cd45_na > 0.3 * n) {
    message(sprintf("  WARNING: CD45 NA in %.0f%% of cells - using CD3/CD68/CD56 as immune proxy",
                    100 * n_cd45_na / n))
    is_immune <- p("CD3") | p("CD68") | p("CD56")
  }

  # T cells
  is_T    <- is_immune & p("CD3")
  is_CD8  <- is_T & p("CD8a")
  is_CD4  <- is_T & p("CD4") & !is_CD8
  is_Toth <- is_T & !is_CD8 & !is_CD4

  celltype[is_CD8 & (p("TIGIT") | p("PD1"))] <- "T_CD8_exhausted"
  celltype[is_CD8 & !(p("TIGIT") | p("PD1"))] <- "T_CD8_effector"
  celltype[is_CD4 & p("FoxP3")] <- "T_reg"
  celltype[is_CD4 & !p("FoxP3")] <- "T_CD4_helper"
  celltype[is_Toth] <- "T_other"

  # Myeloid / NK
  is_nonT_immune <- is_immune & !p("CD3")
  is_NK  <- is_nonT_immune & p("CD56")
  is_My  <- is_nonT_immune & !p("CD56") & p("CD68")

  celltype[is_NK] <- "NK"
  celltype[is_My & p("P2Y12")]                                   <- "Microglia"
  celltype[is_My & !p("P2Y12") & p("CD163")]                     <- "TAM_M2"
  celltype[is_My & !p("P2Y12") & !p("CD163") & p("CD80")]        <- "DC_M1"
  celltype[is_My & !p("P2Y12") & !p("CD163") & !p("CD80")]       <- "TAM_other"
  celltype[is_nonT_immune & !p("CD56") & !p("CD68")]              <- "Other_immune"

  # NON-IMMUNE branch (CD45-)
  is_nonimmune <- !is_immune

  is_tumor <- is_nonimmune & (p("GFAP") | p("BCAN"))
  celltype[is_tumor] <- "Tumor"

  is_endo  <- is_nonimmune & !is_tumor & p("CD31")
  celltype[is_endo] <- "Endothelial"

  is_stro  <- is_nonimmune & !is_tumor & !is_endo & p("Vimentin")
  celltype[is_stro] <- "Stromal"

  celltype[is_nonimmune & !is_tumor & !is_endo & !is_stro] <- "Unclassified_nonimmune"

  celltype
}

message("Applying gating tree...")
mfi_df$celltype_gating <- gate_cell_vectorized()

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
CT_ORDER_GATING <- c(
  "Tumor", "Endothelial", "Stromal", "Unclassified_nonimmune",
  "T_CD8_exhausted", "T_CD8_effector", "T_CD4_helper", "T_reg", "T_other",
  "NK", "Microglia", "TAM_M2", "TAM_other", "DC_M1",
  "Other_immune", "Unclassified"
)
CT_COLORS_GATING <- c(
  Tumor = "#E41A1C", Endothelial = "#A65628", Stromal = "#FDBF6F",
  Unclassified_nonimmune = "#CCCCCC",
  T_CD8_exhausted = "#08519C", T_CD8_effector = "#4292C6",
  T_CD4_helper = "#74ADD1", T_reg = "#ABD9E9", T_other = "#C6DBEF",
  NK = "#FF7F00", Microglia = "#984EA3", TAM_M2 = "#4DAF4A",
  TAM_other = "#B2DF8A", DC_M1 = "#FFFF33",
  Other_immune = "#FB9A99", Unclassified = "#EEEEEE"
)

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
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_comp,
          file.path(FIG_DIR, "clustering", "06_composition_gating.pdf"),
          width = 12, height = 7)

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
         fontsize_row = 11, fontsize_col = 8,
         cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_method = "ward.D2",
         annotation_row = annot_row,
         annotation_colors = annot_colors,
         main = "Patient x cell type fractions - gating",
         filename = file.path(FIG_DIR, "clustering", "06_composition_heatmap_gating.pdf"),
         width = 14, height = 6)

# ── 7. ADD GATING LABELS TO SEURAT ───────────────────────────────────────────
seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = setNames(mfi_df$celltype_gating, mfi_df$cell_id),
  col.name = "celltype_gating"
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
          width = 18, height = 7)

# ── 8. SAVE ──────────────────────────────────────────────────────────────────
saveRDS(seurat_obj, file.path(RDS_DIR, "06_seurat_gated.rds"))
saveRDS(mfi_df,     file.path(RDS_DIR, "06_mfi_gated.rds"))
message("\nStep 6 complete.")
