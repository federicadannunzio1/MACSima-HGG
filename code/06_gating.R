# =============================================================================
# 06_gating.R — Rule-based cell type gating for MACSima imaging data
#
# Gating strategy based on reference diagram (info/gating strategy for the
# analysis.png). Per-patient GMM thresholds for each gating marker.
# Hierarchical decision tree assigns one cell type per cell.
#
# Adaptations from reference:
#   - OLIG2 not in panel → cannot distinguish oligodendrocytes/astrocytes.
#     Tumor identified by GFAP+|BCAN+ (BCAN is brain-specific, good proxy).
#   - Iba1 not in shared panel (9/15 ROI) → Microglia identified by P2Y12+
#     (highly specific for CNS-resident microglia).
#   - CD2 not in panel → omitted (measured marker in reference, not a gate).
#   - T cell activation (HLA-DR+) and proliferation (Ki-67+) stored as
#     metadata columns, not separate cell types (they overlay on
#     differentiation states, not exclusive categories).
#
# Output: cell type per cell + fractions per patient
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

# Gating markers — mapped from logical name to Seurat feature name
# Seurat replaces "_" with "-" in feature names: use the Seurat spelling.
gating_name_map <- c(
  # Immune vs non-immune
  CD45     = "CD45",
  # B cells
  CD20     = "CD20-Cytoplasmic",
  # T cell lineage
  CD3      = "CD3",
  CD56     = "CD56",
  CD4      = "CD4",
  CD8a     = "CD8a",
  FoxP3    = "FoxP3",
  # T cell exhaustion (per reference: PD1, TIM3, LAG-3)
  PD1      = "CD279",
  TIM3     = "TIM3",
  LAG3     = "CD223",
  # T cell differentiation
  CCR7     = "CCR7",
  CD45RA   = "CD45RA",
  CD45RO   = "CD45RO",
  # T cell activation/proliferation (metadata, not exclusive gate)
  HLA_DR   = "HLA-DR",
  Ki67     = "Ki-67",
  # NK subtypes
  CD16     = "CD16",
  # Myeloid
  CD15     = "CD15",
  CD14     = "CD14",
  CD68     = "CD68",
  P2Y12    = "P2Y12",
  CD163    = "CD163",
  CD206    = "CD206",
  CD80     = "CD80",
  # Non-immune
  CD31     = "CD31",
  CD34     = "CD34",
  GFAP     = "GFAP",
  BCAN     = "BCAN",
  Podoplanin = "Podoplanin",
  Vimentin = "Vimentin",
  CD105    = "CD105"
)

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
# Reference: info/gating strategy for the analysis.png
# Adaptations documented at top of file and in deviations table below.

gate_cell_vectorized <- function() {
  n <- nrow(mfi_df)

  # Compute positivity for each gating marker (per-patient thresholds)
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
    flag_col <- paste0("flag_", lname)
    if (flag_col %in% colnames(flag_mat)) {
      is_flagged         <- !is.na(flag_mat[[flag_col]]) & flag_mat[[flag_col]]
      result[is_flagged] <- NA
    }
    pos[[lname]] <- result
  }

  # Helper: treat NA as negative (conservative)
  p <- function(marker) { ifelse(is.na(pos[[marker]]), FALSE, pos[[marker]]) }

  celltype <- rep("Unclassified", n)

  # ═══════════════════════════════════════════════════════════════════════════
  # IMMUNE BRANCH (CD45+)
  # ═══════════════════════════════════════════════════════════════════════════
  is_immune <- if ("CD45" %in% names(pos)) p("CD45") else rep(FALSE, n)

  # Fallback if too many CD45 NA (QC issue)
  n_cd45_na <- sum(is.na(pos[["CD45"]]))
  if (n_cd45_na > 0.3 * n) {
    message(sprintf("  WARNING: CD45 NA in %.0f%% of cells - using CD3/CD68/CD56 as immune proxy",
                    100 * n_cd45_na / n))
    is_immune <- p("CD3") | p("CD68") | p("CD56")
  }

  # ── B cells (CD45+ → CD20+) ──
  is_B <- is_immune & p("CD20")
  celltype[is_B] <- "B_cell"

  # ── T cells & NK (CD45+ → CD20-) ──
  is_cd20neg <- is_immune & !p("CD20")

  # NKT cells: CD3+ & CD56+ (identified before T/NK split per reference)
  is_NKT <- is_cd20neg & p("CD3") & p("CD56")
  celltype[is_NKT] <- "NKT"

  # T cells: CD3+ & CD56- (exclude NKT)
  is_T <- is_cd20neg & p("CD3") & !p("CD56")

  # CD8+ T cells (CD8a+ & CD4-)
  is_CD8 <- is_T & p("CD8a") & !p("CD4")
  # CD4+ T cells (CD4+ & CD8a-)
  is_CD4 <- is_T & p("CD4") & !p("CD8a")
  # Double-negative T cells
  is_Tdn <- is_T & !p("CD8a") & !p("CD4")
  celltype[is_Tdn] <- "T_other"

  # ── CD4+ T cell subtypes ──
  # Priority 1: Treg (FoxP3+)
  is_Treg <- is_CD4 & p("FoxP3")
  celltype[is_Treg] <- "T_reg"

  # Priority 2: Exhausted (PD1+ | TIM3+ | LAG3+), not Treg
  is_CD4_exh <- is_CD4 & !p("FoxP3") & (p("PD1") | p("TIM3") | p("LAG3"))
  celltype[is_CD4_exh] <- "T_CD4_exhausted"

  # Priority 3: Differentiation states for remaining CD4+ (not Treg, not exhausted)
  is_CD4_rest <- is_CD4 & !p("FoxP3") & !(p("PD1") | p("TIM3") | p("LAG3"))

  celltype[is_CD4_rest & p("CCR7") & p("CD45RA")]   <- "T_CD4_naive"
  celltype[is_CD4_rest & p("CCR7") & p("CD45RO")]    <- "T_CD4_CM"
  celltype[is_CD4_rest & !p("CCR7") & p("CD45RO")]   <- "T_CD4_EM"
  celltype[is_CD4_rest & !p("CCR7") & p("CD45RA")]   <- "T_CD4_TEMRA"
  # Remaining CD4+ without clear differentiation
  celltype[is_CD4_rest & celltype == "Unclassified"]  <- "T_CD4_other"

  # ── CD8+ T cell subtypes ──
  # Priority 1: Exhausted
  is_CD8_exh <- is_CD8 & (p("PD1") | p("TIM3") | p("LAG3"))
  celltype[is_CD8_exh] <- "T_CD8_exhausted"

  # Priority 2: Differentiation states for remaining CD8+
  is_CD8_rest <- is_CD8 & !(p("PD1") | p("TIM3") | p("LAG3"))

  celltype[is_CD8_rest & p("CCR7") & p("CD45RA")]   <- "T_CD8_naive"
  celltype[is_CD8_rest & p("CCR7") & p("CD45RO")]    <- "T_CD8_CM"
  celltype[is_CD8_rest & !p("CCR7") & p("CD45RO")]   <- "T_CD8_EM"
  celltype[is_CD8_rest & !p("CCR7") & p("CD45RA")]   <- "T_CD8_TEMRA"
  celltype[is_CD8_rest & celltype == "Unclassified"]  <- "T_CD8_other"

  # ── NK cells (CD3- & CD56+) ──
  is_NK <- is_cd20neg & !p("CD3") & p("CD56")
  celltype[is_NK & p("CD16")]  <- "NK_effector"
  celltype[is_NK & !p("CD16")] <- "NK_immature"

  # ── Myeloid & other immune (CD3- & CD56-) ──
  is_myeloid_pool <- is_cd20neg & !p("CD3") & !p("CD56")

  # Neutrophils (CD15+)
  is_neutro <- is_myeloid_pool & p("CD15")
  celltype[is_neutro] <- "Neutrophil"

  # Non-neutrophil myeloid (CD15-)
  is_cd15neg <- is_myeloid_pool & !p("CD15")

  # Microglia: CD14- & P2Y12+ (adapted from reference: CD14-/CD15-/Iba1+ → P2Y12+)
  is_microglia <- is_cd15neg & !p("CD14") & p("P2Y12")
  celltype[is_microglia] <- "Microglia"

  # CD14+ myeloid cells
  is_cd14pos <- is_cd15neg & p("CD14")
  # CD14+ & CD68+ → macrophages (reference: Iba1+/CD68+ → adapted to CD68+)
  is_macro <- is_cd14pos & p("CD68")
  celltype[is_macro & (p("CD163") | p("CD206"))]  <- "TAM_M2"
  celltype[is_macro & !p("CD163") & !p("CD206")]  <- "TAM_M1"
  # CD14+ & CD68- → monocytes
  celltype[is_cd14pos & !p("CD68")] <- "Monocyte"

  # CD14- & P2Y12- → check for CD68+ (non-classical myeloid)
  is_cd14neg_p2y12neg <- is_cd15neg & !p("CD14") & !p("P2Y12")
  is_other_macro <- is_cd14neg_p2y12neg & p("CD68")
  celltype[is_other_macro & (p("CD163") | p("CD206"))] <- "TAM_M2"
  celltype[is_other_macro & !p("CD163") & !p("CD206") & p("CD80")] <- "TAM_M1"
  celltype[is_other_macro & !p("CD163") & !p("CD206") & !p("CD80")] <- "Macrophage_other"

  # Remaining CD45+ without lineage markers
  celltype[is_cd14neg_p2y12neg & !p("CD68") & celltype == "Unclassified"] <- "Other_immune"

  # ═══════════════════════════════════════════════════════════════════════════
  # NON-IMMUNE BRANCH (CD45-)
  # ═══════════════════════════════════════════════════════════════════════════
  is_nonimmune <- !is_immune

  # Endothelial: CD31+ | CD34+ (reference: CD31+/CD34+)
  is_endo <- is_nonimmune & (p("CD31") | p("CD34"))
  celltype[is_endo] <- "Endothelial"

  # Tumor: GFAP+ | BCAN+ (adapted: reference uses GFAP+/OLIG2+, OLIG2 not available)
  is_tumor <- is_nonimmune & !is_endo & (p("GFAP") | p("BCAN"))
  celltype[is_tumor] <- "Tumor"

  # Stromal subtypes (reference: Podoplanin+/Vimentin+ → mesenchymal/fibroblast)
  is_stromal_pool <- is_nonimmune & !is_endo & !is_tumor

  is_mesench <- is_stromal_pool & p("Podoplanin") & p("Vimentin") & p("CD105")
  celltype[is_mesench] <- "Mesenchymal"

  is_fibro <- is_stromal_pool & p("Podoplanin") & p("Vimentin") & !p("CD105")
  celltype[is_fibro] <- "Fibroblast"

  is_stromal <- is_stromal_pool & p("Vimentin") & !p("Podoplanin")
  celltype[is_stromal] <- "Stromal"

  celltype[is_stromal_pool & !p("Vimentin") & celltype == "Unclassified"] <- "Unclassified_nonimmune"

  celltype
}

message("Applying gating tree...")
mfi_df$celltype_gating <- gate_cell_vectorized()

# ── 4b. T CELL FUNCTIONAL ANNOTATIONS (metadata, not exclusive types) ────────
# Per reference: HLA-DR+ = activated, Ki-67+ = proliferating
# These overlay on differentiation states and are stored as separate columns.
message("Adding T cell functional annotations...")

is_any_T <- grepl("^T_|^T_reg$|^NKT$", mfi_df$celltype_gating)
mfi_df$T_activated     <- NA
mfi_df$T_proliferating <- NA

for (pat in patients) {
  idx_pat <- mfi_df$patient_id == pat & is_any_T

  hla_thresh <- thresholds_per_patient[[pat]][["HLA_DR"]]$threshold
  ki67_thresh <- thresholds_per_patient[[pat]][["Ki67"]]$threshold

  hla_name  <- gating_name_map[["HLA_DR"]]
  ki67_name <- gating_name_map[["Ki67"]]

  mfi_df$T_activated[idx_pat]     <- mfi_df[[hla_name]][idx_pat] > hla_thresh
  mfi_df$T_proliferating[idx_pat] <- mfi_df[[ki67_name]][idx_pat] > ki67_thresh
}

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
CT_ORDER_GATING <- c(
  # Non-immune
  "Tumor", "Endothelial", "Mesenchymal", "Fibroblast", "Stromal",
  "Unclassified_nonimmune",
  # T cells
  "T_CD4_naive", "T_CD4_CM", "T_CD4_EM", "T_CD4_TEMRA",
  "T_CD4_exhausted", "T_CD4_other",
  "T_reg",
  "T_CD8_naive", "T_CD8_CM", "T_CD8_EM", "T_CD8_TEMRA",
  "T_CD8_exhausted", "T_CD8_other",
  "T_other",
  # B, NK, NKT
  "B_cell", "NKT",
  "NK_effector", "NK_immature",
  # Myeloid
  "Neutrophil", "Microglia", "Monocyte",
  "TAM_M1", "TAM_M2", "Macrophage_other",
  "Other_immune",
  # Residual
  "Unclassified"
)

CT_COLORS_GATING <- c(
  # Non-immune
  Tumor = "#E41A1C", Endothelial = "#A65628", Mesenchymal = "#FDBF6F",
  Fibroblast = "#CAB2D6", Stromal = "#D2B48C",
  Unclassified_nonimmune = "#CCCCCC",
  # CD4 T cells (blue scale)
  T_CD4_naive = "#DEEBF7", T_CD4_CM = "#9ECAE1", T_CD4_EM = "#4292C6",
  T_CD4_TEMRA = "#2171B5", T_CD4_exhausted = "#08306B", T_CD4_other = "#6BAED6",
  T_reg = "#9E9AC8",
  # CD8 T cells (green scale)
  T_CD8_naive = "#E5F5E0", T_CD8_CM = "#A1D99B", T_CD8_EM = "#41AB5D",
  T_CD8_TEMRA = "#238B45", T_CD8_exhausted = "#00441B", T_CD8_other = "#74C476",
  T_other = "#C6DBEF",
  # B, NK, NKT
  B_cell = "#FB9A99", NKT = "#E78AC3",
  NK_effector = "#FF7F00", NK_immature = "#FDBF6F",
  # Myeloid
  Neutrophil = "#F0E442", Microglia = "#984EA3", Monocyte = "#BC80BD",
  TAM_M1 = "#FFD92F", TAM_M2 = "#4DAF4A", Macrophage_other = "#B2DF8A",
  Other_immune = "#999999",
  # Residual
  Unclassified = "#EEEEEE"
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
# Document where implementation differs from reference diagram
deviations <- data.frame(
  cell_type = c("Microglia", "Tumor", "Oligodendrocyte", "Astrocyte",
                "Tumor_microglia", "T_exhausted", "Myeloid_entry"),
  reference = c(
    "CD14-/CD15- → Iba1+ → microglia",
    "GFAP+/OLIG2+ → tumor",
    "GFAP-/OLIG2+ → oligodendrocyte",
    "Podoplanin-/Vimentin- → GFAP+/OLIG2- → astrocyte",
    "Non-immune → Iba1+/CD68+ → tumor microglia",
    "PD1+, TIM3+, CD223 (CD4 & CD8)",
    "CD14 as primary myeloid gate"
  ),
  implementation = c(
    "CD14-/CD15- → P2Y12+ → microglia",
    "GFAP+ | BCAN+ → tumor",
    "Not implemented",
    "Not implemented",
    "Not implemented",
    "PD1+ | TIM3+ | LAG3+ (CD4 & CD8)",
    "CD14 as primary myeloid gate (same as reference)"
  ),
  reason = c(
    "Iba1 not in shared panel (9/15 ROI). P2Y12 is highly specific for CNS-resident microglia.",
    "OLIG2 not in panel. BCAN (brevican) is brain-specific, good tumor/glial proxy.",
    "OLIG2 not in panel. Cannot distinguish oligodendrocytes.",
    "OLIG2 not in panel. Cannot distinguish astrocytes from GFAP+ tumor cells.",
    "Iba1 not in shared panel. Tumor-associated microglia not separately identified.",
    "Identical markers. TIGIT available but not included per reference strategy.",
    "Consistent with reference."
  ),
  stringsAsFactors = FALSE
)
write_tsv(deviations, file.path(TAB_DIR, "06_gating_deviations.tsv"))

message("\nStep 6 complete.")
