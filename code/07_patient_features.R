# =============================================================================
# 07_patient_features.R
# Build patient-level feature matrix: cell type fractions + functional scores.
# This is the input for clinical integration (Step 08).
#
# Functional scores (per patient):
#   1. T_CD8 exhaustion (TIGIT, PD-1, TIM-3)
#   2. Myeloid M2/M1 polarization (CD163/CD80)
#   3. Tumor immunogenicity (HLA-ABC, HLA-DR on tumor cells)
#   4. Checkpoint ligand on tumor (B7H3, PD-L1/CD274)
#   5. Proliferation (Ki-67) per compartment
#   6. Tumor antigen load (NKG2D, GD2, MGMT, EZHIP)
#   7. Immune infiltration fraction
#
# All scores are patient-level means (patient = unit of analysis).
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))
source(file.path(getwd(), "code/utils/plot_utils.R"))

library(Seurat)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)

message("=== Step 7: Patient-level feature matrix ===")

seurat_obj <- readRDS(file.path(RDS_DIR, "06_seurat_gated.rds"))

mfi_norm <- t(as.matrix(GetAssayData(seurat_obj, assay = "MFI", layer = "data")))
mfi_df   <- as.data.frame(mfi_norm)
mfi_df$cell_id         <- rownames(mfi_df)
mfi_df$patient_id      <- seurat_obj$patient_id
mfi_df$sample_id       <- seurat_obj$sample_id
mfi_df$celltype_gating <- seurat_obj$celltype_gating

markers_available <- colnames(mfi_norm)

# ── Helpers ──────────────────────────────────────────────────────────────────
safe_mean <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

score_by_patient <- function(subset_df, markers) {
  avail <- intersect(markers, colnames(subset_df))
  if (length(avail) == 0) return(tibble(patient_id = character(), score = numeric()))
  subset_df %>%
    group_by(patient_id) %>%
    summarise(score = safe_mean(rowMeans(across(all_of(avail)), na.rm = TRUE)),
              n_cells = n(), .groups = "drop")
}

# ── 1. Cell type composition (from gating) ────────────────────────────────────
comp_gating <- read_tsv(file.path(TAB_DIR, "06_composition_gating.tsv"), show_col_types = FALSE)

feat_composition <- comp_gating %>%
  select(patient_id, celltype_gating, fraction) %>%
  mutate(celltype_gating = paste0("frac_", celltype_gating)) %>%
  pivot_wider(names_from = celltype_gating, values_from = fraction, values_fill = 0) %>%
  tibble::column_to_rownames("patient_id")

# ── 2. Functional scores ─────────────────────────────────────────────────────

# 2a. T_CD8 exhaustion
EXHAUST_MARKERS <- intersect(c("TIGIT", "CD279", "TIM3", "CTLA-4"), markers_available)
t_cd8_cells <- mfi_df %>% filter(celltype_gating %in% c("T_CD8_exhausted", "T_CD8_effector"))
exhaust_score <- score_by_patient(t_cd8_cells, EXHAUST_MARKERS) %>%
  rename(T_CD8_exhaustion = score, n_TCD8 = n_cells)

# 2b. Myeloid M2/M1
myeloid_cells <- mfi_df %>%
  filter(celltype_gating %in% c("Microglia", "TAM_M2", "TAM_other", "DC_M1", "Other_immune"))
m1_score <- score_by_patient(myeloid_cells, intersect(c("CD80"), markers_available)) %>%
  rename(M1_score = score, n_myeloid = n_cells)
m2_score <- score_by_patient(myeloid_cells, intersect(c("CD163", "CD206"), markers_available)) %>%
  select(patient_id, M2_score = score)
myeloid_score <- m1_score %>%
  left_join(m2_score, by = "patient_id") %>%
  mutate(M2_M1_ratio = M2_score / (M1_score + 0.01))

# 2c. Tumor immunogenicity (HLA on tumor cells)
tumor_cells <- mfi_df %>% filter(celltype_gating == "Tumor")
avail_hla <- intersect(c("HLA-ABC", "HLA-DR-DP-DQ", "HLA-DR"), markers_available)
immunogenicity_score <- score_by_patient(tumor_cells, avail_hla[1:min(2, length(avail_hla))]) %>%
  rename(Tumor_immunogenicity = score, n_tumor = n_cells)

# 2d. Checkpoint ligand on tumor
avail_ckpt <- intersect(c("CD274", "B7H3"), markers_available)
immune_evasion_score <- score_by_patient(tumor_cells, avail_ckpt) %>%
  rename(Tumor_immune_evasion = score)

# 2e. Proliferation (Ki-67)
KI67_MARKER <- intersect(c("Ki-67", "Ki67", "KI-67"), markers_available)
if (length(KI67_MARKER) > 0) {
  ki67_tumor  <- score_by_patient(tumor_cells, KI67_MARKER) %>%
    rename(Ki67_tumor = score) %>% select(patient_id, Ki67_tumor)
  ki67_tcells <- score_by_patient(t_cd8_cells, KI67_MARKER) %>%
    rename(Ki67_Tcell = score) %>% select(patient_id, Ki67_Tcell)
} else {
  all_pats <- unique(mfi_df$patient_id)
  ki67_tumor  <- tibble(patient_id = all_pats, Ki67_tumor = NA_real_)
  ki67_tcells <- tibble(patient_id = all_pats, Ki67_Tcell = NA_real_)
}

# 2f. Tumor antigens
avail_ag <- intersect(c("NKG2D", "GD2", "MGMT", "EZHIP"), markers_available)
tumor_antigen_scores <- lapply(avail_ag, function(m) {
  score_by_patient(tumor_cells, m) %>% rename(!!m := score) %>% select(patient_id, !!m)
})
if (length(tumor_antigen_scores) > 0) {
  tumor_antigens <- Reduce(function(a, b) left_join(a, b, by = "patient_id"), tumor_antigen_scores)
} else {
  tumor_antigens <- tibble(patient_id = unique(mfi_df$patient_id))
}

# 2g. Immune infiltration
immune_ct <- c("T_CD8_exhausted", "T_CD8_effector", "T_CD4_helper", "T_reg", "T_other",
               "NK", "Microglia", "TAM_M2", "TAM_other", "DC_M1", "Other_immune")
immune_fraction <- mfi_df %>%
  mutate(is_immune = celltype_gating %in% immune_ct) %>%
  group_by(patient_id) %>%
  summarise(Immune_fraction = mean(is_immune, na.rm = TRUE), .groups = "drop")

# ── 3. Combine all scores ────────────────────────────────────────────────────
all_patients <- unique(mfi_df$patient_id)
base_df      <- tibble(patient_id = all_patients)

score_table <- base_df %>%
  left_join(exhaust_score %>% select(patient_id, T_CD8_exhaustion, n_TCD8), by = "patient_id") %>%
  left_join(myeloid_score %>% select(patient_id, M1_score, M2_score, M2_M1_ratio, n_myeloid), by = "patient_id") %>%
  left_join(immunogenicity_score %>% select(patient_id, Tumor_immunogenicity, n_tumor), by = "patient_id") %>%
  left_join(immune_evasion_score %>% select(patient_id, Tumor_immune_evasion), by = "patient_id") %>%
  left_join(ki67_tumor,  by = "patient_id") %>%
  left_join(ki67_tcells, by = "patient_id") %>%
  left_join(tumor_antigens, by = "patient_id") %>%
  left_join(immune_fraction, by = "patient_id")

write_tsv(score_table, file.path(TAB_DIR, "07_functional_scores.tsv"))

# ── 4. Immune vs Tumor ratio ─────────────────────────────────────────────────
immune_summary <- comp_gating %>%
  mutate(compartment = case_when(
    celltype_gating %in% immune_ct ~ "Immune",
    celltype_gating == "Tumor"     ~ "Tumor",
    celltype_gating == "Endothelial" ~ "Endothelial",
    celltype_gating == "Stromal"   ~ "Stromal",
    TRUE                            ~ "Unclassified"
  )) %>%
  group_by(patient_id, compartment) %>%
  summarise(fraction = sum(fraction), .groups = "drop")

immune_ratio_wide <- immune_summary %>%
  pivot_wider(names_from = compartment, values_from = fraction, values_fill = 0) %>%
  mutate(immune_tumor_ratio = Immune / (Tumor + 0.001))

write_tsv(immune_ratio_wide, file.path(TAB_DIR, "07_immune_tumor_ratio.tsv"))

# ── 5. Combined feature matrix (composition + scores) ────────────────────────
# This is the key input for Step 08
feat_scores <- score_table %>%
  select(-starts_with("n_")) %>%
  tibble::column_to_rownames("patient_id")

feat_combined <- merge(feat_composition, feat_scores, by = "row.names", all = TRUE) %>%
  tibble::column_to_rownames("Row.names")

# Add immune/tumor ratio
ratio_df <- immune_ratio_wide %>%
  select(patient_id, immune_tumor_ratio) %>%
  tibble::column_to_rownames("patient_id")
feat_combined <- merge(feat_combined, ratio_df, by = "row.names", all = TRUE) %>%
  tibble::column_to_rownames("Row.names")

saveRDS(feat_combined, file.path(RDS_DIR, "07_patient_features.rds"))

message(sprintf("\nPatient feature matrix: %d patients x %d features",
                nrow(feat_combined), ncol(feat_combined)))
message("Features: ", paste(colnames(feat_combined), collapse = ", "))

# ── 6. Heatmap of functional scores ──────────────────────────────────────────
score_cols <- setdiff(colnames(feat_scores),
                      c("n_TCD8", "n_myeloid", "n_tumor"))
score_mat <- feat_scores %>%
  select(all_of(score_cols)) %>%
  select(where(~ !all(is.na(.)))) %>%
  as.matrix()
colnames(score_mat) <- gsub("_", " ", colnames(score_mat))

annot_row <- CLINICAL[, c("patient_id", "pfs_group", "pfs_months")] %>%
  tibble::column_to_rownames("patient_id")
annot_row <- annot_row[rownames(score_mat), , drop = FALSE]

pheatmap(score_mat,
         scale = "column",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_row = 12, fontsize_col = 9,
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_method = "ward.D2",
         annotation_row = annot_row,
         annotation_colors = list(pfs_group = PFS_GROUP_COLORS),
         main = "Functional scores per patient (z-scored)",
         filename = file.path(FIG_DIR, "stratification", "07_functional_scores_heatmap.pdf"),
         width = 13, height = 6)

# ── 7. Barplot summary ───────────────────────────────────────────────────────
score_long <- score_table %>%
  select(patient_id, T_CD8_exhaustion, M2_M1_ratio,
         Tumor_immunogenicity, Tumor_immune_evasion, Immune_fraction) %>%
  pivot_longer(-patient_id, names_to = "score", values_to = "value") %>%
  filter(!is.na(value)) %>%
  left_join(CLINICAL[, c("patient_id", "pfs_group")], by = "patient_id")

p_scores <- ggplot(score_long, aes(x = reorder(patient_id, -value),
                                    y = value, fill = pfs_group)) +
  geom_col(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = PFS_GROUP_COLORS) +
  facet_wrap(~ score, scales = "free_y") +
  labs(x = "Patient", y = "Score", fill = "PFS group",
       title = "Functional immune scores per patient") +
  theme_hgg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

save_plot(p_scores,
          file.path(FIG_DIR, "stratification", "07_functional_scores_barplots.pdf"),
          width = 14, height = 10)

# Immune/tumor ratio plot
p_ratio <- ggplot(immune_ratio_wide %>%
                    left_join(CLINICAL[, c("patient_id", "pfs_group", "pfs_months")], by = "patient_id"),
                  aes(x = reorder(patient_id, -pfs_months), y = immune_tumor_ratio, fill = pfs_group)) +
  geom_col() +
  scale_fill_manual(values = PFS_GROUP_COLORS) +
  labs(x = "Patient (ordered by PFS)", y = "Immune/Tumor ratio",
       fill = "PFS group", title = "Immune vs Tumor ratio") +
  theme_hgg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_ratio,
          file.path(FIG_DIR, "stratification", "07_immune_tumor_ratio.pdf"),
          width = 8, height = 5)

message("\nStep 7 complete.")
