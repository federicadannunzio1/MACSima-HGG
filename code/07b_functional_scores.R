# =============================================================================
# 07b_functional_scores.R
# Score funzionali per paziente — BLIND all'outcome
#
# Per ogni paziente, calcoliamo score di:
#   1. Esaurimento T_CD8 (TIGIT, PD-1, TIM-3)
#   2. Polarizzazione mieloide M1/M2 (CD80 vs CD163)
#   3. Immunogenicità tumorale (HLA-ABC, HLA-DR su cellule tumorali)
#   4. Checkpoint ligand su tumore (B7H3, PD-L1/CD274)
#   5. Proliferazione (Ki-67) per compartimento
#   6. Carico NKG2D/GD2 su tumore (target per immunoterapia)
#
# Tutti gli score sono mediati per paziente (non per cellula), perché
# il paziente è l'unità biologica e clinica di analisi.
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

message("=== Step 7b: Functional immune scores (blind) ===")

seurat_obj <- readRDS(file.path(RDS_DIR, "06_seurat_gated.rds"))

# Estrai dati normalizzati + metadata
mfi_norm <- t(as.matrix(GetAssayData(seurat_obj, assay = "MFI", layer = "data")))
mfi_df   <- as.data.frame(mfi_norm)
mfi_df$cell_id        <- rownames(mfi_df)
mfi_df$patient_id     <- seurat_obj$patient_id
mfi_df$sample_id      <- seurat_obj$sample_id
mfi_df$celltype_gating <- seurat_obj$celltype_gating
mfi_df$celltype        <- seurat_obj$celltype  # annotazione manuale da clustering

markers_available <- colnames(mfi_norm)
message(sprintf("Marker disponibili: %d", length(markers_available)))

# Helper: media sicura (NA se nessuna cellula del subset)
safe_mean <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# Helper: score = media di più marker, per paziente × subset
score_by_patient <- function(subset_df, markers) {
  avail <- intersect(markers, colnames(subset_df))
  if (length(avail) == 0) return(tibble(patient_id = character(), score = numeric()))
  subset_df %>%
    group_by(patient_id) %>%
    summarise(score = safe_mean(rowMeans(across(all_of(avail)), na.rm = TRUE)),
              n_cells = n(), .groups = "drop")
}

# ── 1. T_CD8 exhaustion score ────────────────────────────────────────────────
# Marker di esaurimento: TIGIT, PD-1 (CD279), TIM-3, CTLA-4
EXHAUSTION_MARKERS <- c("TIGIT", "CD279", "TIM3", "CTLA-4")
avail_exhaust <- intersect(EXHAUSTION_MARKERS, markers_available)
message(sprintf("Exhaustion markers disponibili: %s", paste(avail_exhaust, collapse = ", ")))

t_cd8_cells <- mfi_df %>%
  filter(celltype_gating %in% c("T_CD8_exhausted", "T_CD8_effector"))

exhaust_score <- score_by_patient(t_cd8_cells, avail_exhaust) %>%
  rename(T_CD8_exhaustion = score, n_TCD8 = n_cells)

# ── 2. Mieloide M2 score (CD163 / (CD80 + CD163)) ────────────────────────────
M1_MARKERS <- c("CD80")
M2_MARKERS <- c("CD163", "CD206")
avail_m1 <- intersect(M1_MARKERS, markers_available)
avail_m2 <- intersect(M2_MARKERS, markers_available)

myeloid_cells <- mfi_df %>%
  filter(celltype_gating %in% c("Microglia", "TAM_M2", "TAM_other", "DC_M1", "Other_immune"))

m1_score <- score_by_patient(myeloid_cells, avail_m1) %>%
  rename(M1_score = score, n_myeloid = n_cells)
m2_score <- score_by_patient(myeloid_cells, avail_m2) %>%
  select(patient_id, M2_score = score)

myeloid_score <- m1_score %>%
  left_join(m2_score, by = "patient_id") %>%
  mutate(M2_M1_ratio = M2_score / (M1_score + 0.01))

# ── 3. Immunogenicità tumorale ────────────────────────────────────────────────
# HLA-ABC = MHC-I (presentazione a T_CD8), HLA-DR = MHC-II (presentazione a T_CD4)
IMMUNO_MARKERS <- c("HLA-ABC", "HLA-DR-DP-DQ", "HLA-DR")
# Seurat può chiamarlo HLA-DR-DP-DQ o HLA-DR...
avail_hla <- intersect(IMMUNO_MARKERS, markers_available)
# Prendi quello disponibile
avail_hla_simple <- avail_hla[1:min(2, length(avail_hla))]

tumor_cells <- mfi_df %>%
  filter(celltype_gating == "Tumor" | startsWith(as.character(celltype), "Tumor"))

immunogenicity_score <- score_by_patient(tumor_cells, avail_hla_simple) %>%
  rename(Tumor_immunogenicity = score, n_tumor = n_cells)

message(sprintf("HLA markers usati per immunogenicita': %s",
                paste(avail_hla_simple, collapse = ", ")))

# ── 4. Checkpoint ligand su tumore (immune evasion) ──────────────────────────
CHECKPOINT_LIG <- c("CD274", "B7H3")
avail_ckpt <- intersect(CHECKPOINT_LIG, markers_available)
message(sprintf("Checkpoint ligand markers: %s", paste(avail_ckpt, collapse = ", ")))

immune_evasion_score <- score_by_patient(tumor_cells, avail_ckpt) %>%
  rename(Tumor_immune_evasion = score)

# ── 5. Proliferazione per compartimento ──────────────────────────────────────
KI67_MARKER <- intersect(c("Ki-67", "Ki67", "KI-67"), markers_available)

if (length(KI67_MARKER) > 0) {
  ki67_tumor   <- score_by_patient(tumor_cells, KI67_MARKER) %>%
    rename(Ki67_tumor = score) %>% select(patient_id, Ki67_tumor)
  ki67_tcells  <- score_by_patient(t_cd8_cells, KI67_MARKER) %>%
    rename(Ki67_Tcell = score) %>% select(patient_id, Ki67_Tcell)
  message(sprintf("Ki-67 marker usato: %s", KI67_MARKER))
} else {
  message("Ki-67 non trovato — salto")
  ki67_tumor  <- tibble(patient_id = unique(mfi_df$patient_id), Ki67_tumor = NA_real_)
  ki67_tcells <- tibble(patient_id = unique(mfi_df$patient_id), Ki67_Tcell = NA_real_)
}

# ── 6. Antigeni tumorali targetabili: NKG2D, GD2, MGMT, EZHIP ──────────────
TUMOR_ANTIGEN_MARKERS <- c("NKG2D", "GD2", "MGMT", "EZHIP")
avail_ag <- intersect(TUMOR_ANTIGEN_MARKERS, markers_available)
message(sprintf("Tumor antigen markers: %s", paste(avail_ag, collapse = ", ")))

tumor_antigen_scores <- lapply(avail_ag, function(m) {
  score_by_patient(tumor_cells, m) %>%
    rename(!!m := score) %>% select(patient_id, !!m)
})
if (length(tumor_antigen_scores) > 0) {
  tumor_antigens <- Reduce(function(a, b) left_join(a, b, by = "patient_id"),
                            tumor_antigen_scores)
} else {
  tumor_antigens <- tibble(patient_id = unique(mfi_df$patient_id))
}

# ── 7. Immune infiltration density ───────────────────────────────────────────
immune_fraction <- mfi_df %>%
  mutate(is_immune = celltype_gating %in% c(
    "T_CD8_exhausted", "T_CD8_effector", "T_CD4_helper", "T_reg", "T_other",
    "NK", "Microglia", "TAM_M2", "TAM_other", "DC_M1", "Other_immune"
  )) %>%
  group_by(patient_id) %>%
  summarise(Immune_fraction = mean(is_immune, na.rm = TRUE), .groups = "drop")

# ── 8. Combina tutti gli score ────────────────────────────────────────────────
all_patients <- unique(mfi_df$patient_id)
base_df      <- tibble(patient_id = all_patients)

score_table <- base_df %>%
  left_join(exhaust_score     %>% select(patient_id, T_CD8_exhaustion, n_TCD8), by = "patient_id") %>%
  left_join(myeloid_score     %>% select(patient_id, M1_score, M2_score, M2_M1_ratio, n_myeloid), by = "patient_id") %>%
  left_join(immunogenicity_score %>% select(patient_id, Tumor_immunogenicity, n_tumor), by = "patient_id") %>%
  left_join(immune_evasion_score %>% select(patient_id, Tumor_immune_evasion), by = "patient_id") %>%
  left_join(ki67_tumor,  by = "patient_id") %>%
  left_join(ki67_tcells, by = "patient_id") %>%
  left_join(tumor_antigens, by = "patient_id") %>%
  left_join(immune_fraction, by = "patient_id")

write_tsv(score_table, file.path(TAB_DIR, "07b_functional_scores.tsv"))

message("\n--- Functional scores per paziente ---")
score_cols <- setdiff(colnames(score_table), c("patient_id", "n_TCD8", "n_myeloid", "n_tumor"))
for (col in score_cols) {
  vals <- score_table[[col]]
  if (all(is.na(vals))) next
  cat(sprintf("  %-25s: %s\n", col,
              paste(sprintf("%s=%.3f", score_table$patient_id, vals), collapse="  ")))
}

# ── 9. Heatmap degli score ────────────────────────────────────────────────────
score_mat <- score_table %>%
  select(patient_id, all_of(setdiff(score_cols, c("n_TCD8", "n_myeloid", "n_tumor", starts_with("n_"))))) %>%
  select(where(~ !all(is.na(.)))) %>%
  tibble::column_to_rownames("patient_id") %>%
  as.matrix()

# Rinomina per leggibilità
colnames(score_mat) <- gsub("_", " ", colnames(score_mat))

pheatmap(score_mat,
         scale             = "column",
         color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         display_numbers   = TRUE,
         number_format     = "%.2f",
         fontsize_row      = 12,
         fontsize_col      = 9,
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         clustering_method = "ward.D2",
         main = "Score funzionali per paziente (z-scored, blind to outcome)",
         filename = file.path(FIG_DIR, "stratification", "07b_functional_scores_heatmap.pdf"),
         width = 13, height = 5)

# ── 10. Barplot riassuntivo per marker ────────────────────────────────────────
score_long <- score_table %>%
  select(patient_id, T_CD8_exhaustion, M2_M1_ratio,
         Tumor_immunogenicity, Tumor_immune_evasion, Immune_fraction) %>%
  pivot_longer(-patient_id, names_to = "score", values_to = "value") %>%
  filter(!is.na(value))

score_labels <- c(
  T_CD8_exhaustion    = "T CD8\nexhaustion",
  M2_M1_ratio         = "Myeloid\nM2/M1 ratio",
  Tumor_immunogenicity = "Tumor\nimmunog. (HLA)",
  Tumor_immune_evasion = "Tumor\nimmune evasion",
  Immune_fraction     = "Immune\nfraction"
)

p_scores <- ggplot(score_long, aes(x = patient_id, y = value, fill = patient_id)) +
  geom_col(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = PATIENT_COLORS) +
  facet_wrap(~ score, scales = "free_y",
             labeller = labeller(score = score_labels)) +
  labs(x = "Patient", y = "Score (arcsinh units or ratio)",
       title = "Score funzionali immune per paziente (blind to outcome)") +
  theme_hgg() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_scores,
          file.path(FIG_DIR, "stratification", "07b_functional_scores_barplots.pdf"),
          width = 12, height = 8)

message("\nStep 7b complete.")
