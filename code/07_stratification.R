# =============================================================================
# 07_stratification.R
# Stratificazione biologica per paziente — BLIND all'outcome clinico
#
# Input: composizione TME per paziente (da gating Step 6 e clustering Step 5b)
# Metodo: clustering gerarchico + PCA sui profili di cell type per paziente
# Output: dendrogram, biplot PCA, heatmap paziente x feature
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
library(ggrepel)

message("=== Step 7: Patient stratification (blind) ===")

seurat_obj <- readRDS(file.path(RDS_DIR, "06_seurat_gated.rds"))

# ── 1. Feature matrix paziente: composizione gating ──────────────────────────
comp_gating <- read_tsv(file.path(TAB_DIR, "06_composition_gating.tsv"),
                         show_col_types = FALSE) %>%
  mutate(celltype_gating = as.character(celltype_gating))

comp_cluster <- read_tsv(file.path(TAB_DIR, "05b_composition_annotated.tsv"),
                          show_col_types = FALSE) %>%
  mutate(celltype = as.character(celltype))

# Matrice gating: righe = pazienti, colonne = cell types
feat_gating <- comp_gating %>%
  select(patient_id, celltype_gating, fraction) %>%
  pivot_wider(names_from = celltype_gating, values_from = fraction, values_fill = 0) %>%
  tibble::column_to_rownames("patient_id")

# Matrice clustering: righe = pazienti, colonne = cell types (annotati)
feat_cluster <- comp_cluster %>%
  select(patient_id, celltype, fraction) %>%
  pivot_wider(names_from = celltype, values_from = fraction, values_fill = 0) %>%
  tibble::column_to_rownames("patient_id")

message(sprintf("Feature matrix gating: %d pazienti x %d cell types",
                nrow(feat_gating), ncol(feat_gating)))
message(sprintf("Feature matrix clustering: %d pazienti x %d cell types",
                nrow(feat_cluster), ncol(feat_cluster)))

# ── 2. PCA su composizione gating ────────────────────────────────────────────
# Con solo 5 pazienti la PCA è indicativa, non conclusiva
feat_scaled <- scale(feat_gating)
# Rimuovi colonne con varianza 0
feat_scaled <- feat_scaled[, apply(feat_scaled, 2, function(x) !any(is.na(x))), drop = FALSE]

pca_res <- prcomp(feat_scaled, center = FALSE, scale. = FALSE)
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x[, 1:min(3, ncol(pca_res$x))]) %>%
  tibble::rownames_to_column("patient_id")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = patient_id, label = patient_id)) +
  geom_point(size = 6) +
  geom_text_repel(size = 5, fontface = "bold", box.padding = 0.5) +
  scale_color_manual(values = PATIENT_COLORS) +
  labs(
    title = "PCA — Profilo TME per paziente (gating, blind)",
    x = sprintf("PC1 (%s%% var)", var_exp[1]),
    y = sprintf("PC2 (%s%% var)", var_exp[2]),
    color = "Patient"
  ) +
  theme_hgg() +
  theme(legend.position = "none")

save_plot(p_pca,
          file.path(FIG_DIR, "stratification", "07_pca_patients_gating.pdf"),
          width = 7, height = 6)

# Loading plot: quali cell types guidano la separazione
loadings_df <- as.data.frame(pca_res$rotation[, 1:2]) %>%
  tibble::rownames_to_column("celltype") %>%
  mutate(length = sqrt(PC1^2 + PC2^2)) %>%
  arrange(desc(length))

p_load <- ggplot(loadings_df, aes(x = reorder(celltype, PC1), y = PC1, fill = PC1 > 0)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"), guide = "none") +
  coord_flip() +
  labs(title = "PC1 loadings — cell types che separano i pazienti",
       x = "", y = "Loading PC1") +
  theme_hgg()

save_plot(p_load,
          file.path(FIG_DIR, "stratification", "07_pca_loadings.pdf"),
          width = 7, height = 6)

# ── 3. Clustering gerarchico ─────────────────────────────────────────────────
dist_mat   <- dist(feat_scaled, method = "euclidean")
hclust_res <- hclust(dist_mat, method = "ward.D2")

pdf(file.path(FIG_DIR, "stratification", "07_dendrogram_patients.pdf"),
    width = 7, height = 5)
plot(hclust_res,
     main = "Clustering gerarchico pazienti — profilo TME gating (blind)",
     xlab = "Patient", ylab = "Ward distance",
     sub  = "", hang = -1)
dev.off()
message("Dendrogram saved.")

# ── 4. Heatmap paziente x cell type (gating) ─────────────────────────────────
CT_ORDER_GATING <- c(
  "Tumor", "Endothelial", "Stromal", "Unclassified_nonimmune",
  "T_CD8_exhausted", "T_CD8_effector", "T_CD4_helper", "T_reg", "T_other",
  "NK",
  "Microglia", "TAM_M2", "TAM_other", "DC_M1",
  "Other_immune", "Unclassified"
)
feat_ordered <- feat_gating[, intersect(CT_ORDER_GATING, colnames(feat_gating)), drop = FALSE]

pheatmap(as.matrix(feat_ordered),
         color             = colorRampPalette(c("white", "#08519C"))(100),
         display_numbers   = TRUE,
         number_format     = "%.2f",
         fontsize_row      = 12,
         fontsize_col      = 9,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         clustering_method = "ward.D2",
         main = "Profilo TME per paziente — gating (blind to outcome)",
         filename = file.path(FIG_DIR, "stratification", "07_heatmap_patients_gating.pdf"),
         width = 14, height = 5)

# Heatmap clustering annotation (secondo approccio, convergenza dei metodi)
feat_cl_ordered <- feat_cluster[, intersect(
  c("Tumor", "Tumor_G34", "Tumor_inflamed",
    "T_CD8_exhausted", "T_CD4",
    "TAM_M2", "TAM_Microglia", "Microglia", "DC_activated", "Endothelial"),
  colnames(feat_cluster)), drop = FALSE]

pheatmap(as.matrix(feat_cl_ordered),
         color             = colorRampPalette(c("white", "#08519C"))(100),
         display_numbers   = TRUE,
         number_format     = "%.2f",
         fontsize_row      = 12,
         fontsize_col      = 9,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         clustering_method = "ward.D2",
         main = "Profilo TME per paziente — clustering manuale (blind to outcome)",
         filename = file.path(FIG_DIR, "stratification", "07_heatmap_patients_clustering.pdf"),
         width = 13, height = 5)

# ── 5. Immune vs Tumor ratio per paziente (metrica riassuntiva chiave) ────────
immune_ct_gating <- c(
  "T_CD8_exhausted", "T_CD8_effector", "T_CD4_helper", "T_reg", "T_other",
  "NK", "Microglia", "TAM_M2", "TAM_other", "DC_M1", "Other_immune"
)

immune_summary <- comp_gating %>%
  mutate(compartment = case_when(
    celltype_gating %in% immune_ct_gating    ~ "Immune",
    celltype_gating == "Tumor"               ~ "Tumor",
    celltype_gating == "Endothelial"         ~ "Endothelial",
    celltype_gating == "Stromal"             ~ "Stromal",
    TRUE                                      ~ "Unclassified"
  )) %>%
  group_by(patient_id, compartment) %>%
  summarise(fraction = sum(fraction), .groups = "drop")

pat_order <- immune_summary %>%
  filter(compartment == "Immune") %>%
  arrange(desc(fraction)) %>%
  pull(patient_id)

p_immune_ratio <- ggplot(immune_summary,
                          aes(x = factor(patient_id, levels = pat_order),
                              y = fraction, fill = compartment)) +
  geom_col(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    Immune       = "#377EB8",
    Tumor        = "#E41A1C",
    Endothelial  = "#A65628",
    Stromal      = "#FDBF6F",
    Unclassified = "#CCCCCC"
  )) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Patient", y = "Fraction", fill = "Compartment",
       title = "Immune vs Tumor ratio per paziente (gating, blind)") +
  theme_hgg()

save_plot(p_immune_ratio,
          file.path(FIG_DIR, "stratification", "07_immune_tumor_ratio.pdf"),
          width = 8, height = 5)

# Salva tabella riassuntiva
immune_ratio_wide <- immune_summary %>%
  pivot_wider(names_from = compartment, values_from = fraction, values_fill = 0) %>%
  mutate(immune_tumor_ratio = Immune / (Tumor + 0.001))

write_tsv(immune_ratio_wide, file.path(TAB_DIR, "07_immune_tumor_ratio.tsv"))
message("\nImmune/Tumor ratio per paziente:")
print(as.data.frame(immune_ratio_wide[, c("patient_id", "Immune", "Tumor", "immune_tumor_ratio")]))

message("\nStep 7 complete.")
