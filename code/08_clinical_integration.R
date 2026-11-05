# =============================================================================
# 08_clinical_integration.R
# Clinical integration & patient stratification
#
# a. Integrated heatmap: patients x TME features + clinical annotations
# b. Patient PCA/dendrogram with PFS group overlay
# c. Group comparison: responder vs non-responder (PFS >= 6m vs < 6m)
# d. Spearman correlations: TME features vs PFS/OS (continuous)
# e. Patient profiles & candidate biomarkers
#
# NOTE: n=8, all results are descriptive/hypothesis-generating.
#       No formal p-values — report effect sizes and rho only.
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))
source(file.path(getwd(), "code/utils/plot_utils.R"))

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

message("=== Step 8: Clinical integration ===")

# ── Load data ────────────────────────────────────────────────────────────────
feat_combined  <- readRDS(file.path(RDS_DIR, "07_patient_features.rds"))
score_table    <- read_tsv(file.path(TAB_DIR, "07_functional_scores.tsv"), show_col_types = FALSE)
comp_gating    <- read_tsv(file.path(TAB_DIR, "06_composition_gating.tsv"), show_col_types = FALSE)
immune_ratio   <- read_tsv(file.path(TAB_DIR, "07_immune_tumor_ratio.tsv"), show_col_types = FALSE)

strat_dir <- file.path(FIG_DIR, "stratification")

# ── 1. INTEGRATED HEATMAP: patients x features + clinical annotations ────────
message("1. Integrated heatmap...")

# Select biologically meaningful features (drop near-zero-variance)
feat_mat <- as.matrix(feat_combined)
feat_var <- apply(feat_mat, 2, var, na.rm = TRUE)
feat_mat <- feat_mat[, feat_var > 1e-6 & !is.na(feat_var), drop = FALSE]

# Row annotation: clinical variables
annot_row <- CLINICAL %>%
  select(patient_id, pfs_group, os_months, pfs_months, mgmt_status, histology) %>%
  tibble::column_to_rownames("patient_id")
annot_row <- annot_row[rownames(feat_mat), , drop = FALSE]

annot_colors <- list(
  pfs_group   = PFS_GROUP_COLORS,
  mgmt_status = c(methyliert = "#2166AC", nicht_methyliert = "#B2182B", unknown = "#999999")
)

pheatmap(feat_mat,
         scale = "column",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         fontsize_row = 11, fontsize_col = 7,
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_method = "ward.D2",
         annotation_row = annot_row,
         annotation_colors = annot_colors,
         main = "Patient x TME features (z-scored) with clinical annotations",
         filename = file.path(strat_dir, "08_integrated_heatmap.pdf"),
         width = 18, height = 7)

# ── 2. PCA + DENDROGRAM ON TME FEATURES ─────────────────────────────────────
message("2. Patient PCA & dendrogram...")

feat_scaled <- scale(feat_mat)
# Remove columns with NA (from scaling zero-variance)
feat_scaled <- feat_scaled[, !apply(feat_scaled, 2, function(x) any(is.na(x))), drop = FALSE]

pca_res <- prcomp(feat_scaled, center = FALSE, scale. = FALSE)
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x[, 1:min(3, ncol(pca_res$x))]) %>%
  tibble::rownames_to_column("patient_id") %>%
  left_join(CLINICAL, by = "patient_id")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = pfs_group, label = patient_id)) +
  geom_point(size = 6) +
  geom_text_repel(size = 4, fontface = "bold", box.padding = 0.5) +
  scale_color_manual(values = PFS_GROUP_COLORS) +
  labs(
    title = "PCA - TME profile per patient (colored by PFS group)",
    x = sprintf("PC1 (%s%% var)", var_exp[1]),
    y = sprintf("PC2 (%s%% var)", var_exp[2]),
    color = "PFS group"
  ) +
  theme_hgg()

save_plot(p_pca,
          file.path(strat_dir, "08_pca_patients_pfs.pdf"),
          width = 8, height = 7)

# PCA loadings
loadings_df <- as.data.frame(pca_res$rotation[, 1:2]) %>%
  tibble::rownames_to_column("feature") %>%
  mutate(length = sqrt(PC1^2 + PC2^2)) %>%
  arrange(desc(length)) %>%
  head(20)

p_load <- ggplot(loadings_df, aes(x = reorder(feature, PC1), y = PC1, fill = PC1 > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"), guide = "none") +
  coord_flip() +
  labs(title = "Top 20 PC1 loadings - features driving patient separation",
       x = "", y = "Loading PC1") +
  theme_hgg()

save_plot(p_load,
          file.path(strat_dir, "08_pca_loadings.pdf"),
          width = 8, height = 7)

# Dendrogram
dist_mat   <- dist(feat_scaled, method = "euclidean")
hclust_res <- hclust(dist_mat, method = "ward.D2")

# Color dendrogram labels by PFS group
label_colors <- PFS_GROUP_COLORS[CLINICAL$pfs_group[match(hclust_res$labels, CLINICAL$patient_id)]]

pdf(file.path(strat_dir, "08_dendrogram_patients.pdf"), width = 8, height = 5)
plot(hclust_res,
     main = "Hierarchical clustering of patients by TME profile",
     xlab = "Patient", ylab = "Ward distance", sub = "", hang = -1)
dev.off()

# ── 3. GROUP COMPARISON: responder vs non-responder ──────────────────────────
message("3. Group comparison...")

feat_df <- as.data.frame(feat_mat) %>%
  tibble::rownames_to_column("patient_id") %>%
  left_join(CLINICAL[, c("patient_id", "pfs_group")], by = "patient_id")

# For each feature: mean in each group + Cohen's d (effect size)
group_comparison <- feat_df %>%
  pivot_longer(-c(patient_id, pfs_group), names_to = "feature", values_to = "value") %>%
  group_by(feature, pfs_group) %>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            sd_val   = sd(value, na.rm = TRUE),
            n        = n(), .groups = "drop") %>%
  pivot_wider(names_from = pfs_group,
              values_from = c(mean_val, sd_val, n)) %>%
  mutate(
    diff = mean_val_responder - mean_val_non_responder,
    pooled_sd = sqrt((sd_val_responder^2 + sd_val_non_responder^2) / 2),
    cohens_d = ifelse(pooled_sd > 0, diff / pooled_sd, 0)
  ) %>%
  arrange(desc(abs(cohens_d)))

write_tsv(group_comparison, file.path(TAB_DIR, "08_group_comparison.tsv"))

# Dot plot of top features by effect size
top_features <- group_comparison %>%
  filter(!is.na(cohens_d)) %>%
  head(20)

p_effect <- ggplot(top_features, aes(x = cohens_d,
                                      y = reorder(feature, abs(cohens_d)),
                                      color = cohens_d > 0)) +
  geom_point(size = 4) +
  geom_segment(aes(xend = 0, yend = reorder(feature, abs(cohens_d))), linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                     labels = c("Higher in non-responder", "Higher in responder"),
                     name = "Direction") +
  labs(title = "TME features: responder vs non-responder (Cohen's d)",
       x = "Cohen's d (effect size)", y = "") +
  theme_hgg()

save_plot(p_effect,
          file.path(strat_dir, "08_effect_size_comparison.pdf"),
          width = 10, height = 8)

# ── 4. SPEARMAN CORRELATIONS: features vs PFS/OS ────────────────────────────
message("4. Spearman correlations...")

feat_with_clinical <- as.data.frame(feat_mat) %>%
  tibble::rownames_to_column("patient_id") %>%
  left_join(CLINICAL[, c("patient_id", "pfs_months", "os_months")], by = "patient_id")

# Compute Spearman rho for each feature vs PFS and OS
feature_names <- colnames(feat_mat)

cor_results <- lapply(feature_names, function(f) {
  vals <- feat_with_clinical[[f]]
  if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0) {
    return(data.frame(feature = f, rho_PFS = NA, rho_OS = NA))
  }
  rho_pfs <- cor(vals, feat_with_clinical$pfs_months, method = "spearman", use = "complete.obs")
  rho_os  <- cor(vals, feat_with_clinical$os_months,  method = "spearman", use = "complete.obs")
  data.frame(feature = f, rho_PFS = rho_pfs, rho_OS = rho_os)
}) %>% bind_rows() %>%
  filter(!is.na(rho_PFS)) %>%
  arrange(desc(abs(rho_PFS)))

write_tsv(cor_results, file.path(TAB_DIR, "08_spearman_correlations.tsv"))

message("\nTop 15 features by |rho| with PFS:")
print(head(cor_results, 15))

# Correlation heatmap
cor_mat <- cor_results %>%
  head(25) %>%
  select(feature, rho_PFS, rho_OS) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix()

pheatmap(cor_mat,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_row = 9, fontsize_col = 11,
         cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = seq(-1, 1, length.out = 101),
         main = "Spearman rho: top TME features vs PFS/OS (n=8, descriptive)",
         filename = file.path(strat_dir, "08_spearman_heatmap.pdf"),
         width = 7, height = 10)

# ── 5. SCATTER PLOTS: top features vs PFS ────────────────────────────────────
message("5. Scatter plots...")

top_corr_features <- cor_results %>% head(6)

scatter_plots <- lapply(seq_len(nrow(top_corr_features)), function(i) {
  f   <- top_corr_features$feature[i]
  rho <- top_corr_features$rho_PFS[i]

  ggplot(feat_with_clinical, aes(x = .data[[f]], y = pfs_months,
                                  color = patient_id, label = patient_id)) +
    geom_point(size = 4) +
    geom_text_repel(size = 3, box.padding = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = PATIENT_COLORS) +
    labs(title = sprintf("%s (rho=%.2f)", f, rho),
         x = f, y = "PFS (months)") +
    theme_hgg() +
    theme(legend.position = "none", plot.title = element_text(size = 9))
})

save_plot(wrap_plots(scatter_plots, ncol = 3),
          file.path(strat_dir, "08_scatter_top_features_vs_pfs.pdf"),
          width = 16, height = 10)

# ── 6. CANDIDATE BIOMARKERS SUMMARY ──────────────────────────────────────────
message("6. Candidate biomarkers summary...")

# Combine effect size + correlation into a single ranking
candidates <- group_comparison %>%
  select(feature, cohens_d) %>%
  left_join(cor_results, by = "feature") %>%
  mutate(
    abs_d   = abs(cohens_d),
    abs_rho = abs(rho_PFS),
    # Combined rank: average of rank by |d| and rank by |rho|
    rank_d   = rank(-abs_d),
    rank_rho = rank(-abs_rho),
    combined_rank = (rank_d + rank_rho) / 2
  ) %>%
  arrange(combined_rank) %>%
  select(feature, cohens_d, rho_PFS, rho_OS, combined_rank)

write_tsv(candidates, file.path(TAB_DIR, "08_candidate_biomarkers.tsv"))

message("\n=== TOP 10 CANDIDATE BIOMARKERS ===")
message("(ranked by combined effect size + Spearman correlation with PFS)")
print(head(candidates, 10))

# ── 7. PATIENT PROFILES ─────────────────────────────────────────────────────
message("\n7. Patient profiles...")

profile_df <- CLINICAL %>%
  select(patient_id, pfs_months, os_months, pfs_group, mgmt_status, histology) %>%
  left_join(immune_ratio %>% select(patient_id, Immune, Tumor, immune_tumor_ratio), by = "patient_id") %>%
  left_join(score_table %>% select(patient_id, T_CD8_exhaustion, M2_M1_ratio,
                                    Tumor_immunogenicity, Tumor_immune_evasion, Immune_fraction),
            by = "patient_id")

write_tsv(profile_df, file.path(TAB_DIR, "08_patient_profiles.tsv"))

message("\n--- Patient profiles ---")
print(as.data.frame(profile_df))

message("\n=== Step 8 complete. ===")
message("Key outputs:")
message("  - Tables: 08_group_comparison.tsv, 08_spearman_correlations.tsv, 08_candidate_biomarkers.tsv")
message("  - Figures: 08_integrated_heatmap.pdf, 08_pca_patients_pfs.pdf, 08_effect_size_comparison.pdf")
message("  - Figures: 08_spearman_heatmap.pdf, 08_scatter_top_features_vs_pfs.pdf")
