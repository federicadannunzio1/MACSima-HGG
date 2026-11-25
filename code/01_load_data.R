# =============================================================================
# 01_load_data.R
# Load all CSV files, apply column harmonization, create Seurat object.
# Attaches clinical metadata (OS, PFS, MGMT, histology, PFS group) from start.
# Output: RDS with raw Seurat object (counts = raw MFI, no normalization yet).
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))
source(file.path(getwd(), "code/utils/plot_utils.R"))

library(Seurat)
library(readr)
library(dplyr)
library(stringr)

message("=== Step 1: Load data ===")

harm <- readRDS(file.path(RDS_DIR, "00_harmonization.rds"))
mapping_df     <- harm$mapping_df
shared_markers <- harm$shared_markers

# ── 1. Load each CSV, rename to canonical, add metadata ──────────────────────
load_roi <- function(sample_id, csv_path) {
  message("  Loading: ", sample_id)
  df <- as.data.frame(read_csv(csv_path, show_col_types = FALSE))

  # Mapping for this sample: raw_col → canonical (biological markers only)
  sm <- mapping_df[mapping_df$sample_id == sample_id & !mapping_df$is_control, ]
  old_to_new <- setNames(sm$canonical, sm$raw_col)

  for (old_name in intersect(sm$raw_col, colnames(df))) {
    colnames(df)[colnames(df) == old_name] <- old_to_new[[old_name]]
  }

  # Keep only shared biological markers
  bio_cols <- intersect(colnames(df), shared_markers)
  n_rows   <- nrow(df)

  pat_id <- SAMPLE_META$patient_id[SAMPLE_META$sample_id == sample_id]
  roi    <- SAMPLE_META$roi_id[SAMPLE_META$sample_id == sample_id]

  df_out <- df[, bio_cols, drop = FALSE]
  df_out$cell_id    <- paste0(sample_id, "_cell", seq_len(n_rows))
  df_out$sample_id  <- sample_id
  df_out$patient_id <- pat_id
  df_out$roi_id     <- roi

  # QC flag columns
  flags <- QC_FLAGS[[sample_id]]
  if (!is.null(flags)) {
    for (m in intersect(flags, shared_markers)) {
      df_out[[paste0("flag_", m)]] <- TRUE
    }
  }

  df_out
}

roi_list <- mapply(load_roi, names(CSV_FILES), CSV_FILES, SIMPLIFY = FALSE)

# ── 2. Merge all ROIs ─────────────────────────────────────────────────────────
all_cells <- do.call(dplyr::bind_rows, roi_list)

# Fill missing flag columns with FALSE
flag_cols <- grep("^flag_", colnames(all_cells), value = TRUE)
for (fc in flag_cols) {
  all_cells[[fc]] <- replace(all_cells[[fc]], is.na(all_cells[[fc]]), FALSE)
}

message(sprintf("Total cells loaded:  %d", nrow(all_cells)))
message(sprintf("Shared markers:      %d", length(shared_markers)))
message(sprintf("QC flag columns:     %d", length(flag_cols)))

# ── 3. Build MFI matrix (markers × cells) ────────────────────────────────────
mfi_matrix <- t(as.matrix(all_cells[, shared_markers]))
colnames(mfi_matrix) <- all_cells$cell_id

# ── 4. Create Seurat object ──────────────────────────────────────────────────
seurat_obj <- CreateSeuratObject(counts = mfi_matrix, assay = "MFI")

# Capture Seurat's internal feature names (after _ → - conversion)
seurat_markers <- rownames(seurat_obj)
message(sprintf("Seurat feature names: %d", length(seurat_markers)))

canonical_to_seurat <- setNames(seurat_markers, shared_markers)
saveRDS(canonical_to_seurat, file.path(RDS_DIR, "01_canonical_to_seurat_names.rds"))

# ── 5. Add sample + clinical metadata ────────────────────────────────────────
meta_df <- as.data.frame(all_cells[, c("sample_id", "patient_id", "roi_id", flag_cols)])
rownames(meta_df) <- all_cells$cell_id
seurat_obj <- AddMetaData(seurat_obj, metadata = meta_df)

# Attach clinical data per patient (OS, PFS, MGMT, histology, PFS group)
clinical_cols <- c("os_months", "pfs_months", "mgmt_status", "histology",
                   "cohort", "pfs_group", "os_group")
for (col in clinical_cols) {
  lookup <- setNames(CLINICAL[[col]], CLINICAL$patient_id)
  values <- lookup[seurat_obj$patient_id]
  names(values) <- colnames(seurat_obj)
  seurat_obj <- AddMetaData(seurat_obj, metadata = values, col.name = col)
}

# Verify metadata attached
stopifnot("patient_id not in metadata" = "patient_id" %in% colnames(seurat_obj@meta.data))
stopifnot("pfs_group not in metadata"  = "pfs_group"  %in% colnames(seurat_obj@meta.data))

# ── 6. Cell count summary ─────────────────────────────────────────────────────
cell_counts <- seurat_obj@meta.data %>%
  group_by(patient_id, sample_id, roi_id, pfs_group) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(patient_id, roi_id)

message("\n--- Cell counts per ROI ---")
print(as.data.frame(cell_counts))
message(sprintf("\nTotal cells: %d", sum(cell_counts$n_cells)))
message(sprintf("Patients: %d (responder=%d, non-responder=%d)",
                length(unique(cell_counts$patient_id)),
                sum(CLINICAL$pfs_group == "responder"),
                sum(CLINICAL$pfs_group == "non_responder")))

write_tsv(cell_counts, file.path(TAB_DIR, "01_cell_counts.tsv"))

# ── 7. Save ──────────────────────────────────────────────────────────────────
saveRDS(seurat_obj, file.path(RDS_DIR, "01_seurat_raw.rds"))
message("\nStep 1 complete. Saved: output/data/01_seurat_raw.rds")
