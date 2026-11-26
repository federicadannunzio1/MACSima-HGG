# =============================================================================
# constants.R — Paths, metadata, marker synonyms, QC flags, clinical data
# HGG MACSima analysis — 8 patients, 15 ROIs
# =============================================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
# Project root (contains data/ and output/): HGG_BASE_DIR if set, otherwise the
# working directory (run the pipeline from the project root).
BASE_DIR  <- Sys.getenv("HGG_BASE_DIR", unset = getwd())
DATA_DIR  <- file.path(BASE_DIR, "data")
OUT_DIR   <- file.path(BASE_DIR, "output")
FIG_DIR   <- file.path(OUT_DIR, "figures")
RDS_DIR   <- file.path(OUT_DIR, "data")
TAB_DIR   <- file.path(OUT_DIR, "tables")

# ── Sample and clinical metadata (data/, not versioned) ──────────────────────
# Schema: docs/data.md; synthetic example: example/. Never put names or other
# identifiers in these files: patient_id only.
read_metadata <- function(file, numeric_cols = character()) {
  path <- file.path(DATA_DIR, file)
  if (!file.exists(path)) {
    stop("Missing ", path, ": see docs/data.md (format) and example/ (synthetic example)")
  }
  df <- read.delim(path, colClasses = "character", check.names = FALSE)
  for (col in numeric_cols) df[[col]] <- as.numeric(df[[col]])
  df
}

# One row per ROI: sample_id (N<year>_<case>_R<roi>), patient_id, roi_id, csv_file.
# csv_file = instrument export ("none" = raw MFI), relative to data/.
SAMPLE_META_FILE <- read_metadata("sample_metadata.tsv")
stopifnot(!anyDuplicated(SAMPLE_META_FILE$sample_id))
CSV_FILES <- setNames(as.list(file.path(DATA_DIR, SAMPLE_META_FILE$csv_file)),
                      SAMPLE_META_FILE$sample_id)
SAMPLE_META <- SAMPLE_META_FILE[, c("sample_id", "patient_id", "roi_id")]

# One row per patient: OS and PFS in months from surgery, MGMT, histology, cohort.
CLINICAL <- read_metadata("clinical.tsv", numeric_cols = c("os_months", "pfs_months"))
CLINICAL <- CLINICAL[, c("patient_id", "os_months", "pfs_months",
                         "mgmt_status", "histology", "cohort")]
stopifnot(!anyDuplicated(CLINICAL$patient_id),
          all(SAMPLE_META$patient_id %in% CLINICAL$patient_id))

# PFS group: >= cutoff = "responder", otherwise "non_responder" (cutoff and marker signatures: config/pipeline.yaml)
source(file.path(getwd(), "R/config.R"))
PIPELINE   <- read_pipeline_config()
SIGNATURES <- PIPELINE$signatures
PFS_CUTOFF <- PIPELINE$pfs_cutoff_months  # months
CLINICAL$pfs_group <- ifelse(CLINICAL$pfs_months >= PFS_CUTOFF, "responder", "non_responder")

# Median split for OS (secondary)
CLINICAL$os_group <- ifelse(CLINICAL$os_months >= median(CLINICAL$os_months),
                            "OS_high", "OS_low")

# ── Raw name corrections ─────────────────────────────────────────────────────
# Same protein, different raw names across MACSima batches.
# Applied AFTER stripping "Biomarker Exp" suffix and AFTER generic clone removal.
# Maps inconsistent cleaned names → unified canonical name.
# This is NOT a synonym dictionary — it only fixes batch-level naming inconsistencies.
RAW_NAME_CORRECTIONS <- c(
  # TIM-3: old batch = "TIM_3", new batch = "CD366_Isoform2"
  "CD366_Isoform2" = "TIM3",
  "TIM_3"          = "TIM3",

  # B7-H3: old batch = "B7H3_CD276", new batch = "CD276_B7_H3"
  "B7H3_CD276"     = "B7H3",
  "CD276_B7_H3"    = "B7H3",

  # CD45: some batches retain clone info after generic cleanup
  "CD45_PD7_26"    = "CD45",

  # CTLA-4: catalog suffix remains after cleanup
  "CTLA_4_CTLA4_6864R" = "CTLA_4",

  # ALK2/ACVR1: backslash in raw export
  "ALK2_ACVR1"     = "ALK2_ACVR1",

  # Sialic acid: typo in MACSima export ("siasilic")
  "siasilic_acid"  = "sialic_acid"
)

# ── QC flags: markers excluded per ROI (from collaborator's QC file) ──────────
# Source: info/Normalization - Excluded markers .xlsx
# These markers stay in the dataset but are flagged for interpretation.
# New samples (N21_838, N22_170, N23_818): no QC exclusion info available.
QC_FLAGS <- list(
  # --- Old samples (QC from Miltenyi normalization report) ---
  N17_135_R1 = c("FoxP3", "PAX5", "CD16", "CD133", "CD223"),
  N17_135_R2 = c("H3_3G34V", "CIITA", "PAX5", "CD20_Cytoplasmic", "AN2"),
  N21_898_R1 = c("H3_3G34R", "BCAN", "CD3", "Arginase_1", "CD34", "PAX5", "CD223", "BATF"),
  N22_1043_R1 = c("CD133", "CXCL10", "CD209", "PAX5", "AN2", "HLA_DR_DP_DQ"),
  N23_573_R1 = c("CCR7", "CD45", "HLA_ABC", "CD209", "PAX5", "CD274", "AN2", "Ki_67"),
  N23_573_R2 = c("TIGIT", "IgG", "PAX5", "CD274", "CD20_Cytoplasmic", "CD4", "AN2"),
  N23_573_R3 = c("CCR7", "B7H3", "H3_3G34V", "BCAN", "PAX5", "CD223", "CD15", "AN2"),
  N24_1017_R1 = c("CD105", "CIITA", "IgG", "CD56")

  # --- New samples: no QC exclusion info available ---
  # N21_838_R1, N21_838_R2, N21_838_R3: empty
  # N22_170_R1, N22_170_R2, N22_170_R3: empty
  # N23_818_R1: empty
)

# ── Markers to EXCLUDE from biological analysis (technical/control columns) ───
# DAPI channels = nuclear stain repeated each cycle (not a protein marker)
# IgG = isotype control (used for QC only, not biology)
CONTROL_MARKER_PATTERNS <- c("^DAPI", "^IgG$")

# ── Cell type gating rules (reference: info/gating strategy for the analysis.png)
# Full gating tree implemented in 06_gating.R. This summary is for documentation.
# Markers are canonical names. Gating is hierarchical (order matters).
CELLTYPE_RULES <- list(
  # Non-immune (CD45-)
  Endothelial  = list(gate = c("CD31", "CD34"),        measured = c()),
  Tumor        = list(gate = c("GFAP", "BCAN"),        measured = c("H3_3G34V", "H3_3G34R", "EZHIP", "p53", "GD2", "MGMT", "HLA_ABC", "HLA_DR", "NKG2D", "CIITA")),
  Mesenchymal  = list(gate = c("Podoplanin", "Vimentin", "CD105"), measured = c()),
  Fibroblast   = list(gate = c("Podoplanin", "Vimentin"), measured = c("CD105")),
  Stromal      = list(gate = c("Vimentin"),             measured = c()),
  # Immune - lymphoid (CD45+)
  B_cell       = list(gate = c("CD20_Cytoplasmic"),     measured = c()),
  NKT          = list(gate = c("CD3", "CD56"),          measured = c("NKG2D")),
  T_reg        = list(gate = c("CD3", "CD4", "FoxP3"), measured = c("CTLA_4")),
  T_CD4_exh    = list(gate = c("CD3", "CD4", "CD279", "TIM3", "CD223"), measured = c()),
  T_CD4_diff   = list(gate = c("CD3", "CD4", "CCR7", "CD45RA", "CD45RO"), measured = c("HLA_DR", "Ki_67")),
  T_CD8_exh    = list(gate = c("CD3", "CD8a", "CD279", "TIM3", "CD223"), measured = c()),
  T_CD8_diff   = list(gate = c("CD3", "CD8a", "CCR7", "CD45RA", "CD45RO"), measured = c("HLA_DR", "Ki_67")),
  NK_effector  = list(gate = c("CD56", "CD16"),         measured = c("NKG2D", "TIGIT")),
  NK_immature  = list(gate = c("CD56"),                 measured = c("NKG2D", "TIGIT")),
  # Immune - myeloid (CD45+ / CD3- / CD56-)
  Neutrophil   = list(gate = c("CD15"),                 measured = c()),
  Microglia    = list(gate = c("P2Y12"),                measured = c("TMEM119", "Arginase_1")),
  Monocyte     = list(gate = c("CD14"),                 measured = c()),
  TAM_M1       = list(gate = c("CD14", "CD68"),         measured = c("CD80")),
  TAM_M2       = list(gate = c("CD14", "CD68", "CD163", "CD206"), measured = c("Arginase_1")),
  Other_immune = list(gate = c("CD45"),                 measured = c())
)

# ── Analysis parameters ───────────────────────────────────────────────────────
# Parameters of steps 02-06 (PCA, UMAP, Louvain, QC, cofactor, plots): config/pipeline.yaml, key params
PARAMS <- PIPELINE$params

# ── Color palettes ────────────────────────────────────────────────────────────
PATIENT_COLORS <- c(
  N17_135  = "#E41A1C",
  N21_838  = "#377EB8",
  N21_898  = "#4DAF4A",
  N22_170  = "#FF7F00",
  N22_1043 = "#984EA3",
  N23_573  = "#A65628",
  N23_818  = "#F781BF",
  N24_1017 = "#66C2A5"
)

ROI_COLORS <- c(
  ROI1 = "#8DD3C7", ROI2 = "#FFFFB3", ROI3 = "#BEBADA"
)

PFS_GROUP_COLORS <- c(
  responder     = "#2166AC",
  non_responder = "#B2182B"
)

CELLTYPE_COLORS <- c(
  # Non-immune
  Tumor = "#E41A1C", Endothelial = "#A65628", Mesenchymal = "#FDBF6F",
  Fibroblast = "#CAB2D6", Stromal = "#D2B48C",
  Unclassified_nonimmune = "#CCCCCC",
  # CD4 T cells
  T_CD4_naive = "#DEEBF7", T_CD4_CM = "#9ECAE1", T_CD4_EM = "#4292C6",
  T_CD4_TEMRA = "#2171B5", T_CD4_exhausted = "#08306B", T_CD4_other = "#6BAED6",
  T_reg = "#9E9AC8",
  # CD8 T cells
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
  Unclassified = "#EEEEEE"
)
