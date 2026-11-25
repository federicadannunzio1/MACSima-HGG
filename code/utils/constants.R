# =============================================================================
# constants.R — Paths, metadata, marker synonyms, QC flags, clinical data
# HGG MACSima analysis — 8 patients, 15 ROIs
# =============================================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR  <- "/Users/federicadannunzio/Documents/projects/collaborations/IgnazioCaruana/HGG"
DATA_DIR  <- file.path(BASE_DIR, "data")
OUT_DIR   <- file.path(BASE_DIR, "output")
FIG_DIR   <- file.path(OUT_DIR, "figures")
RDS_DIR   <- file.path(OUT_DIR, "data")
TAB_DIR   <- file.path(OUT_DIR, "tables")

# ── CSV file paths (none = raw MFI, no normalization) ─────────────────────────
# sample_id naming: N<year>_<case>_R<roi>
# NOTE: file names from MACSima are inconsistent (spaces, "Norm", etc.)
#       Each path is verified against the actual file on disk.
CSV_FILES <- list(
  # --- N17_135 — 2 ROIs ---
  N17_135_R1 = file.path(DATA_DIR, "HGG N17 135 /ROI 1/N17 135 ROI 1 none.csv"),
  N17_135_R2 = file.path(DATA_DIR, "HGG N17 135 /ROI 2/N 17 135 ROI 2 none.csv"),

  # --- N21_838 — 3 ROIs ---
  N21_838_R1 = file.path(DATA_DIR, "HGG N21 838/ROI 1/Normalization HGG 21 838 ROI 1 none.csv"),
  N21_838_R2 = file.path(DATA_DIR, "HGG N21 838/ROI 2/HGG N21 838 ROI 2 Norm none.csv"),
  N21_838_R3 = file.path(DATA_DIR, "HGG N21 838/ROI 3/HGG N21 838 ROI 3 Norm none.csv"),

  # --- N21_898 — 1 ROI ---
  N21_898_R1 = file.path(DATA_DIR, "HGG N21 898 1-1-3/ROI 1/N21 898 ROI 1 none.csv"),

  # --- N22_170 — 3 ROIs ---
  N22_170_R1 = file.path(DATA_DIR, "HGG N22 170/ROI 1/HGG N22 170 ROI 1 none.csv"),
  N22_170_R2 = file.path(DATA_DIR, "HGG N22 170/ROI 2/HGG N22 170 ROI 2 none.csv"),
  N22_170_R3 = file.path(DATA_DIR, "HGG N22 170/ROI 3/HGG N22 170 ROI 3 none.csv"),

  # --- N22_1043 — 1 ROI ---
  N22_1043_R1 = file.path(DATA_DIR, "HGG N22 1043 1-1-13/ROI1/N22 1043 none.csv"),

  # --- N23_573 — 3 ROIs ---
  N23_573_R1 = file.path(DATA_DIR, "HGG N23 573 1-1-6/ROI 1/N23 573 ROI 1 none.csv"),
  N23_573_R2 = file.path(DATA_DIR, "HGG N23 573 1-1-6/ROI2/N23 573 ROI2 none.csv"),
  N23_573_R3 = file.path(DATA_DIR, "HGG N23 573 1-1-6/ROI 3/N23 573 ROI 3 none.csv"),

  # --- N23_818 — 1 ROI ---
  # NB: filename has trailing space before .csv
  N23_818_R1 = file.path(DATA_DIR, "HGG N23 818/ROI 1/HGG N23 818 ROI 1 Norm none .csv"),

  # --- N24_1017 — 1 ROI ---
  N24_1017_R1 = file.path(DATA_DIR, "HGG N24 1017/ROI 1/N24 1017 ROI 1 none.csv")
)

# ── Sample metadata ──────────────────────────────────────────────────────────
# patient_id = neuropath code (unique per patient)
SAMPLE_META <- data.frame(
  sample_id  = names(CSV_FILES),
  patient_id = c(
    "N17_135", "N17_135",                          # 2 ROIs
    "N21_838", "N21_838", "N21_838",               # 3 ROIs
    "N21_898",                                      # 1 ROI
    "N22_170", "N22_170", "N22_170",               # 3 ROIs
    "N22_1043",                                     # 1 ROI
    "N23_573", "N23_573", "N23_573",               # 3 ROIs
    "N23_818",                                      # 1 ROI
    "N24_1017"                                      # 1 ROI
  ),
  roi_id = c(
    "ROI1", "ROI2",
    "ROI1", "ROI2", "ROI3",
    "ROI1",
    "ROI1", "ROI2", "ROI3",
    "ROI1",
    "ROI1", "ROI2", "ROI3",
    "ROI1",
    "ROI1"
  ),
  stringsAsFactors = FALSE
)

# ── Clinical data per patient ────────────────────────────────────────────────
# Source: info/Patienten_Neuropatho_Immunovac.xlsx
# OS and PFS in months from surgery (OP date)
# PFS_group: ≥6 months = "responder", <6 months = "non-responder"
#   Natural gap in PFS data between 5.18 and 6.59 months → clean split 4 vs 4
CLINICAL <- data.frame(
  patient_id = c("N17_135", "N21_838", "N21_898", "N22_170",
                  "N22_1043", "N23_573", "N23_818", "N24_1017"),
  os_months  = c(32.89, 10.98, 8.36, 8.98, 12.33, 9.67, 7.08, 7.48),
  pfs_months = c(16.82,  4.03, 8.30, 3.44,  6.59, 7.08, 4.26, 5.18),
  mgmt_status = c("unknown", "methyliert", "methyliert", "nicht_methyliert",
                   "methyliert", "nicht_methyliert", "nicht_methyliert", "methyliert"),
  histology = c(
    "GBM_IV",
    "diffuse_hemispheric_HGG_H3_G34",
    "infant_HGG",
    "DMG_H3_K27",
    "diffuse_hemispheric_HGG_H3_G34",
    "diffuse_pediatric_HGG_H3_IDH_wt",
    "diffuse_pediatric_HGG_H3_IDH_wt",
    "diffuse_hemispheric_HGG_H3_G34"
  ),
  cohort = c("Vorphase", "Studienpatienten", "Studienpatienten", "Studienpatienten",
             "Studienpatienten", "Studienpatienten", "Studienpatienten", "Studienpatienten"),
  stringsAsFactors = FALSE
)

PFS_CUTOFF <- 6  # months
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
PARAMS <- list(
  # Arcsinh normalization
  cofactor_quantile  = 0.50,   # P50 per patient per marker

  # Cell QC
  min_total_mfi_pct  = 0.01,   # remove cells below 1st percentile of total MFI

  # PCA
  n_pcs              = 25,

  # Harmony
  harmony_group      = "patient_id",  # batch = patient

  # UMAP
  umap_neighbors     = 30,
  umap_mindist       = 0.3,

  # Clustering
  louvain_resolutions = c(0.05, 0.1, 0.15, 0.2, 0.3),  # added 0.3 for more granularity with 8 patients

  # Plots
  point_size         = 0.2,    # smaller points for larger dataset
  dpi                = 150
)

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
