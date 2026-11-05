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
  # --- N17/135 (Schmidt-Tudl) — 2 ROIs ---
  N17_135_R1 = file.path(DATA_DIR, "HGG N17 135 /ROI 1/N17 135 ROI 1 none.csv"),
  N17_135_R2 = file.path(DATA_DIR, "HGG N17 135 /ROI 2/N 17 135 ROI 2 none.csv"),

  # --- N21/838 (Gippert) — 3 ROIs ---
  N21_838_R1 = file.path(DATA_DIR, "HGG N21 838/ROI 1/Normalization HGG 21 838 ROI 1 none.csv"),
  N21_838_R2 = file.path(DATA_DIR, "HGG N21 838/ROI 2/HGG N21 838 ROI 2 Norm none.csv"),
  N21_838_R3 = file.path(DATA_DIR, "HGG N21 838/ROI 3/HGG N21 838 ROI 3 Norm none.csv"),

  # --- N21/898 (Möschl) — 1 ROI ---
  N21_898_R1 = file.path(DATA_DIR, "HGG N21 898 1-1-3/ROI 1/N21 898 ROI 1 none.csv"),

  # --- N22/170 (Eckert) — 3 ROIs ---
  N22_170_R1 = file.path(DATA_DIR, "HGG N22 170/ROI 1/HGG N22 170 ROI 1 none.csv"),
  N22_170_R2 = file.path(DATA_DIR, "HGG N22 170/ROI 2/HGG N22 170 ROI 2 none.csv"),
  N22_170_R3 = file.path(DATA_DIR, "HGG N22 170/ROI 3/HGG N22 170 ROI 3 none.csv"),

  # --- N22/1043 (Backes) — 1 ROI ---
  N22_1043_R1 = file.path(DATA_DIR, "HGG N22 1043 1-1-13/ROI1/N22 1043 none.csv"),

  # --- N23/573 (Zeun) — 3 ROIs ---
  N23_573_R1 = file.path(DATA_DIR, "HGG N23 573 1-1-6/ROI 1/N23 573 ROI 1 none.csv"),
  N23_573_R2 = file.path(DATA_DIR, "HGG N23 573 1-1-6/ROI2/N23 573 ROI2 none.csv"),
  N23_573_R3 = file.path(DATA_DIR, "HGG N23 573 1-1-6/ROI 3/N23 573 ROI 3 none.csv"),

  # --- N23/818 (Schütte) — 1 ROI ---
  # NB: filename has trailing space before .csv
  N23_818_R1 = file.path(DATA_DIR, "HGG N23 818/ROI 1/HGG N23 818 ROI 1 Norm none .csv"),

  # --- N24/1017 (Treitinger) — 1 ROI ---
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
  patient_name = c("Schmidt-Tudl", "Gippert", "Moeschl", "Eckert",
                    "Backes", "Zeun", "Schuette", "Treitinger"),
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

# ── Marker synonym mapping: raw column name → canonical name ──────────────────
# Applied AFTER stripping "Biomarker Exp" suffix.
# The generic cleaner handles most clone/catalog removal automatically;
# these entries handle edge cases the regex misses.
MARKER_SYNONYMS <- c(
  # TIM-3: two different raw names across sample batches
  "TIM_3 REAL818"            = "TIM3",
  "CD366_Isoform2 REAL818"   = "TIM3",

  # B7-H3 / CD276: two different raw names
  "B7H3 CD276"               = "B7H3",
  "CD276/B7-H3"              = "B7H3",

  # NKG2D: with or without clone suffix
  "NKG2D EPR24072-342"       = "NKG2D",
  "NKG2D"                    = "NKG2D",

  # CD45: clone name not auto-detected
  "CD45 2B11 PD7/26"         = "CD45",

  # CTLA-4: catalog suffix not auto-detected
  "CTLA-4 CTLA4/6864R"       = "CTLA_4",

  # ALK2/ACVR1: backslash in name from MACSima export
  "ALK2\\ ACVR1"             = "ALK2_ACVR1",

  # Sialic acid: typo in MACSima export
  "siasilic acid"            = "sialic_acid"
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

# ── Cell type annotation rules ────────────────────────────────────────────────
# Each entry: list(primary = required markers, secondary = supporting markers)
# Markers are canonical names. Annotation is hierarchical (order matters).
CELLTYPE_RULES <- list(
  Tumor        = list(primary = c("GFAP", "BCAN"),        secondary = c("H3_3K27M", "H3_3G34V", "H3_3G34R", "Vimentin")),
  T_CD8        = list(primary = c("CD3", "CD8a"),         secondary = c("CD45")),
  T_CD4        = list(primary = c("CD3", "CD4"),          secondary = c("CD45", "FoxP3")),
  T_reg        = list(primary = c("CD3", "CD4", "FoxP3"), secondary = c("CD45")),
  NK           = list(primary = c("CD56"),                secondary = c("CD16", "CD57")),
  Macrophage   = list(primary = c("CD68"),                secondary = c("CD163", "CD206", "Iba1")),
  Microglia    = list(primary = c("CD68", "P2Y12"),       secondary = c("TMEM119", "Iba1")),
  B_cell       = list(primary = c("CD20_Cytoplasmic"),    secondary = c("PAX5")),
  Endothelial  = list(primary = c("CD31"),                secondary = c("CD34")),
  Myeloid_other= list(primary = c("CD14"),                secondary = c("CD15", "CD16", "CD45"))
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
  Tumor         = "#E41A1C",
  T_CD8         = "#377EB8",
  T_CD4         = "#4DAF4A",
  T_reg         = "#984EA3",
  NK            = "#FF7F00",
  Macrophage    = "#A65628",
  Microglia     = "#F781BF",
  B_cell        = "#999999",
  Endothelial   = "#66C2A5",
  Myeloid_other = "#FC8D62",
  Unknown       = "#CCCCCC"
)
