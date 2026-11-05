# =============================================================================
# 00_harmonize_columns.R
# Read all CSVs, resolve column name synonyms, identify shared marker set,
# report column discrepancies. Does NOT modify data — outputs a harmonization
# map saved as RDS and a summary table as TSV.
# =============================================================================

source(file.path(getwd(), "code/utils/constants.R"))

library(readr)
library(stringr)
library(dplyr)

message("=== Step 0: Column harmonization ===")

# ── 1. Read headers only ──────────────────────────────────────────────────────
raw_headers <- lapply(CSV_FILES, function(path) {
  line <- readLines(path, n = 1)
  cols <- str_split(line, ",")[[1]]
  str_remove_all(cols, '"')
})

# ── 2. Strip "Biomarker Exp" suffix → intermediate name ───────────────────────
strip_suffix <- function(x) str_trim(str_remove(x, "Biomarker Exp$"))

stripped <- lapply(raw_headers, function(cols) {
  setNames(strip_suffix(cols), cols)
})

# ── 3. Apply synonym mapping → canonical name ─────────────────────────────────
clean_to_canonical <- function(stripped_name) {
  # Apply explicit synonym map first
  if (stripped_name %in% names(MARKER_SYNONYMS)) {
    return(MARKER_SYNONYMS[[stripped_name]])
  }
  # Remove clone catalog numbers: words starting with REA, REAL, EPR, HD, etc.
  tokens   <- str_split(stripped_name, "\\s+")[[1]]
  is_clone <- str_detect(tokens, "^(REA|REAL|EPR|HD|BLR|RM|VIMC|POLYCLONAL|[0-9]|L50)")
  is_clone <- is_clone | str_detect(tokens, "^[0-9]")
  clean    <- str_c(tokens[!is_clone], collapse = "_")
  # Replace problematic characters
  clean    <- str_replace_all(clean, "[./\\\\-]", "_")
  clean    <- str_replace_all(clean, "__+", "_")
  str_remove(clean, "_$")
}

# Build full mapping table
build_mapping <- function(sample_id, raw_cols, stripped_cols) {
  canonical <- sapply(stripped_cols, clean_to_canonical, USE.NAMES = FALSE)
  data.frame(
    sample_id   = sample_id,
    raw_col     = raw_cols,
    stripped    = stripped_cols,
    canonical   = canonical,
    stringsAsFactors = FALSE
  )
}

mapping_list <- mapply(function(sid, raw_h, strip_h) {
  build_mapping(sid, raw_h, strip_h)
}, names(raw_headers), raw_headers, stripped, SIMPLIFY = FALSE)

mapping_df <- bind_rows(mapping_list)

# ── 4. Identify control vs biological markers ─────────────────────────────────
is_control <- function(canonical) {
  any(sapply(CONTROL_MARKER_PATTERNS, function(pat) str_detect(canonical, pat)))
}
mapping_df <- mapping_df %>%
  mutate(is_control = sapply(canonical, is_control))

bio_mapping  <- mapping_df %>% filter(!is_control)
ctrl_mapping <- mapping_df %>% filter(is_control)

# ── 5. Find shared markers (present in ALL ROIs) ─────────────────────────────
canonical_per_sample <- bio_mapping %>%
  group_by(sample_id) %>%
  summarise(markers = list(unique(canonical)), .groups = "drop")

shared_markers <- Reduce(intersect, canonical_per_sample$markers)
all_markers    <- unique(bio_mapping$canonical)
subset_markers <- setdiff(all_markers, shared_markers)

message(sprintf("Total unique biological markers (canonical): %d", length(all_markers)))
message(sprintf("Markers present in ALL %d ROIs (shared panel): %d",
                length(CSV_FILES), length(shared_markers)))
message(sprintf("Markers present only in subset of ROIs:       %d", length(subset_markers)))

# ── 6. Report which ROIs have each subset marker ─────────────────────────────
if (length(subset_markers) > 0) {
  message("\n--- Subset-specific markers ---")
  for (m in sort(subset_markers)) {
    present_in <- bio_mapping %>%
      filter(canonical == m) %>%
      pull(sample_id) %>% unique()
    n_present <- length(present_in)
    message(sprintf("  %-30s -> %d/%d ROIs: %s",
                    m, n_present, length(CSV_FILES),
                    paste(present_in, collapse = ", ")))
  }
}

# ── 7. Report QC-flagged markers ─────────────────────────────────────────────
message("\n--- QC-flagged markers per ROI ---")
for (sid in names(QC_FLAGS)) {
  flags <- QC_FLAGS[[sid]]
  message(sprintf("  %s: %s", sid, paste(flags, collapse = ", ")))
}
samples_no_flags <- setdiff(names(CSV_FILES), names(QC_FLAGS))
if (length(samples_no_flags) > 0) {
  message(sprintf("  [No QC flags]: %s", paste(samples_no_flags, collapse = ", ")))
}

# ── 8. Save outputs ──────────────────────────────────────────────────────────
saveRDS(list(
  mapping_df     = mapping_df,
  shared_markers = shared_markers,
  subset_markers = subset_markers,
  all_markers    = all_markers,
  ctrl_markers   = unique(ctrl_mapping$canonical)
), file = file.path(RDS_DIR, "00_harmonization.rds"))

write_tsv(mapping_df, file.path(TAB_DIR, "00_column_mapping.tsv"))
write_tsv(
  data.frame(canonical = shared_markers),
  file.path(TAB_DIR, "00_shared_markers.tsv")
)

message("\nStep 0 complete. Outputs saved to output/data/ and output/tables/")
