# =============================================================================
# run_pipeline.R — Master script: runs all steps sequentially
# Run from the HGG project root:
#   Rscript code/run_pipeline.R
# =============================================================================

# Must be run from the project root (folder that contains code/)
if (!file.exists("code/utils/constants.R")) stop("Run from the project root: Rscript code/run_pipeline.R")
set.seed(42)

steps <- c(
  "code/00_harmonize_columns.R",
  "code/01_load_data.R",
  "code/02_qc.R",
  "code/03_normalize.R",
  "code/04_dimred_harmony.R",
  "code/05_clustering.R",
  # NOTE: 05b_annotation.R is manual — run separately after inspecting heatmaps
  "code/06_gating.R",
  "code/07_patient_features.R",
  "code/08_clinical_integration.R"
)

for (step in steps) {
  message("\n", strrep("=", 60))
  message("RUNNING: ", step)
  message(strrep("=", 60))
  source(step, local = new.env(parent = globalenv()))
}

source("code/utils/constants.R")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

message("\n", strrep("=", 60))
message("PIPELINE COMPLETE")
message(strrep("=", 60))
