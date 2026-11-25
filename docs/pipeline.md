# Pipeline: inputs and outputs of each step

All paths are relative to `BASE_DIR` (project root, or `HGG_BASE_DIR`). RDS files are in
`output/data/`, tables in `output/tables/`, figures in `output/figures/`. Each step reads
`code/utils/constants.R` (paths, metadata, marker names, QC flags, parameters).

| Step | Script | Reads | Writes | What it does |
|---|---|---|---|---|
| 00 | `00_harmonize_columns.R` | CSV headers (`CSV_FILES`) | `00_harmonization.rds`; `00_column_mapping.tsv`, `00_shared_markers.tsv` | Strips the "Biomarker Exp" suffix and clone names, applies name corrections, keeps the markers shared by all ROIs |
| 01 | `01_load_data.R` | CSVs, `00_harmonization.rds`, `SAMPLE_META`, `CLINICAL` | `01_seurat_raw.rds`, `01_canonical_to_seurat_names.rds`; `01_cell_counts.tsv` | Builds the Seurat object (counts = raw MFI), adds sample and clinical metadata and the `flag_<marker>` columns from `QC_FLAGS` |
| 02 | `02_qc.R` | `01_seurat_raw.rds` | `02_seurat_qc.rds`; `02_cell_counts_postQC.tsv`; `figures/qc/02_qc_panel.pdf` | Removes cells below the 1st percentile of total MFI |
| 03 | `03_normalize.R` | `02_seurat_qc.rds` | `03_seurat_normalized.rds`; `03_arcsinh_cofactors.tsv`; `figures/qc/03_*.pdf` | arcsinh(MFI / cofactor), cofactor = median over patients of the per-patient P50 |
| 04 | `04_dimred_harmony.R` | `03_seurat_normalized.rds` | `04_seurat_dimred.rds`; `figures/dimred/` | PCA, Harmony by patient, UMAP |
| 05 | `05_clustering.R` | `04_seurat_dimred.rds` | `05_seurat_clustered.rds`; `05_composition_per_patient.tsv`; `figures/clustering/05_*.pdf` | Louvain at several resolutions, clustree, automatic labels |
| 05b | `05b_annotation.R` (manual, run separately) | `05_seurat_clustered.rds` | `05b_seurat_annotated.rds`; `05b_composition_annotated.tsv`; `figures/clustering/05b_*.pdf` | Manual cluster labels (cluster IDs written in the script) |
| 06 | `06_gating.R` | `05_seurat_clustered.rds` | `06_seurat_gated.rds`, `06_mfi_gated.rds`; `06_gating_summary.tsv`, `06_composition_gating.tsv`, `06_gating_deviations.tsv`; `figures/qc/06_*.pdf`, `figures/clustering/06_*.pdf` | Per-patient two-component Gaussian mixture threshold per marker (P90 if the fit fails, P85 if more than 60% of cells are positive); gating tree of `info/gating strategy for the analysis.png`, one cell type per cell; T-cell activation (HLA-DR) and proliferation (Ki-67) as extra columns. `Mclust` starts from a random subset of cells: results change between runs unless a seed is set |
| 07 | `07_patient_features.R` | `06_seurat_gated.rds` | `07_patient_features.rds`; `07_functional_scores.tsv`, `07_immune_tumor_ratio.tsv` | Patient-level cell-type fractions and functional scores |
| 07b | `07b_functional_scores.R` (not in `run_pipeline.R`) | `06_seurat_gated.rds` | `07b_functional_scores.tsv` | Earlier version of the functional scores |
| 07s | `07_stratification.R` (not in `run_pipeline.R`) | `06_seurat_gated.rds` | `07_immune_tumor_ratio.tsv` (same name as step 07); `figures/stratification/07_*.pdf` | Outcome-blind patient clustering and PCA |
| 08 | `08_clinical_integration.R` | `07_patient_features.rds`, `CLINICAL` | `08_group_comparison.tsv`, `08_spearman_correlations.tsv`, `08_candidate_biomarkers.tsv`, `08_patient_profiles.tsv`; `figures/stratification/08_*.pdf` | Responder vs non-responder effect sizes, Spearman's rho with PFS and OS, ranking of features |

Notes

- Step 06 reads the step-05 object, not the 05b one: manual labels do not enter gating.
- `07_stratification.R` and step 07 write a table with the same name; running both
  overwrites it.
- `run_pipeline.R` runs 00–05 and 06–08.
