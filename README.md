# MACSima-HGG

Analysis pipeline for MACSima multiplex imaging of paediatric high-grade gliomas:
per-cell marker intensities exported from the instrument are harmonised, quality
controlled, normalised and gated into cell types; per-patient composition and
functional scores are then related to clinical outcome.

The cohort is small (8 patients, 15 regions of interest). Every comparison with outcome
is descriptive: effect sizes and Spearman's rho, no p-values.

## Data

No data are in this repository. The pipeline expects, under the project root:

| Path | Content | Where it comes from |
|---|---|---|
| `data/<patient folder>/<ROI>/*.csv` | per-cell export of the instrument, `none` (raw MFI) | project shared drive (ask the maintainers) |
| `data/sample_metadata.tsv` | one row per ROI: `sample_id`, `patient_id`, `roi_id`, `csv_file` | written by the maintainers; schema in [docs/data.md](docs/data.md) |
| `data/clinical.tsv` | one row per patient: outcome and annotations | written by the maintainers; schema in [docs/data.md](docs/data.md) |

`data/` and `output/` are ignored by git. A synthetic example of the two metadata files
is in [`example/`](example/): it shows the format and contains no real patient data.

## Requirements

R ≥ 4.3 with Seurat 5, harmony, mclust, tidyverse, pheatmap, clustree, patchwork.
Exact versions will be pinned with `renv` (see [CONTRIBUTING](CONTRIBUTING.md)).

## How to run

From the project root (the folder that contains `code/`):

```sh
Rscript code/run_pipeline.R                         # data/ and output/ in the project root
HGG_BASE_DIR=/path/to/project Rscript code/run_pipeline.R   # data/ and output/ elsewhere
```

`run_pipeline.R` runs steps 00–05 and 06–08. `05b_annotation.R` (manual cluster labels)
is run separately after inspecting the step-05 heatmaps.

## Pipeline

```
CSV per ROI ─▶ 00 harmonise columns ─▶ 01 load (Seurat) ─▶ 02 cell QC ─▶ 03 arcsinh
            ─▶ 04 PCA / Harmony / UMAP ─▶ 05 Louvain clusters ─▶ (05b manual labels)
            ─▶ 06 marker gating ─▶ 07 patient features ─▶ 08 clinical integration
```

Inputs and outputs of each step: [docs/pipeline.md](docs/pipeline.md).
Analysis decisions and their status: [docs/decisions/](docs/decisions/).

## Outputs

- `output/data/*.rds`: intermediate objects, one per step.
- `output/tables/*.tsv`: cell counts, cofactors, gating thresholds and composition,
  patient features, group comparisons and correlations.
- `output/figures/`: QC, clustering, gating and clinical figures (PDF).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md): branches and pull requests, privacy rules, and how
to check that a change leaves the results unchanged.
