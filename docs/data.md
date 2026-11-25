# Input metadata: schema

Two tab-separated files in `data/` (ignored by git), read by `code/utils/constants.R`.
Synthetic examples with the same columns: `example/sample_metadata.tsv`,
`example/clinical.tsv`.

## `data/sample_metadata.tsv` — one row per ROI

| Column | Type | Meaning |
|---|---|---|
| `sample_id` | text, unique | `<patient_id>_R<n>`, e.g. `N17_135_R1` |
| `patient_id` | text | neuropathology code, e.g. `N17_135` |
| `roi_id` | text | `ROI1`, `ROI2`, ... |
| `csv_file` | text | path of the ROI export, relative to `data/` |

Row order is the order in which ROIs are loaded.

## `data/clinical.tsv` — one row per patient

| Column | Type | Meaning |
|---|---|---|
| `patient_id` | text, unique | as in `sample_metadata.tsv` |
| `os_months` | number | overall survival from surgery, months |
| `pfs_months` | number | progression-free survival from surgery, months |
| `mgmt_status` | text | `methyliert`, `nicht_methyliert`, `unknown` |
| `histology` | text | diagnosis label |
| `cohort` | text | study arm |

Derived in the code, not in the file: `pfs_group` (`responder` if `pfs_months` ≥
`PFS_CUTOFF`, else `non_responder`) and `os_group` (median split).

No names, dates of birth or other identifiers: `patient_id` only.
