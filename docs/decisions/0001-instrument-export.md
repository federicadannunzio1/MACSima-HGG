# 0001 — Instrument export: fixed list of per-cell fields

Status: Proposed (2026-11-19, list narrowed 2026-11-20)

## Context
The per-ROI CSV exports differ between acquisition batches: different columns, some
files with only a subset of markers ("only" files whose cells cannot be matched to the
main file), no coordinates and no cell size. The shared-marker intersection (step 00)
therefore drops markers that exist in the panel, and cell-level QC can only use total MFI.

## Decision
Ask for one uniform export of every ROI with: nucleus centroid (Nuc X, Nuc Y),
Nucleus Size, Cell Size, and the complete set of "Biomarker Exp" columns of the panel
(including OLIG2). Quality is used only as a margin check. ROI index and Convexity are
not requested.

## Consequences
- Step 00 can be replaced by a fixed panel with a column check.
- Coordinates allow checking thresholds on images and computing densities.
- Until the new export arrives, the current files are used as they are.
