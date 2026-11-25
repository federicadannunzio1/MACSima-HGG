# 0003 — QC-flagged markers become missing for their ROI

Status: Proposed (2026-11-19)

## Context
The collaborators' QC file lists, per ROI, markers excluded for technical reasons.
The pipeline stores them as `flag_<marker>` but the gating treats a flagged marker as
negative, so every population that depends on it is silently lost in that ROI.

## Decision
A flagged marker is missing (NA) in its ROI. Every population and score that depends on
it is NA in that ROI; patient-level fractions use only the ROIs where they can be
evaluated. The dependencies are declared in one table in the code.

## Consequences
Some ROIs and patients lose some populations, and the number of patients per feature is
reported. The exclusion lists of the newer batch must be added to `QC_FLAGS` (marker names
translated to the panel names, e.g. PD-1 as `CD279`).
