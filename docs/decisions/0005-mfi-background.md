# 0005 — MFI scores subtract the ROI background (mode)

Status: Proposed (2026-11-19)

## Context
Functional scores average marker intensities over cells. Those averages follow the
background level of each ROI (low percentiles, isotype control) as much as the biology,
and some files were exported with a different instrument normalisation.

## Decision
Before averaging, subtract for each ROI and marker the mode of the negative population
(same estimate as decision 0002). In the longer term, replace intensity scores with
fractions of positive cells once the thresholds are validated.

## Consequences
The subtraction reduces but does not remove the dependence on background; with eight
patients no variant can be shown to be unbiased. ROIs from differently normalised exports
are flagged and not compared on MFI.
