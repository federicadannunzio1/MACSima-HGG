# 0002 — Positivity thresholds per ROI, from the negative population

Status: Proposed (2026-11-19)

## Context
Step 06 fits a two-component Gaussian mixture per patient and marker, with a P85 fallback.
Per-patient pooling mixes ROIs acquired in different conditions; the fit can place the
threshold in the tail of the distribution, and its result changes with package versions
(the same input gives different thresholds in another R environment). No threshold has
been checked on the images.

## Decision
- Thresholds are computed **per ROI** with a deterministic rule: mode of the negative
  population (kernel density) + 3 × MAD of the part below the mode. Default for markers
  expected to be rare (CD3, CD8a, CD4, FoxP3, CD56, CD68, P2Y12, CD163, CD80, TIGIT, PD-1,
  CD31).
- For markers that can be the majority of a ROI (CD45, GFAP, BCAN, Vimentin) the threshold
  is chosen and checked on the images.
- Kernel bandwidth and subsampling are fixed; ROIs with two peaks of similar height are
  flagged, because there the mode depends on the bandwidth.
- Thresholds are saved in a table (ROI, marker, threshold, method, % positive).

## Consequences
Rejected alternatives: mixture per ROI (unstable between bootstrap samples); isotype
(IgG) threshold (not comparable across fluorochromes). The rule fails when positives are
the majority: hence the separate handling of those four markers. Validation needs images
or the instrument's own per-cell population calls.
