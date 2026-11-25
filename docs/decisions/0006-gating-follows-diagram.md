# 0006 — Gating follows the reference diagram literally

Status: **Accepted provisionally** (2026-11-27, maintainer; to confirm with the owner).

## Context
The gating of commit 432d1c3 differed from `info/gating strategy for the analysis.png` in ways not written in
`06_gating_deviations.tsv`: tumor without the podoplanin+/vimentin+ gate, "or" instead of "and" for endothelial,
M2 and exhausted cells, branches not in the diagram. With those, endothelial cells were 34% of all cells
(lab's own counts: 1.3%) and 78% of CD8 T cells were "exhausted".

## Decision
Read the diagram literally, with one rule for its notation:
- "A+/B+" = both positive (as "CD3+/CD56+" for NKT cells); "A+, B+" = both (as "HLA-DR+, Ki-67+").
- The branches of a gateway partition the parent (the lab's counts add up at every gateway): a "−/−" branch
  takes every cell the "+/+" branch does not take.
- Hence: tumor = podoplanin+ and vimentin+, then GFAP+ or BCAN+ (OLIG2 not in the panel); endothelial =
  CD31+ and CD34+; M2 = CD163+ and CD206+ (mixed: `Macrophage_other`); exhausted = PD-1+, TIM-3+ and LAG-3+.
- Populations of the diagram implemented with the panel: astrocytes (GFAP+ in podoplanin−/vimentin−),
  tumor microglia (CD68+ in CD105−). Branches not in the diagram removed (`Stromal`, CD80-based M1,
  `Other_immune` → `Myeloid_other`); `T_other` kept (owner, 2026-11-25).
- NK cells stay CD45+ (without it they would match the lab's 49% of immune cells, but 56% of CD56+ CD45−
  cells are GFAP+ or BCAN+).
- Every difference from the diagram is a `deviation` field in `config/gating_tree.yaml`, written to
  `06_gating_deviations.tsv` and `docs/gating_tree.md`; a test fails if one is missing.

## Alternatives measured (seed 1, mean over patients)
| choice | tumor | endothelial | CD8 exhausted / CD8 |
|---|---:|---:|---:|
| 432d1c3 + fixes G1, G2 | 15% | 34% | 78% |
| this decision | 6% | 7% | 15% |
| podoplanin **or** vimentin | 18% | 7% | 15% |
| exhausted ≥ 2 of 3 receptors | 6% | 7% | 52% |
| lab (their counts) | 52% | 1.3% | — |

## Consequences
Composition changes a lot; no outcome signal is robust before or after (seeds, thresholds, flags, batch).
Half of the non-immune cells reach no leaf with automatic thresholds: the distance from the lab's counts comes
from the thresholds (decision 0002), not from the tree. Open: podoplanin/vimentin "and" vs "or", and the
number of receptors for exhaustion, to confirm with the owner.
