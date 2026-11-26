# Gating tree

Generated from `config/gating_tree.yaml` by `write_gating_mermaid()` (code/utils/gating_tree.R): do not edit by hand.
Edge text = gate applied to the cells of the parent (`+` positive, `-` negative, `rest` = cells not taken by an earlier sibling).

```mermaid
flowchart LR
  root["all cells"]
  immune["immune<br/><i>diagram: immune_cells</i>"]
  B_cell["B_cell<br/><i>diagram: B_cells</i>"]
  CD20neg["CD20neg<br/><i>diagram: CD20neg</i>"]
  NKT["NKT<br/><i>diagram: NKT_cells</i>"]
  T["T<br/><i>diagram: T_cells</i>"]
  T_CD8["T_CD8<br/><i>diagram: CD8_T_cells</i>"]
  T_CD8_exhausted["T_CD8_exhausted<br/><i>diagram: exhausted_CD8</i>"]
  T_CD8_rest["T_CD8_rest → T_CD8_other"]
  T_CD8_naive["T_CD8_naive<br/><i>diagram: naive_CD8</i>"]
  T_CD8_CM["T_CD8_CM<br/><i>diagram: CM_CD8</i>"]
  T_CD8_EM["T_CD8_EM<br/><i>diagram: EM_CD8</i>"]
  T_CD8_TEMRA["T_CD8_TEMRA<br/><i>diagram: TEMRA_CD8</i>"]
  T_CD4["T_CD4<br/><i>diagram: CD4_T_cells</i>"]
  T_reg["T_reg<br/><i>diagram: Tregs</i>"]
  T_CD4_nonreg["T_CD4_nonreg"]
  T_CD4_exhausted["T_CD4_exhausted<br/><i>diagram: exhausted_CD4</i>"]
  T_CD4_rest["T_CD4_rest → T_CD4_other"]
  T_CD4_naive["T_CD4_naive<br/><i>diagram: naive_CD4</i>"]
  T_CD4_CM["T_CD4_CM<br/><i>diagram: CM_CD4</i>"]
  T_CD4_EM["T_CD4_EM<br/><i>diagram: EM_CD4</i>"]
  T_CD4_TEMRA["T_CD4_TEMRA<br/><i>diagram: TEMRA_CD4</i>"]
  T_other["T_other"]
  T_DP["T_DP"]
  NK["NK<br/><i>diagram: NK_cells</i>"]
  NK_effector["NK_effector<br/><i>diagram: effector_NK</i>"]
  NK_immature["NK_immature<br/><i>diagram: immature_NK</i>"]
  myeloid["myeloid → Myeloid_other<br/><i>diagram: myeloids</i>"]
  Neutrophil["Neutrophil<br/><i>diagram: neutrophils</i>"]
  CD15neg["CD15neg"]
  Microglia["Microglia<br/><i>diagram: normal_microglia</i>"]
  CD14pos["CD14pos"]
  macrophage["macrophage → Macrophage_other<br/><i>diagram: macrophages</i>"]
  TAM_M2["TAM_M2<br/><i>diagram: M2_TAM</i>"]
  TAM_M1["TAM_M1<br/><i>diagram: M1</i>"]
  Monocyte["Monocyte<br/><i>diagram: monocytes</i>"]
  nonimmune["nonimmune<br/><i>diagram: CD45neg</i>"]
  Endothelial["Endothelial<br/><i>diagram: endothelial_cells</i>"]
  CD31neg_CD34neg["CD31neg_CD34neg → Unclassified_nonimmune<br/><i>diagram: CD31neg_CD34neg</i>"]
  PodVim["PodVim<br/><i>diagram: podoplanin_vimentin_pos</i>"]
  Tumor["Tumor<br/><i>diagram: tumor_cells</i>"]
  PodVim_GFAPneg["PodVim_GFAPneg"]
  Mesenchymal["Mesenchymal<br/><i>diagram: mesenchymal_cells</i>"]
  Fibroblast["Fibroblast<br/><i>diagram: fibroblasts</i>"]
  Tumor_microglia["Tumor_microglia<br/><i>diagram: tumor_microglia</i>"]
  PodVim_neg["PodVim_neg<br/><i>diagram: podoplanin_vimentin_neg</i>"]
  Astrocyte["Astrocyte<br/><i>diagram: astrocytes</i>"]
  root -->|"CD45+"| immune
  immune -->|"CD20+"| B_cell
  immune -->|"CD20-"| CD20neg
  CD20neg -->|"CD3+ CD56+"| NKT
  CD20neg -->|"CD3+ CD56-"| T
  T -->|"CD8a+ CD4-"| T_CD8
  T_CD8 -->|"PD1+ TIM3+ LAG3+"| T_CD8_exhausted
  T_CD8 -->|"rest"| T_CD8_rest
  T_CD8_rest -->|"CCR7+ CD45RA+"| T_CD8_naive
  T_CD8_rest -->|"CCR7+ CD45RO+"| T_CD8_CM
  T_CD8_rest -->|"CD45RO+ CCR7-"| T_CD8_EM
  T_CD8_rest -->|"CD45RA+ CCR7-"| T_CD8_TEMRA
  T -->|"CD4+ CD8a-"| T_CD4
  T_CD4 -->|"FoxP3+"| T_reg
  T_CD4 -->|"FoxP3-"| T_CD4_nonreg
  T_CD4_nonreg -->|"PD1+ TIM3+ LAG3+"| T_CD4_exhausted
  T_CD4_nonreg -->|"rest"| T_CD4_rest
  T_CD4_rest -->|"CCR7+ CD45RA+"| T_CD4_naive
  T_CD4_rest -->|"CCR7+ CD45RO+"| T_CD4_CM
  T_CD4_rest -->|"CD45RO+ CCR7-"| T_CD4_EM
  T_CD4_rest -->|"CD45RA+ CCR7-"| T_CD4_TEMRA
  T -->|"CD8a- CD4-"| T_other
  T -->|"CD4+ CD8a+"| T_DP
  CD20neg -->|"CD56+ CD3-"| NK
  NK -->|"CD16+"| NK_effector
  NK -->|"CD16-"| NK_immature
  CD20neg -->|"CD3- CD56-"| myeloid
  myeloid -->|"CD15+"| Neutrophil
  myeloid -->|"CD15-"| CD15neg
  CD15neg -->|"P2Y12+ CD14-"| Microglia
  CD15neg -->|"CD14+"| CD14pos
  CD14pos -->|"CD68+"| macrophage
  macrophage -->|"CD163+ CD206+"| TAM_M2
  macrophage -->|"CD163- CD206-"| TAM_M1
  CD14pos -->|"CD68-"| Monocyte
  root -->|"rest"| nonimmune
  nonimmune -->|"CD31+ CD34+"| Endothelial
  nonimmune -->|"rest"| CD31neg_CD34neg
  CD31neg_CD34neg -->|"(Podoplanin+ or Vimentin+)"| PodVim
  PodVim -->|"(GFAP+ or BCAN+)"| Tumor
  PodVim -->|"GFAP- BCAN-"| PodVim_GFAPneg
  PodVim_GFAPneg -->|"CD105+"| Mesenchymal
  PodVim_GFAPneg -->|"CD105-"| Fibroblast
  Fibroblast -->|"CD68+"| Tumor_microglia
  CD31neg_CD34neg -->|"rest"| PodVim_neg
  PodVim_neg -->|"GFAP+"| Astrocyte
```

## Populations of the diagram not implemented

| Diagram | Reason |
|---|---|
| oligodendrocytes | OLIG2 not in the panel |

## Deviations from the diagram

Same rows as `06_gating_deviations.tsv` (step 06). Implementation = gates from the root (`>` = child of).

| Item | Diagram | Implementation | Reason |
|---|---|---|---|
| T | CD45+ -> CD3+ -> T cells | CD45+ > CD20- > CD3+ CD56- | Under CD20- and CD56- so that B, T and NKT cells (CD20- CD3+ CD56+) do not overlap |
| T_CD8 | CD8+ -> CD8+ T cells | CD45+ > CD20- > CD3+ CD56- > CD8a+ CD4- | CD4- so that CD8 and CD4 T cells do not overlap; CD4+ CD8+ T cells: node T_DP |
| T_CD8_other | CD8+ T cells -> CCR7 -> CD45RA/CD45RO (parallel to the exhausted branch) | CD45+ > CD20- > CD3+ CD56- > CD8a+ CD4- > rest | One label per cell: exhausted has priority, differentiation states only for the other CD8 T cells; CD8 T cells in no state are T_CD8_other. CCR7+ CD45RA+ CD45RO+: CM wins over naive; CCR7- CD45RA+ CD45RO+: TEMRA wins over EM |
| T_CD4 | CD4+ -> CD4+ T cells | CD45+ > CD20- > CD3+ CD56- > CD4+ CD8a- | CD8a- so that CD4 and CD8 T cells do not overlap; CD4+ CD8+ T cells: node T_DP |
| T_CD4_nonreg | CD4+ T cells -> Tregs \| activated \| exhausted, and (parallel) CCR7 -> CD45RA/CD45RO | CD45+ > CD20- > CD3+ CD56- > CD4+ CD8a- > FoxP3- | One label per cell: Tregs first, then exhausted, then differentiation states for the other FoxP3- CD4 T cells (Tregs get no differentiation state) |
| T_CD4_other | not in the diagram | CD45+ > CD20- > CD3+ CD56- > CD4+ CD8a- > FoxP3- > rest | FoxP3- non-exhausted CD4 T cells in no differentiation state |
| T_other | not in the diagram | CD45+ > CD20- > CD3+ CD56- > CD8a- CD4- | CD4- CD8- T cells (owner, 25/11) |
| T_DP | not in the diagram | CD45+ > CD20- > CD3+ CD56- > CD4+ CD8a+ | CD4+ CD8+ T cells: no label (Unclassified), as in 432d1c3 |
| NK | CD20- -> CD56+ -> NK cells | CD45+ > CD20- > CD56+ CD3- | CD3- so that NK and NKT cells do not overlap. NK cells must be CD45+: without it they would be 49% of the immune cells, as in the lab's counts, but 56% of CD56+ CD45- cells are GFAP+ or BCAN+ (glial) |
| Myeloid_other | CD3-/CD56- -> myeloids -> neutrophils \| microglia \| CD14+ | CD45+ > CD20- > CD3- CD56- | Myeloid cells in no branch of the diagram (CD15- CD14- P2Y12-) are Myeloid_other. 432d1c3 split their CD68+ part into M1 (by CD80, not in the diagram), M2 and Macrophage_other, and called the rest Other_immune |
| Microglia | CD14-/CD15- -> Iba1+ -> normal microglia | CD45+ > CD20- > CD3- CD56- > CD15- > P2Y12+ CD14- | Iba1 not in the shared panel (9/15 ROI): P2Y12+ instead (owner, 25/11); P2Y12 is specific for CNS-resident microglia |
| Macrophage_other | CD14+ -> Iba1+/CD68+ -> macrophages -> M1 \| M2 | CD45+ > CD20- > CD3- CD56- > CD15- > CD14+ > CD68+ | Iba1 not in the shared panel: CD68+ only; CD15- because neutrophils are taken first. Macrophages that are neither M1 (CD163- CD206-) nor M2 (CD163+ CD206+) are Macrophage_other |
| Monocyte | CD14+ -> Iba1-/CD68- -> monocytes | CD45+ > CD20- > CD3- CD56- > CD15- > CD14+ > CD68- | Iba1 not in the shared panel: CD68- only; CD15- because neutrophils are taken first |
| nonimmune | CD45- -> non-immune cells | rest | All cells not taken by `immune`: CD45- and CD45-flagged cells (if CD45 is flagged in > 30% of all cells, immune = CD3+ or CD68+ or CD56+, see immune_fallback) |
| Unclassified_nonimmune | CD31-/CD34- | rest > rest | All the other CD45- cells, also CD31+ CD34- and CD31- CD34+ (the lab's counts partition every gateway: their leaves add up to the parent). CD45- cells in no leaf are Unclassified_nonimmune |
| PodVim | podoplanin+/vimentin+ | rest > rest > (Podoplanin+ or Vimentin+) | Podoplanin+ OR vimentin+ (owner, 26/11; decision 0006): with AND tumor cells were 1-15% of the cells per patient; the podoplanin-/vimentin- branch is then exactly the complement. Reactive astrocytes (GFAP+ vimentin+) cannot be told apart from tumor without OLIG2 |
| Tumor | podoplanin+/vimentin+ -> GFAP+/OLIG2+ -> tumor cells | rest > rest > (Podoplanin+ or Vimentin+) > (GFAP+ or BCAN+) | OLIG2 not in the panel: GFAP+ OR BCAN+ (owner, 25/11); BCAN (brevican) is brain-specific |
| PodVim_GFAPneg | podoplanin+/vimentin+ -> GFAP-/OLIG2- -> CD105 | rest > rest > (Podoplanin+ or Vimentin+) > GFAP- BCAN- | OLIG2 not in the panel: GFAP- and BCAN- (complement of the tumor gate) |
| Fibroblast | CD105- -> fibroblasts; CD105- -> Iba1+/CD68+ -> tumor microglia | rest > rest > (Podoplanin+ or Vimentin+) > GFAP- BCAN- > CD105- | CD105- cells are fibroblasts unless tumor microglia (child node) |
| Tumor_microglia | CD105- -> Iba1+/CD68+ -> tumor microglia | rest > rest > (Podoplanin+ or Vimentin+) > GFAP- BCAN- > CD105- > CD68+ | Iba1 not in the shared panel: CD68+ only (as for macrophages) |
| PodVim_neg | podoplanin-/vimentin- | rest > rest > rest | All the other CD45- CD31-/CD34- cells: podoplanin- and vimentin- (complement of the OR gate above), and cells with both markers QC-flagged |
| Astrocyte | podoplanin-/vimentin- -> GFAP+/OLIG2- -> astrocytes | rest > rest > rest > GFAP+ | OLIG2 not in the panel: GFAP+ only (OLIG2- not checked: may include GFAP+ OLIG2+ tumor cells). GFAP- cells here (oligodendrocytes or other) are Unclassified_nonimmune |
| thresholds | not specified (the lab gates by hand in its software) | automatic: two-component GMM per patient (P85 fallback) or negative population per ROI (`thresholds`) | Not yet validated against the lab's populations (F4a, F14) |
| QC flags | not specified | a QC-flagged marker counts as negative in the tree and is NA in the T annotations | Marker not measurable in that ROI (constants.R QC_FLAGS) |
| T annotations | activated (HLA-DR+, Ki-67+ or both) for CD4 non-Tregs and CD8, alternative to exhausted | columns T_activated (HLA-DR+) and T_proliferating (Ki-67+) for every T label and NKT, independent of the label | Readouts that overlay the differentiation states |

