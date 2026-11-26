# R environment

Reference versions (those of the repository owner's machine):

| package | version |
|---|---|
| R | 4.4.2 |
| Seurat | 5.3.1 |
| SeuratObject | 5.2.0 |
| harmony | 1.2.4 |
| mclust | 6.1.2 |
| clustree | 0.5.1 |
| pheatmap | 1.0.13 |
| ggplot2 | 4.0.3 |
| dplyr | 1.1.4 |
| tidyr | 1.3.2 |
| patchwork | 1.3.2 |
| readr | 2.2.0 |
| ggrepel | 0.9.8 |
| RColorBrewer | 1.1-3 |
| openxlsx | 4.2.8.1 |
| scales | 1.4.0 |
| stringr | 1.6.0 |

Also needed: `ggrastr` (step 05b), `yaml` (gating tree) and, for the tests, `testthat`.

## Pinning the versions with `renv` (owner, once)

`renv.lock` is the single source of the versions. It must be generated on the machine that runs
the pipeline, so that it records the versions that produced the results:

```r
# from the repository root, in R 4.4.2
install.packages("renv")
renv::init()        # creates renv/, .Rprofile and renv.lock from the packages the code uses
renv::snapshot()    # after installing anything new
```

Commit `renv.lock`, `renv/activate.R`, `renv/settings.json` and `.Rprofile` (not `renv/library/`). Anyone else
runs `renv::restore()` from the repository root.

## Known version constraints

- **harmony ≥ 2.0** renamed the `assay.use` argument of `RunHarmony()`, which step 04
  (`code/04_dimred_harmony.R`) uses: with harmony 2.x step 04 stops. Use harmony 1.2.x until the call is updated.
- **ggrastr** is not a dependency of Seurat; step 05b stops without it.

## Check done

Steps 06–08, same code and seed (`set.seed(1)`), run once in an environment with the versions above
and once with R 4.5.3 / Seurat 5.5.1 / harmony 2.0.5 (for 06–08, which do not call harmony): the tables and the
per-cell gating labels are identical (262 240 cells, 0 differences).
