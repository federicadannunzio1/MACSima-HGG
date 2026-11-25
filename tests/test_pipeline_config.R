# Tests of config/pipeline.yaml and R/config.R (no data). Run from the project root:
#   Rscript tests/test_pipeline_config.R
n_fail <- 0
check <- function(ok, what) { cat(if (ok) "PASS " else "FAIL ", what, "\n", sep = ""); if (!ok) n_fail <<- n_fail + 1 }
source("R/config.R")
cfg <- read_pipeline_config("config/pipeline.yaml")
check(identical(cfg$pfs_cutoff_months, 6L) || identical(cfg$pfs_cutoff_months, 6), "PFS cut-off is 6 months")
check(identical(cfg$signatures$exhaustion, c("TIGIT", "CD279", "TIM3", "CTLA-4")), "exhaustion signature")
check(identical(cfg$signatures$m2, c("CD163", "CD206")) && identical(cfg$signatures$m1, "CD80"), "M1/M2 signatures")
check(identical(cfg$signatures$tumor_antigens, c("NKG2D", "GD2", "MGMT", "EZHIP")), "tumour antigens")
bad <- tempfile(fileext = ".yaml"); writeLines(c("pfs_cutoff_months: 6", "signatures:", "  m1: [CD80]"), bad)
check(inherits(try(read_pipeline_config(bad), silent = TRUE), "try-error"), "incomplete config stops")
# the code no longer hard-codes the signatures
src <- unlist(lapply(c("code/07_patient_features.R", "code/07b_functional_scores.R"), readLines))
check(!any(grepl('"CD163"|"TIGIT"|"NKG2D"|"HLA-ABC"', src)), "07/07b contain no hard-coded signature markers")
if (n_fail > 0) stop(n_fail, " test(s) failed")
cat("all tests passed\n")
