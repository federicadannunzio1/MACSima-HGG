# =============================================================================
# Tests of R/scores.R (no data). Run from the project root:
#   Rscript tests/test_scores.R
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tibble) })
source("R/scores.R")

n_fail <- 0
check <- function(ok, what) { cat(if (ok) "PASS " else "FAIL ", what, "\n", sep = ""); if (!ok) n_fail <<- n_fail + 1 }

check(is.na(safe_mean(numeric(0))), "safe_mean: empty vector gives NA")
check(is.na(safe_mean(c(NA_real_, NA_real_))), "safe_mean: all NA gives NA")
check(identical(safe_mean(c(1, NA, 3)), 2), "safe_mean: NA ignored")

d <- tibble(patient_id = c("A", "A", "B"), m1 = c(1, 3, 5), m2 = c(3, 5, NA), other = 9)
s <- score_by_patient(d, c("m1", "m2", "absent"))
check(identical(names(s), c("patient_id", "score", "n_cells")), "score_by_patient: columns")
check(isTRUE(all.equal(s$score, c(3, 5))), "score_by_patient: mean of per-cell marker means (absent marker ignored)")
check(identical(s$n_cells, c(2L, 1L)), "score_by_patient: cell counts per patient")
e <- score_by_patient(d, "absent")
check(nrow(e) == 0 && identical(names(e), c("patient_id", "score")), "score_by_patient: no available marker gives an empty tibble")

if (n_fail > 0) stop(n_fail, " test(s) failed")
cat("all tests passed\n")
