# =============================================================================
# Tests of R/io.R and R/markers.R (no data). Run from the project root:
#   Rscript tests/test_io_markers.R
# =============================================================================
n_fail <- 0
check <- function(ok, what) { cat(if (ok) "PASS " else "FAIL ", what, "\n", sep = ""); if (!ok) n_fail <<- n_fail + 1 }

RDS_DIR <- "some/dir"
source("R/io.R"); source("R/markers.R")

check(identical(rds_path("seurat_gated"), file.path("some/dir", "06_seurat_gated.rds")), "rds_path: known key")
check(inherits(try(rds_path("seurat_gatd"), silent = TRUE), "try-error"), "rds_path: unknown key stops")
check(!anyDuplicated(RDS_FILES) && all(grepl("\\.rds$", RDS_FILES)), "RDS_FILES: unique .rds names")
check(identical(to_seurat_name(c("HLA_DR", "CD3", "Ki_67")), c("HLA-DR", "CD3", "Ki-67")), "to_seurat_name")
check(identical(to_canonical_name(c("HLA-DR", "CD3")), c("HLA_DR", "CD3")), "to_canonical_name")
check(identical(to_canonical_name(to_seurat_name("CD20_Cytoplasmic")), "CD20_Cytoplasmic"), "round trip")

# every step reads/writes only registered keys
used <- unlist(regmatches(unlist(lapply(list.files("code", "\\.R$", full.names = TRUE), readLines)),
                          gregexpr('rds_path\\("[a-z_]+"\\)', unlist(lapply(list.files("code", "\\.R$", full.names = TRUE), readLines)))))
keys <- unique(gsub('rds_path\\("|"\\)', "", used))
check(length(keys) > 0 && all(keys %in% names(RDS_FILES)), "code/ uses only registered RDS keys")

if (n_fail > 0) stop(n_fail, " test(s) failed")
cat("all tests passed\n")
