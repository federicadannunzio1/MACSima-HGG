# One-off migration (PR "Repository structure"): write data/sample_metadata.tsv and
# data/clinical.tsv from the constants.R of main BEFORE this change, where the metadata
# were written in the code. Keeps only the columns of docs/data.md (no names).
# Usage, from the project root, before pulling the new main:
#   git show main:code/utils/constants.R > /tmp/constants_old.R
#   Rscript tools/write_metadata_tsv.R /tmp/constants_old.R data/
args <- commandArgs(trailingOnly = TRUE)
env <- new.env()
src <- readLines(args[1])
src <- sub('^BASE_DIR  <- .*', 'BASE_DIR <- "__BASE__"', src)          # main: literal path
src <- src[!grepl('^ *unset = ', src)]                                  # fix-maintenance: 2-line form
src <- sub('^BASE_DIR  <- Sys.getenv\\("HGG_BASE_DIR",$', 'BASE_DIR <- "__BASE__"', src)
eval(parse(text = src), envir = env)
prefix <- file.path("__BASE__", "data", "")
meta <- env$SAMPLE_META
meta$csv_file <- sub(prefix, "", unlist(env$CSV_FILES[meta$sample_id]), fixed = TRUE)
stopifnot(identical(names(env$CSV_FILES), meta$sample_id))
clin <- env$CLINICAL[, c("patient_id", "os_months", "pfs_months", "mgmt_status", "histology", "cohort")]
write.table(meta[, c("sample_id", "patient_id", "roi_id", "csv_file")], file.path(args[2], "sample_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(clin, file.path(args[2], "clinical.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("written: ", nrow(meta), " ROIs, ", nrow(clin), " patients")
