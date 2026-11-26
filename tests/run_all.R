# Runs every tests/test_*.R and the syntax check of all R files. Same command locally and in CI:
#   Rscript tests/run_all.R        (from the repository root; no data needed)
# Exit status 1 if anything fails.
failed <- character()

for (f in list.files(c("code", "R", "tests"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)) {
  ok <- tryCatch({ parse(f); TRUE }, error = function(e) { message("SYNTAX ", f, ": ", conditionMessage(e)); FALSE })
  if (!ok) failed <- c(failed, paste("syntax", f))
}
cat("syntax check done\n")

for (t in list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE)) {
  cat("\n== ", t, "\n", sep = "")
  rc <- system2(file.path(R.home("bin"), "Rscript"), t)
  if (rc != 0) failed <- c(failed, t)
}

if (length(failed)) { cat("\nFAILED:\n", paste0("  ", failed, collapse = "\n"), "\n", sep = ""); quit(status = 1) }
cat("\nall checks passed\n")
