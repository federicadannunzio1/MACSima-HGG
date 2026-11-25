# =============================================================================
# Tests of config/gating_tree.yaml. Run from the project root:
#   Rscript tests/test_gating_tree.R
# Needs only the yaml package (no data).
# =============================================================================

source("code/utils/gating_tree.R")

n_fail <- 0
check <- function(ok, what, detail = character()) {
  cat(if (ok) "PASS " else "FAIL ", what, "\n", sep = "")
  if (!ok) { n_fail <<- n_fail + 1; if (length(detail)) cat("     ", paste(detail, collapse = "\n      "), "\n") }
}

cfg   <- read_gating_tree("config/gating_tree.yaml")   # stops if malformed
nodes <- gating_nodes(cfg)
check(TRUE, "gating_tree.yaml is well formed")

# 1. Every population of the diagram is in the YAML (node, annotation or not_implemented),
#    and the YAML names no population that is not in the diagram.
png <- read.delim("tests/fixtures/png_nodes.tsv", colClasses = "character")
in_yaml <- unique(c(nodes$png[!is.na(nodes$png)],
                    unlist(lapply(cfg$annotations, `[[`, "png")),
                    vapply(cfg$not_implemented, `[[`, "", "png")))
check(length(setdiff(png$id, in_yaml)) == 0, "diagram populations all in the YAML",
      setdiff(png$id, in_yaml))
check(length(setdiff(in_yaml, png$id)) == 0, "YAML names only diagram populations",
      setdiff(in_yaml, png$id))

# 2. The MFI readouts of a node are those of the diagram.
mfi_set <- function(x) sort(strsplit(x, ",")[[1]])
withpng <- nodes[!is.na(nodes$png), ]
bad <- character()
for (i in seq_len(nrow(withpng))) {
  ref <- png$mfi[png$id == withpng$png[i]]
  if (!identical(mfi_set(withpng$mfi[i]), mfi_set(ref))) {
    bad <- c(bad, sprintf("%s: YAML [%s] diagram [%s]", withpng$id[i], withpng$mfi[i], ref))
  }
}
check(length(bad) == 0, "MFI readouts match the diagram", bad)

# 3. Every label in the tree has a colour, and every colour is a label of the tree.
labs <- c(cfg$default_label, unique(nodes$label[!is.na(nodes$label)]))
check(setequal(labs, names(cfg$labels)), "labels of the tree = labels with a colour",
      c(setdiff(labs, names(cfg$labels)), setdiff(names(cfg$labels), labs)))

# 4. Same labels as the hand-written tree of 06_gating.R at commit 432d1c3, on synthetic
#    cells (random positivity with 5% NA, and a case with CD45 NA in > 30% of cells).
fx  <- read.delim("tests/fixtures/gating_expected_432d1c3.tsv", check.names = FALSE)
ok  <- TRUE; detail <- character()
for (cs in unique(fx$case)) {
  d   <- fx[fx$case == cs, ]
  pos <- lapply(setNames(names(cfg$markers), names(cfg$markers)), function(m) as.logical(d[[m]]))
  got <- suppressMessages(apply_gating_tree(cfg, pos, nrow(d)))
  nd  <- sum(got != d$expected)
  if (nd > 0) { ok <- FALSE; detail <- c(detail, sprintf("%s: %d of %d cells differ", cs, nd, nrow(d))) }
}
check(ok, "labels identical to 06_gating.R of 432d1c3 (synthetic cells)", detail)

# 5. docs/gating_tree.md is up to date with the YAML.
tmp <- tempfile(fileext = ".md")
write_gating_mermaid(cfg, tmp)
check(identical(readLines(tmp), readLines("docs/gating_tree.md")), "docs/gating_tree.md up to date",
      "regenerate: Rscript -e 'source(\"code/utils/gating_tree.R\"); write_gating_mermaid(read_gating_tree(), \"docs/gating_tree.md\")'")

if (n_fail > 0) quit(status = 1)
cat("all tests passed\n")
