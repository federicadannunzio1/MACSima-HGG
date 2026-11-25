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

# 4. Expected labels of hand-written cell profiles (tests/fixtures/gating_cases.tsv: the
#    markers listed are positive, all the others negative), read from the diagram.
cases <- read.delim("tests/fixtures/gating_cases.tsv", colClasses = "character")
pos <- lapply(setNames(names(cfg$markers), names(cfg$markers)), function(m)
  vapply(strsplit(cases$positive, ","), function(p) m %in% p, TRUE))
got <- suppressMessages(apply_gating_tree(cfg, pos, nrow(cases)))
bad <- which(got != cases$expected)
check(length(bad) == 0, sprintf("expected labels of %d cell profiles", nrow(cases)),
      sprintf("%s: got %s, expected %s", cases$case[bad], got[bad], cases$expected[bad]))

# 4b. Every difference from the diagram has a `deviation` (so 06_gating_deviations.tsv is complete).
#     Node with `png`: its gate, together with the gates of the nodes without `png` between it and
#     its nearest ancestor with `png`, must equal the diagram's gate (fixture column `gate`) and
#     that ancestor must be the diagram's parent; otherwise the node or one of those nodes needs a
#     `deviation`. Node without `png` that gives a label or ends a branch, and node with `png`
#     that gives a label to cells of no child: own `deviation`.
tokens <- function(id) {
  i <- match(id, nodes$id); g <- nodes$gate[i]
  if (g == "rest") return("rest")
  gsub("[()]", "", unlist(strsplit(gsub(" or ", "|", g), " ")))
}
png_anc <- function(id) {                    # nearest ancestor with png, and the nodes in between
  between <- character(); p <- nodes$parent[match(id, nodes$id)]
  while (p != "root" && is.na(nodes$png[match(p, nodes$id)])) { between <- c(between, p); p <- nodes$parent[match(p, nodes$id)] }
  list(anc = if (p == "root") "root" else nodes$png[match(p, nodes$id)], between = between)
}
has_dev <- function(ids) any(!is.na(nodes$dev_reason[match(ids, nodes$id)]))
missing <- character()
for (k in seq_len(nrow(nodes))) {
  id <- nodes$id[k]
  if (!is.na(nodes$png[k])) {
    a   <- png_anc(id)
    yt  <- sort(unique(unlist(lapply(c(id, a$between), tokens))))
    ref <- png[png$id == nodes$png[k], ]
    dt  <- sort(strsplit(ref$gate, " ")[[1]])
    if ((!identical(yt, dt) || a$anc != ref$parent) && !has_dev(c(id, a$between)))
      missing <- c(missing, sprintf("%s: [%s] vs diagram [%s] under %s", id, paste(yt, collapse = " "), paste(dt, collapse = " "), ref$parent))
    if (!is.na(nodes$label[k]) && nodes$has_children[k] && !has_dev(id))
      missing <- c(missing, sprintf("%s: label %s for cells of no child", id, nodes$label[k]))
  } else if ((!is.na(nodes$label[k]) || !nodes$has_children[k]) && !has_dev(id)) {
    missing <- c(missing, sprintf("%s: not in the diagram", id))
  }
}
check(length(missing) == 0, "every difference from the diagram has a deviation", missing)
dv <- gating_deviations(cfg)
check(nrow(dv) == sum(!is.na(nodes$dev_reason)) + length(cfg$general_deviations) + length(cfg$not_implemented) &&
        !anyNA(dv$reason) && all(nzchar(dv$reason)), "06_gating_deviations.tsv: one row per deviation, all with a reason")

# 4c. Label groups used by steps 07/07b.
imm <- gating_labels_under(cfg, "immune")
check(all(c("B_cell", "NKT", "T_CD8_exhausted", "Myeloid_other") %in% imm) && !any(c("Tumor", "Endothelial") %in% imm),
      "gating_labels_under: immune labels")

# 4d. Threshold method.
check(cfg$thresholds$method %in% c("gmm2", "negpop"), "thresholds$method is gmm2 or negpop")

# 5. docs/gating_tree.md is up to date with the YAML.
tmp <- tempfile(fileext = ".md")
write_gating_mermaid(cfg, tmp)
check(identical(readLines(tmp), readLines("docs/gating_tree.md")), "docs/gating_tree.md up to date",
      "regenerate: Rscript -e 'source(\"code/utils/gating_tree.R\"); write_gating_mermaid(read_gating_tree(), \"docs/gating_tree.md\")'")

if (n_fail > 0) quit(status = 1)
cat("all tests passed\n")
