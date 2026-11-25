# =============================================================================
# gating_tree.R — read and apply config/gating_tree.yaml (rules: see that file)
# =============================================================================

GATE_KEYS <- c("pos", "neg", "any", "rest")

# Read the gating tree and check that it is well formed. Stops on any error.
read_gating_tree <- function(path = file.path(getwd(), "config/gating_tree.yaml")) {
  cfg   <- yaml::read_yaml(path)
  nodes <- gating_nodes(cfg)
  errs  <- character()
  if (anyDuplicated(nodes$id)) errs <- c(errs, paste("duplicated node id:", nodes$id[duplicated(nodes$id)]))
  gate_markers <- unique(unlist(strsplit(nodes$markers[nodes$markers != ""], ",")))
  miss <- setdiff(gate_markers, names(cfg$markers))
  if (length(miss)) errs <- c(errs, paste("gate marker not in `markers`:", miss))
  ann_markers <- vapply(cfg$annotations, `[[`, "", "marker")
  miss <- setdiff(ann_markers, names(cfg$markers))
  if (length(miss)) errs <- c(errs, paste("annotation marker not in `markers`:", miss))
  labs <- c(cfg$default_label, unique(nodes$label[!is.na(nodes$label)]))
  miss <- setdiff(labs, names(cfg$labels))
  if (length(miss)) errs <- c(errs, paste("label without colour in `labels`:", miss))
  bad <- nodes$id[nodes$gate == ""]
  if (length(bad)) errs <- c(errs, paste("node without gate:", bad))
  if (length(errs)) stop("config/gating_tree.yaml:\n  ", paste(errs, collapse = "\n  "))
  cfg
}

# Flatten the tree: one row per node (id, parent, depth, gate text, markers, label, png, mfi).
gating_nodes <- function(cfg) {
  rows <- list()
  walk <- function(nodes, parent, depth) {
    for (nd in nodes) {
      g <- nd$gate
      if (!is.null(g) && length(setdiff(names(g), GATE_KEYS))) {
        stop("node ", nd$id, ": unknown gate key ", paste(setdiff(names(g), GATE_KEYS), collapse = ", "))
      }
      rows[[length(rows) + 1]] <<- data.frame(
        id      = nd$id,
        parent  = parent,
        depth   = depth,
        gate    = gate_text(g),
        markers = paste(c(g$pos, g$neg, g$any), collapse = ","),
        label   = if (is.null(nd$label)) NA_character_ else nd$label,
        png     = if (is.null(nd$png)) NA_character_ else nd$png,
        mfi     = paste(nd$mfi, collapse = ","),
        has_children = length(nd$children) > 0,
        dev_diagram  = if (is.null(nd$deviation$diagram)) NA_character_ else nd$deviation$diagram,
        dev_reason   = if (is.null(nd$deviation$reason)) NA_character_ else nd$deviation$reason,
        stringsAsFactors = FALSE
      )
      walk(nd$children, nd$id, depth + 1)
    }
  }
  walk(cfg$tree, "root", 1)
  do.call(rbind, rows)
}

# Labels given by a node and by all its descendants, e.g. every immune label:
# gating_labels_under(cfg, "immune"). Used by steps 07/07b instead of hand-written lists.
gating_labels_under <- function(cfg, ids) {
  nodes <- gating_nodes(cfg)
  miss  <- setdiff(ids, nodes$id)
  if (length(miss)) stop("gating_labels_under: no node ", paste(miss, collapse = ", "))
  keep <- ids
  repeat {
    more <- union(keep, nodes$id[nodes$parent %in% keep])
    if (length(more) == length(keep)) break
    keep <- more
  }
  unique(nodes$label[nodes$id %in% keep & !is.na(nodes$label)])
}

gate_text <- function(g) {
  if (is.null(g)) return("")
  if (isTRUE(g$rest)) return("rest")
  parts <- c(if (length(g$pos)) paste0(g$pos, "+"),
             if (length(g$neg)) paste0(g$neg, "-"),
             if (length(g$any)) paste0("(", paste0(g$any, "+", collapse = " or "), ")"))
  paste(parts, collapse = " ")
}

# Gate path from the root to a node, e.g. "CD45+ > CD20- > CD3+ CD56-".
gating_path <- function(nodes, id) {
  out <- character()
  while (id != "root") {
    i   <- match(id, nodes$id)
    out <- c(nodes$gate[i], out)
    id  <- nodes$parent[i]
  }
  paste(out, collapse = " > ")
}

# Every difference from the reference diagram, from the YAML: node `deviation` fields,
# `general_deviations` and `not_implemented` (written by step 06 as 06_gating_deviations.tsv).
gating_deviations <- function(cfg) {
  nd  <- gating_nodes(cfg)
  dv  <- nd[!is.na(nd$dev_reason), ]
  rows <- data.frame(
    item           = ifelse(is.na(dv$label), dv$id, dv$label),
    node           = dv$id,
    diagram        = dv$dev_diagram,
    implementation = vapply(dv$id, function(i) gating_path(nd, i), ""),
    reason         = dv$dev_reason,
    stringsAsFactors = FALSE
  )
  gen <- do.call(rbind, lapply(cfg$general_deviations, function(x) data.frame(
    item = x$item, node = NA_character_, diagram = x$diagram, implementation = x$implementation,
    reason = x$reason, stringsAsFactors = FALSE)))
  ni  <- do.call(rbind, lapply(cfg$not_implemented, function(x) data.frame(
    item = x$png, node = NA_character_, diagram = x$png, implementation = "not implemented",
    reason = x$reason, stringsAsFactors = FALSE)))
  out <- rbind(rows, gen, ni)
  rownames(out) <- NULL
  out
}

# Apply the tree. pos = named list (logical marker name -> TRUE/FALSE/NA per cell; NA = flagged,
# counted as negative). Returns one label per cell.
apply_gating_tree <- function(cfg, pos, n) {
  p <- function(m) {
    v <- pos[[m]]
    if (is.null(v)) stop("gating marker without positivity: ", m)
    !is.na(v) & v
  }
  fb <- cfg$immune_fallback
  celltype <- rep(cfg$default_label, n)

  node_cells <- function(nd, parent, taken) {
    g <- nd$gate
    if (isTRUE(g$rest)) return(parent & !taken)
    if (!is.null(fb) && identical(nd$id, fb$node)) {
      n_na <- sum(is.na(pos[[fb$marker]]))
      if (n_na > fb$max_na_fraction * n) {
        message(sprintf("  WARNING: %s NA in %.0f%% of cells - using %s as immune proxy",
                        fb$marker, 100 * n_na / n, paste(fb$proxy_any, collapse = "/")))
        return(parent & Reduce(`|`, lapply(fb$proxy_any, p)))
      }
    }
    s <- parent
    for (m in g$pos) s <- s & p(m)
    for (m in g$neg) s <- s & !p(m)
    if (length(g$any)) s <- s & Reduce(`|`, lapply(g$any, p))
    s
  }
  walk <- function(nodes, parent) {
    taken <- rep(FALSE, n)
    for (nd in nodes) {
      s <- node_cells(nd, parent, taken)
      if (!is.null(nd$label)) celltype[s] <<- nd$label
      if (length(nd$children)) walk(nd$children, s)
      taken <- taken | s
    }
  }
  walk(cfg$tree, rep(TRUE, n))
  celltype
}

# Mermaid flowchart of the tree (docs/gating_tree.md).
write_gating_mermaid <- function(cfg, path) {
  nd  <- gating_nodes(cfg)
  esc <- function(x) gsub('"', "'", x)
  box <- function(i) {
    lab <- nd$id[i]
    if (!is.na(nd$label[i]) && nd$label[i] != nd$id[i]) lab <- paste0(lab, " → ", nd$label[i])
    if (!is.na(nd$png[i])) lab <- paste0(lab, "<br/><i>diagram: ", nd$png[i], "</i>")
    sprintf('  %s["%s"]', nd$id[i], esc(lab))
  }
  edges <- sprintf('  %s -->|"%s"| %s', nd$parent, esc(nd$gate), nd$id)
  ni <- vapply(cfg$not_implemented, function(x) sprintf("| %s | %s |", x$png, x$reason), "")
  dv <- gating_deviations(cfg)
  dv <- dv[dv$implementation != "not implemented", ]
  cell <- function(x) gsub("|", "\\|", ifelse(is.na(x), "", x), fixed = TRUE)
  dl <- sprintf("| %s | %s | %s | %s |", cell(dv$item), cell(dv$diagram), cell(dv$implementation), cell(dv$reason))
  lines <- c(
    "# Gating tree",
    "",
    "Generated from `config/gating_tree.yaml` by `write_gating_mermaid()` (code/utils/gating_tree.R): do not edit by hand.",
    "Edge text = gate applied to the cells of the parent (`+` positive, `-` negative, `rest` = cells not taken by an earlier sibling).",
    "",
    "```mermaid",
    "flowchart LR",
    '  root["all cells"]',
    vapply(seq_len(nrow(nd)), box, ""),
    edges,
    "```",
    "",
    "## Populations of the diagram not implemented",
    "",
    "| Diagram | Reason |",
    "|---|---|",
    ni,
    "",
    "## Deviations from the diagram",
    "",
    "Same rows as `06_gating_deviations.tsv` (step 06). Implementation = gates from the root (`>` = child of).",
    "",
    "| Item | Diagram | Implementation | Reason |",
    "|---|---|---|---|",
    dl,
    ""
  )
  writeLines(lines, path)
  invisible(path)
}
