# Annotate and plot branch points / branches on the graph.

# --- Convenience wrapper ------------------------------------------------------
#' @export
label_branches <- function(cds, file = "branches.pdf",
                           drop_non_binary = TRUE, ...) {
  annot <- .annotate_branches(cds, drop_non_binary = drop_non_binary)
  if (nrow(annot$arm_table) == 0) {
    warning("No binary branch points found to label.", call. = FALSE)
  }
  .plot_labeled_branches(annot, file = file, ...)
  message("Wrote ", file)
  invisible(annot)
}

# --- Core annotation ----------------------------------------------------------
.annotate_branches <- function(cds, drop_non_binary = TRUE) {
  graph_list <- cds@graphs
  start      <- .find_start_point(graph_list)
  bps        <- .get_branch_points(graph_list)

  combined <- .build_combined_graph(graph_list)

  # Order branch points by distance from start (matches .find_branches).
  distances <- sapply(bps, function(bp) {
    sp <- suppressWarnings(shortest.paths(combined, v = start, to = bp))
    sp[1, 1]
  })
  bps <- bps[order(distances)]

  arm_rows   <- list()
  bp_label_of <- setNames(rep("", length(bps)), bps)  # vertex -> label
  k <- 0
  for (bv in bps) {
    keep_idx <- vapply(graph_list, function(g) bv %in% V(g)$name, logical(1))
    gf       <- graph_list[keep_idx]
    gt       <- lapply(gf, .get_subgraph_opposite_to_start,
                       start = start, branch_vertex = bv)
    arms     <- .group_graphs_by_vertex_overlap(gt)

    if (length(arms) != 2 && drop_non_binary) next
    k <- k + 1
    lab <- paste0("BP_", k)
    bp_label_of[bv] <- lab

    for (a in seq_along(arms)) {
      lineages_in_arm <- arms[[a]]
      verts <- unique(unlist(lapply(lineages_in_arm,
                                    function(ln) V(gt[[ln]])$name)))
      arm_rows[[length(arm_rows) + 1]] <- data.frame(
        bp_label      = lab,
        branch_vertex = bv,
        arm           = paste0(lab, "_B", a),
        lineages      = paste(lineages_in_arm, collapse = ","),
        vertices      = paste(verts, collapse = ","),
        stringsAsFactors = FALSE)
    }
  }

  arm_table <- if (length(arm_rows)) do.call(rbind, arm_rows) else
    data.frame(bp_label = character(), branch_vertex = character(),
               arm = character(), lineages = character(),
               vertices = character(), stringsAsFactors = FALSE)

  # Write vertex attributes onto the combined graph.
  V(combined)$is_branch_point <- V(combined)$name %in% bps
  V(combined)$bp_label <- bp_label_of[V(combined)$name]
  V(combined)$bp_label[is.na(V(combined)$bp_label)] <- ""

  layout <- .make_layout(combined, graph_list, start)

  list(combined_graph = combined,
       layout         = layout,
       start          = start,
       branch_points  = bps,
       arm_table      = arm_table,
       drop_non_binary = drop_non_binary)
}

# --- Plotting -----------------------------------------------------------------
.plot_labeled_branches <- function(annot, file = NULL, which_bp = NULL,
                                  label_arms = TRUE,
                                  label_arms_on_overview = FALSE,
                                  vertex_label = FALSE,
                                  vertex_size = 4, width = 8, height = 8) {
  g      <- annot$combined_graph
  nlay   <- .norm_layout(annot$layout)      # normalized; used for nodes + text
  start  <- annot$start
  vnames <- V(g)$name
  bp_labels <- unique(annot$arm_table$bp_label)
  if (!is.null(which_bp)) bp_labels <- intersect(bp_labels, which_bp)

  if (!is.null(file)) {
    grDevices::pdf(file, width = width, height = height, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  # a little headroom so edge labels aren't clipped
  lim <- c(-1.18, 1.18)
  arm_palette <- c("#D55E00", "#0072B2")     # colour-blind safe orange / blue
  arm_text_col <- c("#8A3B00", "#004E7A")    # darker, for readable text

  .plot_base <- function(vcol, vsz, vlab, main) {
    plot(g, layout = nlay, rescale = FALSE, xlim = lim, ylim = lim,
         vertex.color = vcol, vertex.size = vsz, vertex.frame.color = NA,
         vertex.label = vlab, vertex.label.cex = 0.85,
         vertex.label.color = "black", vertex.label.dist = 0.7,
         vertex.label.font = 2,
         edge.color = "grey65", edge.width = 1, main = main)
  }

  # ---- Overview: mark every branch point ----
  is_bp <- nchar(V(g)$bp_label) > 0
  vcol <- rep("grey80", length(vnames))
  vcol[is_bp] <- "red"; vcol[vnames == start] <- "forestgreen"
  vsz <- rep(vertex_size, length(vnames))
  vsz[is_bp] <- vertex_size * 2.4; vsz[vnames == start] <- vertex_size * 2.4
  vlab <- rep(NA_character_, length(vnames))
  vlab[is_bp] <- V(g)$bp_label[is_bp]; vlab[vnames == start] <- "start"
  if (vertex_label) vlab <- vnames

  .plot_base(vcol, vsz, vlab, "Trajectory: branch points")

  if (label_arms_on_overview) {
    for (lab in bp_labels) {
      rows <- annot$arm_table[annot$arm_table$bp_label == lab, , drop = FALSE]
      bv   <- rows$branch_vertex[1]
      for (p in .arm_label_positions(rows, nlay, vnames, bv)) {
        text(p$x, p$y, labels = p$arm, cex = 0.7, font = 2, col = "grey20")
      }
    }
  }

  # ---- One page per branch point: colour + label its two arms ----
  for (lab in bp_labels) {
    rows <- annot$arm_table[annot$arm_table$bp_label == lab, , drop = FALSE]
    bv   <- rows$branch_vertex[1]

    vcol <- rep("grey85", length(vnames))
    for (a in seq_len(nrow(rows))) {
      verts <- strsplit(rows$vertices[a], ",", fixed = TRUE)[[1]]
      vcol[vnames %in% verts] <- arm_palette[((a - 1) %% length(arm_palette)) + 1]
    }
    vcol[vnames == start] <- "forestgreen"; vcol[vnames == bv] <- "black"

    vsz <- rep(vertex_size, length(vnames))
    vsz[vnames == bv] <- vertex_size * 2.6; vsz[vnames == start] <- vertex_size * 2.2

    vlab <- rep(NA_character_, length(vnames))
    vlab[vnames == bv] <- lab
    if (vertex_label) vlab <- vnames

    .plot_base(vcol, vsz, vlab, paste0(lab, "  (vertex ", bv, ")"))

    # arm labels drawn directly on the graph, at each arm's tip
    positions <- .arm_label_positions(rows, nlay, vnames, bv)
    for (a in seq_along(positions)) {
      p <- positions[[a]]
      col_a <- arm_text_col[((a - 1) %% length(arm_text_col)) + 1]
      text(p$x, p$y, labels = p$arm, cex = 0.9, font = 2, col = col_a)
    }

    if (label_arms) {
      legend("topleft", bty = "n",
             legend = paste0(rows$arm, ":  ", rows$lineages),
             fill = arm_palette[seq_len(nrow(rows))], cex = 0.8)
    }
  }
  invisible(annot)
}

# --- Build the union graph the same way .find_branches does -------------------
.build_combined_graph <- function(graph_list) {
  combined <- graph.empty(directed = FALSE)
  for (g in graph_list) combined <- igraph::union(combined, g)
  simplify(combined)
}

# --- Coordinate detection: return a name -> c(x, y) map, or NULL -------------
# Checks the per-lineage graphs (not the union, whose attrs may be suffixed)
# for common coordinate attribute name pairs.
.detect_coord_map <- function(graph_list) {
  candidate_pairs <- list(c("x", "y"), c("X", "Y"),
                          c("coord_x", "coord_y"),
                          c("umap_1", "umap_2"), c("UMAP_1", "UMAP_2"),
                          c("dim1", "dim2"))
  coord <- list()
  for (g in graph_list) {
    va <- vertex_attr_names(g)
    hit <- Filter(function(p) all(p %in% va), candidate_pairs)
    if (length(hit) == 0) return(NULL)     # a graph without coords -> give up
    p  <- hit[[1]]
    nm <- V(g)$name
    xs <- as.numeric(vertex_attr(g, p[1]))
    ys <- as.numeric(vertex_attr(g, p[2]))
    for (k in seq_along(nm)) coord[[nm[k]]] <- c(xs[k], ys[k])
  }
  if (length(coord) == 0) return(NULL)
  coord
}

# --- A layout matrix aligned to V(combined)$name ------------------------------
.make_layout <- function(combined, graph_list, start, seed = 1L) {
  vnames <- V(combined)$name
  cmap   <- .detect_coord_map(graph_list)

  if (!is.null(cmap) && all(vnames %in% names(cmap))) {
    lay <- t(vapply(vnames, function(n) cmap[[n]], numeric(2)))
    return(lay)
  }

  # No usable coordinates: compute a deterministic layout.
  set.seed(seed)
  root_idx <- which(vnames == start)
  lay <- tryCatch(
    layout_as_tree(combined, root = if (length(root_idx)) root_idx else 1L),
    error = function(e) NULL)
  if (is.null(lay)) lay <- layout_with_kk(combined)
  lay
}

# Normalize a layout to [-1, 1] on both axes so we can plot with
# rescale = FALSE and place text() at the same coordinates the vertices use.
.norm_layout <- function(lay) {
  f <- function(v) {
    r <- range(v, na.rm = TRUE)
    if (diff(r) == 0) rep(0, length(v)) else 2 * (v - r[1]) / diff(r) - 1
  }
  cbind(f(lay[, 1]), f(lay[, 2]))
}

# Position + label for each arm: placed at the arm tip (vertex furthest from
# the branch vertex), nudged a little further outward so it clears the node.
.arm_label_positions <- function(rows, nlay, vnames, bv, nudge = 0.10) {
  bpos <- nlay[which(vnames == bv), , drop = TRUE]
  out  <- list()
  for (a in seq_len(nrow(rows))) {
    verts <- strsplit(rows$vertices[a], ",", fixed = TRUE)[[1]]
    idx   <- which(vnames %in% verts)
    if (length(idx) == 0) next
    d    <- sqrt((nlay[idx, 1] - bpos[1])^2 + (nlay[idx, 2] - bpos[2])^2)
    tip  <- idx[which.max(d)]
    pos  <- nlay[tip, ]
    dir  <- pos - bpos
    nrm  <- sqrt(sum(dir^2))
    if (nrm > 0) pos <- pos + nudge * dir / nrm
    out[[length(out) + 1]] <- list(arm = rows$arm[a], x = pos[1], y = pos[2])
  }
  out
}
