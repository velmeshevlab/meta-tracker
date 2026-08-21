# Trajectory graph plotting.

#' @export
plot_combined_graph <- function(cds, reduction_method = "UMAP", color_cells_by = NULL, alpha = 0.4) {
  library(igraph)
  library(ggplot2)
  library(dplyr)
  
  # Extract all graphs and merge into one
  graph_list <- cds@graphs
  combined_graph <- Reduce(igraph::union, graph_list)
  
  # Get graph node names
  graph_nodes <- V(combined_graph)$name
  
  # Get graph node coordinates from dp_mst
  node_coords <- as.data.frame(t(cds@principal_graph_aux[[reduction_method]]$dp_mst))
  node_coords$node <- rownames(node_coords)
  colnames(node_coords)[1:2] <- c("x", "y")
  
  node_coords <- node_coords %>% filter(node %in% graph_nodes)
  
  # Get edge list with coordinates
  el <- igraph::get.edgelist(combined_graph)
  edges_df <- data.frame(
    from = el[,1],
    to = el[,2],
    x = node_coords$x[match(el[,1], node_coords$node)],
    y = node_coords$y[match(el[,1], node_coords$node)],
    xend = node_coords$x[match(el[,2], node_coords$node)],
    yend = node_coords$y[match(el[,2], node_coords$node)]
  )
  
  # Get cell coordinates
  cell_coords <- as.data.frame(reducedDims(cds)[[reduction_method]], stringsAsFactors = FALSE)
  colnames(cell_coords) <- c("x", "y")
  
  # Prepare cell coordinates with optional color
  cell_coords$cell_id <- rownames(cell_coords)
  if (!is.null(color_cells_by) && color_cells_by %in% colnames(colData(cds))) {
    cell_coords$color <- colData(cds)[cell_coords$cell_id, color_cells_by]
  } else {
    cell_coords$color <- NA
  }
  
  # Plot
  ggplot() +
    geom_point(data = cell_coords, aes(x = x, y = y, color = color), size = 0.01, alpha = alpha) +
    geom_segment(data = edges_df, aes(x = x, y = y, xend = xend, yend = yend), color = "black", size = 1) +
    geom_point(data = node_coords, aes(x = x, y = y), color = "blue", size = 0.2) +
    theme_minimal() +
    coord_fixed()
}
#' Plot cells of an isolated lineage over the full embedding extent
#'
#' Draws the cells (and trajectory graph) of a single isolated lineage with
#' \code{monocle3::plot_cells}, but framed to the axis limits of the *full*
#' object so lineages are directly comparable across plots.
#'
#' @param cds Full \code{metatracker_data_set}: supplies the embedding limits,
#'   and (unless \code{cds_sub} is given) the lineage to extract.
#' @param lineage Name of the lineage to plot (e.g. "VIP"). Ignored when
#'   \code{cds_sub} is supplied; required otherwise.
#' @param cds_sub Optional pre-extracted lineage object. If \code{NULL}, it is
#'   obtained with \code{get_lineage_object(cds, lineage)}.
#' @param reduction Reduced-dimension name used for the axis limits (default "UMAP").
#' @param color_cells_by Column passed to \code{plot_cells} (default "pseudotime").
#' @param cell_size,graph_label_size,trajectory_graph_color,trajectory_graph_segment_size
#'   Passed through to \code{plot_cells}.
#' @return A \code{ggplot} object.
#' @export
plot_lineage_cells <- function(cds, lineage = NULL, cds_sub = NULL,
                               reduction = "UMAP",
                               color_cells_by = "pseudotime",
                               cell_size = 0.1,
                               graph_label_size = 1.5,
                               trajectory_graph_color = "cyan",
                               trajectory_graph_segment_size = 1.5) {
  emb <- reducedDims(cds)[[reduction]]
  if (is.null(emb))
    stop("Reduction '", reduction, "' not found in reducedDims(cds).")
  xlim <- range(emb[, 1], na.rm = TRUE)
  ylim <- range(emb[, 2], na.rm = TRUE)

  if (is.null(cds_sub)) {
    if (is.null(lineage))
      stop("Provide either `lineage` (to extract) or `cds_sub` (already extracted).")
    cds_sub <- get_lineage_object(cds, lineage)
  }

  umap_theme <- theme(
    plot.title       = element_blank(),
    legend.position  = "none",
    panel.border     = element_blank(),
    axis.text.x      = element_text(size = 16),
    axis.text.y      = element_text(size = 16),
    axis.title.x     = element_text(size = 18, face = "bold"),
    axis.title.y     = element_text(size = 18, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black"),
    panel.background = element_blank(),
    legend.title     = element_text(size = 16))

  plot_cells(cds_sub,
             color_cells_by                = color_cells_by,
             label_cell_groups             = FALSE,
             label_leaves                  = FALSE,
             label_branch_points           = FALSE,
             graph_label_size              = graph_label_size,
             cell_size                     = cell_size,
             trajectory_graph_color        = trajectory_graph_color,
             trajectory_graph_segment_size = trajectory_graph_segment_size) +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    umap_theme
}


# ============================================================================
# Accessors for the post-compression cds structure, used by the plotters below.
#   @expression[[lin]] : data.frame (method "sum") OR list(sum=, mean=) ("all")
#   @expectation[[lin]]: matrix (grid x genes), gene column names
#   @pseudotime[[lin]] : list(real=, scaled=)
#   @lineages[[lin]]   : list(name=, updated_pt=)  (see .lineage_cells)
# ============================================================================
.cmp_expr_sum <- function(cds, lineage) {
  e <- cds@expression[[lineage]]
  if (is.list(e) && !is.null(e[["sum"]])) e[["sum"]] else e
}
.cmp_expr_mean <- function(cds, lineage) {
  e <- cds@expression[[lineage]]
  if (is.list(e) && !is.null(e[["mean"]])) e[["mean"]] else e
}
.cmp_pt_real <- function(cds, lineage) {
  p <- cds@pseudotime[[lineage]]
  if (is.list(p) && !is.null(p[["real"]])) p[["real"]]
  else if (is.data.frame(p) || is.matrix(p)) p[, 1] else p
}
.cmp_pt_scaled <- function(cds, lineage) {
  p <- cds@pseudotime[[lineage]]
  if (is.list(p) && !is.null(p[["scaled"]])) p[["scaled"]]
  else if (is.data.frame(p) || is.matrix(p)) p[, 1] else p
}
# Fitted (expectation) values for one gene; NULL if the gene is absent.
.cmp_fit <- function(cds, lineage, gene) {
  em <- cds@expectation[[lineage]]
  if (is.matrix(em) || is.data.frame(em)) {
    if (gene %in% colnames(em)) as.numeric(em[, gene]) else NULL
  } else if (is.list(em) && !is.null(em[[gene]])) {   # legacy nested-list layout
    p <- em[[gene]]
    if (is.list(p) && !is.null(p[["prediction"]])) as.numeric(p[["prediction"]]) else as.numeric(p)
  } else NULL
}

# Minimal monocle-style theme (ggplot2 is imported, so bare names resolve).
.monocle_theme_opts <- function() {
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(linewidth = 0.25, color = "black")) +
    theme(axis.line.y = element_line(linewidth = 0.25, color = "black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.key = element_blank())
}

# Apply FUN over sliding windows of `data` (evobiR::SlidingWindow equivalent).
.sliding_window <- function(FUN, data, window, step) {
  FUN    <- match.fun(FUN)
  total  <- length(data)
  window <- max(1L, floor(window))
  step   <- max(1L, floor(step))
  spots  <- seq(1, max(1, total - window), by = step)
  vapply(spots, function(s) FUN(data[s:(s + window - 1)]), numeric(1))
}

#' Overlay expression/fit curves for a gene across lineages
#'
#' Plots per-lineage fitted expectation curves (and optionally the meta-cell
#' expression points) for one gene, on a shared pseudotime axis. Works with the
#' post-compression cds structure (\code{@expression} sum/mean, list-form
#' \code{@pseudotime}, matrix \code{@expectation}).
#'
#' @param cds A compressed \code{metatracker_data_set}.
#' @param gene Gene to plot.
#' @param lineages Character vector of lineage names to overlay.
#' @param meta Optional metadata data.frame (rows = cells) for \code{age.scale}.
#' @param points If TRUE, draw meta-cell expression points under each fit.
#' @param age.scale,scale.lineage,age.points,breaks.labels Optional age-axis mapping.
#' @param point_size,line_size,text.size,plot.title.size,legend.key.size,legend.text.size Sizes.
#' @param colors Per-lineage colours (recycled in order).
#' @param N Unused legacy arg (grid length is taken from the data).
#' @param legend_position ggplot legend position.
#' @return A ggplot object.
#' @export
plot_multiple <- function(cds, gene, lineages, meta = NULL, points = TRUE,
                          age.scale = FALSE, scale.lineage = NULL,
                          age.points = c("3rd trimester", "0-1 years", "2-4 years", "4-10 years"),
                          breaks.labels = c("2nd", "3rd", "birth", "1y", "4y"),
                          point_size = 0.1, line_size = 1, text.size = 14,
                          plot.title.size = 36, legend.key.size = 0.5,
                          legend.text.size = 10,
                          colors = c("red","blue","green","cyan","magenta","purple","orange","black","yellow","tan"),
                          N = 500, legend_position = "none"){
  N <- nrow(.cmp_expr_sum(cds, lineages[1]))
  if (length(scale.lineage) == 0) {
    max.pt <- max(unlist(lapply(lineages, function(l) .cmp_pt_real(cds, l))))
  } else {
    max.pt <- max(.cmp_pt_real(cds, scale.lineage))
  }
  step_pt <- max.pt / (N - 1)
  if (isTRUE(points)) {
    dd <- data.frame(pseudotime = seq(0, max.pt, by = step_pt))
    fits <- c()
    for (lineage in lineages) {
      exp_df <- .cmp_expr_sum(cds, lineage)
      if (gene %in% colnames(exp_df)) {
        expv <- as.numeric(exp_df[, gene])
        fitv <- .cmp_fit(cds, lineage, gene); if (is.null(fitv)) fitv <- rep(0, N)
      } else { expv <- rep(0, N); fitv <- rep(0, N) }
      dd[[paste0("exp_", lineage)]] <- expv
      dd[[paste0("fit_", lineage)]] <- fitv
      fits <- c(fits, fitv)
    }
    ymax <- max(fits)
  } else {
    fits <- c(); rows <- list()
    for (lineage in lineages) {
      fitv <- .cmp_fit(cds, lineage, gene); if (is.null(fitv)) fitv <- rep(0, N)
      fits <- c(fits, fitv)
      rows[[lineage]] <- data.frame(pseudotime = seq(0, max.pt, by = step_pt),
                                    fit = fitv, lineage = lineage)
    }
    dd <- do.call(rbind, rows)
    dd$lineage <- factor(dd$lineage, levels = lineages)
    ymax <- max(fits)
  }
  q <- ggplot(data = dd)
  if (isTRUE(points)) {
    for (M in seq_along(lineages)) {
      q <- q +
        geom_point(aes_string(x = "pseudotime", y = paste0("exp_", lineages[M]),
                              color = "pseudotime"), size = I(point_size)) +
        scale_color_gradient2(lineages[M], low = "grey", high = colors[M]) +
        ggnewscale::new_scale_color() +
        geom_line(aes_string(x = "pseudotime", y = paste0("fit_", lineages[M])),
                  linewidth = I(line_size), color = colors[M])
    }
  } else {
    q <- q + geom_line(aes(x = pseudotime, y = fit, color = lineage),
                       linewidth = I(line_size)) +
      scale_color_manual(values = colors)
  }
  q <- q + scale_y_log10()
  if (isTRUE(age.scale)) {
    if (length(scale.lineage) == 1) {
      cells <- .lineage_cells(cds@lineages[[scale.lineage]])
      age <- meta[cells, c("age_num", "age_range")]
    } else {
      age <- meta[, c("age_num", "age_range")]
    }
    age <- age[order(age$age_num), ]
    window <- nrow(age) / N
    step   <- (nrow(age) - window) / N
    age.comp <- .sliding_window("mean", age$age_num, window, step)
    d <- cbind(as.data.frame(seq(0, max.pt, by = step_pt)), age.comp)
    breaks.list <- c(0)
    for (age.point in age.points) {
      age.break <- stats::quantile(age[age$age_range == age.point, ]$age_num, 0.95)
      age.break <- d[which.min(abs(d[, 2] - age.break)), 1]
      breaks.list <- append(breaks.list, age.break)
    }
    q <- q + scale_x_continuous(breaks = breaks.list, labels = breaks.labels)
  }
  q <- q + ylim(c(0, ymax)) + .monocle_theme_opts() +
    ylab("Expression") + xlab("Pseudotime") + ggtitle(gene) +
    theme(legend.key.size = grid::unit(legend.key.size, "cm"),
          plot.title = element_text(size = plot.title.size, face = "bold", hjust = 0.5),
          axis.text = element_text(size = text.size),
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.title = element_blank(),
          legend.text = element_text(size = legend.text.size),
          legend.title = element_text(size = text.size, face = "bold"),
          legend.position = legend_position)
  q
}

#' Plot a gene across all lineages (faceted or overlaid)
#'
#' For every lineage in \code{cds@lineages}, plots the mean meta-cell expression
#' points and the fitted expectation curve for \code{gene}. Uses the mean matrix
#' (\code{@expression[[lin]]$mean}) and matrix \code{@expectation}; lineages
#' missing the gene are skipped.
#'
#' @param cds A compressed \code{metatracker_data_set} (use \code{method = "all"}
#'   so the mean matrix is available).
#' @param gene Gene to plot.
#' @param overlay If FALSE (default) facet by lineage; if TRUE overlay on one panel.
#' @param custom_lineage_colors Optional named vector overriding lineage colours.
#' @return A ggplot object.
#' @export
plot_all <- function(cds, gene, overlay = FALSE, custom_lineage_colors = NULL) {
  lineages <- names(cds@lineages)
  n_lin    <- length(lineages)
  numeric_suffix <- suppressWarnings(as.numeric(sub("^[^0-9]*", "", lineages)))
  ordered_lineages <- if (!all(is.na(numeric_suffix))) lineages[order(numeric_suffix)] else lineages

  base_colors  <- RColorBrewer::brewer.pal(9, "Blues")[-c(1, 2, 3)]
  lineage_cols <- grDevices::colorRampPalette(base_colors)(n_lin)
  names(lineage_cols) <- ordered_lineages
  if (!is.null(custom_lineage_colors)) {
    for (lin in names(custom_lineage_colors))
      if (lin %in% lineages) lineage_cols[lin] <- custom_lineage_colors[lin]
  }
  lineage_cols_points <- rev(lineage_cols); names(lineage_cols_points) <- ordered_lineages

  first_lin <- lineages[1]
  N <- length(.cmp_pt_scaled(cds, first_lin)); if (N == 0) N <- 100
  pt_grid <- seq(0, 1, length.out = N)

  df_list <- list()
  for (lin in lineages) {
    fit <- .cmp_fit(cds, lin, gene)
    if (is.null(fit) || all(is.na(fit))) next
    fit <- log(fit + 1)
    mean_df <- .cmp_expr_mean(cds, lin)
    if (!gene %in% colnames(mean_df)) next
    expr <- log(as.numeric(mean_df[[gene]]) + 1)
    pt   <- .cmp_pt_scaled(cds, lin)
    keep <- !is.na(expr) & !is.na(pt)
    if (sum(keep) == 0) next
    df_points <- data.frame(gene = gene, lineage = lin, pt = pt[keep],
                            expr = expr[keep], type = "raw", stringsAsFactors = FALSE)
    df_fit    <- data.frame(gene = gene, lineage = lin, pt = pt_grid,
                            expr = fit, type = "fit", stringsAsFactors = FALSE)
    df_list[[lin]] <- rbind(df_points, df_fit)
  }
  if (length(df_list) == 0) stop("No valid data found for gene: ", gene)
  df_all <- do.call(rbind, df_list)
  df_all$lineage <- factor(df_all$lineage, levels = ordered_lineages)

  df_points <- subset(df_all, type == "raw")
  df_points$color_mirror <- lineage_cols_points[as.character(df_points$lineage)]

  p <- ggplot(df_all, aes(x = pt, y = expr)) +
    geom_point(data = df_points, aes(color = lineage),
               color = df_points$color_mirror, alpha = 0.8, size = 1) +
    geom_line(data = subset(df_all, type == "fit"), aes(color = lineage), linewidth = 2)
  if (!overlay) p <- p + facet_wrap(~ lineage, scales = "free_y")
  p <- p + scale_color_manual(values = lineage_cols) +
    theme_classic(base_size = 12) +
    theme(legend.position = "right", strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
    labs(title = gene, x = "Pseudotime", y = "log(Expression + 1)", color = "Lineage")
  p
}
