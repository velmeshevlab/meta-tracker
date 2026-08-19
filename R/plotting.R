umap_theme <- theme(plot.title = element_blank(), legend.position="none", panel.border = element_blank(), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title.x = element_text(size=18, face="bold"), axis.title.y = element_text(size=18, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.title=element_text(size=16))

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
