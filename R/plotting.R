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
