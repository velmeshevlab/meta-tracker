umap_theme <- theme(plot.title = element_blank(), legend.position="none", panel.border = element_blank(), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title.x = element_text(size=18, face="bold"), axis.title.y = element_text(size=18, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.title=element_text(size=16))

find_start_point <-function(graph_list){
  start_ends = c()
  for(graph_name in names(graph_list)){
    graph = graph_list[[graph_name]]
    start_end = V(graph)[degree(graph) == 1]$name
    start_ends = c(start_ends, start_end)
  }
  start = names(sort(table(start_ends),decreasing=TRUE)[1])
  start
}

get_subgraph_opposite_to_start <- function(graph, branch_vertex, start) {
    # Remove the branch vertex
    g_split <- delete_vertices(graph, branch_vertex)
    # Get connected components
    comps <- components(g_split)
    # Identify component containing `start`
    start_component <- comps$membership[start]
    # Get vertices not in the same component as `start`
    other_vertices <- names(comps$membership[comps$membership != start_component])
    # Return subgraph of those vertices
    subgraph = induced_subgraph(g_split, vids = other_vertices)
    subgraph
}

get_branch_points <- function(graph_list){
  combined_graph <- graph.empty(directed = FALSE)
  for (g in graph_list) {
    combined_graph <- igraph::union(combined_graph, g)
  }
  combined_graph <- simplify(combined_graph)
  branch_points <- V(combined_graph)[degree(combined_graph) >= 3]
  branch_points = branch_points$name
  branch_points
}

group_graphs_by_vertex_overlap <- function(graph_list) {
  graph_names <- names(graph_list)
  
  if (is.null(graph_names)) {
    graph_names <- paste0("graph_", seq_along(graph_list))
    names(graph_list) <- graph_names
  }
  
  # Ensure all graphs have vertex names
  for (i in seq_along(graph_list)) {
    g <- graph_list[[i]]
    if (is.null(V(g)$name)) {
      V(g)$name <- as.character(V(g))
    }
  }
  
  # Prepare vertex sets
  vertex_sets <- lapply(graph_list, function(g) V(g)$name)
  
  # Build overlap matrix
  n <- length(graph_list)
  overlap_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(graph_names, graph_names))
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (length(intersect(vertex_sets[[i]], vertex_sets[[j]])) > 0) {
        overlap_matrix[i, j] <- 1
        overlap_matrix[j, i] <- 1
      }
    }
  }
  
  # Create graph of overlaps and find components
  overlap_graph <- graph_from_adjacency_matrix(overlap_matrix, mode = "undirected", diag = FALSE)
  comps <- components(overlap_graph)
  
  # Group graph names by component
  grouped_graphs <- split(names(comps$membership), comps$membership)
  
  return(grouped_graphs)
}

find_branches <- function(cds){
graph_list = cds@graphs
start = find_start_point(graph_list)
branch_points = get_branch_points(graph_list)
branch_list_full = list()
for(branch_vertex in branch_points){
  graph_list_f = list()
  i <- 1
  for(graph in graph_list){
    if(branch_vertex %in% V(graph)$name)
      graph_list_f[[names(graph_list)[i]]] <- graph
      i <- i + 1
  }
  graph_list_trunc = lapply(graph_list_f, get_subgraph_opposite_to_start, start = start, branch_vertex = branch_vertex)
  branch_list = group_graphs_by_vertex_overlap(graph_list_trunc)
  names(branch_list) <- c("B1", "B2")
  branch_list_full[[branch_vertex]] <- branch_list
}
combined_graph <- graph.empty(directed = FALSE)
for (g in graph_list) {
  combined_graph <- igraph::union(combined_graph, g)
}
combined_graph <- simplify(combined_graph)
distances <- sapply(branch_points, function(bp) {
  sp <- suppressWarnings(shortest.paths(combined_graph, v = start, to = bp))
  return(sp[1, 1])
})
branch_points_sorted <- branch_points[order(distances)]
branch_list_full = branch_list_full[branch_points_sorted]
names = c()
for(i in 1:length(branch_list_full)){
  names = c(names, paste0("BP_", i))
  i <- i + 1
}
names(branch_list_full) <- names
branch_list_full
}
                        
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

#' @export
plot_multiple <- function(cds, gene, lineages, meta = NULL, points = T, age.scale = F, scale.lineage = NULL, age.points = c("3rd trimester", "0-1 years", "2-4 years", "4-10 years"), breaks.labels = c("2nd", "3rd", "birth", "1y", "4y"), point_size = 0.1, line_size = 1, text.size = 14, plot.title.size = 36, legend.key.size = 0.5, legend.text.size = 10, colors = c("red", "blue", "green", "cyan", "magenta", "purple", "orange", "black", "yellow", "tan"), N = 500, legend_position = "none"){
  cds_name = deparse(substitute(cds))
  input = paste0(cds_name,"@expression$", lineages[1])
  N = nrow(eval(parse(text = input)))
  pts = c()
  if(length(scale.lineage) == 0){
  for(lineage in lineages){
    input = paste0(cds_name,"@pseudotime$", lineage)
    pt = eval(parse(text = input))[,1]
    pts = c(pts, pt)
  }
  max.pt = max(pts)
  }
  else{
  input = paste0(cds_name,"@pseudotime$", scale.lineage)
  pt = eval(parse(text = input))[,1]
  max.pt = max(pt)
  }
  print(max.pt)
  if(points == T){
    dd = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    cols = c("pseudotime")
    fits = c()
    exps = c()
    for(lineage in lineages){
      input = paste0("exp = ",cds_name,"@expression$", lineage)
      eval(parse(text=input))
      if(gene %in% colnames(exp)){
        input = paste0("exp = ",cds_name,"@expression$", lineage,"[,'",gene,"']")
        eval(parse(text=input))
        input = paste0("fit = ",cds_name,"@expectation$", lineage,"[,'",gene,"']")
        eval(parse(text=input))
      }
      else{
        exp = rep(0, N)
        fit = rep(0, N)
      }
      dd = cbind(dd, as.numeric(exp), as.numeric(fit))
      cols = append(cols, paste0("exp_", lineage))
      cols = append(cols, paste0("fit_", lineage))
      fits = c(fits, fit)
      exps = c(exps, exp)
    }
    colnames(dd) <- cols
    ymax = max(fits)
  }
  else{
    fits = c()
    dd = matrix(ncol = 3, nrow = 0,)
    for(lineage in lineages){
      input = paste0("exp = ",cds_name,"@expression$", lineage)
      eval(parse(text=input))
      if(gene %in% colnames(exp)){
        input = paste0("fit = ",cds_name,"@expectation$", lineage,"[,'",gene,"']")
        eval(parse(text=input))
      }
      else{
        fit = rep(0, N)
      }
      fits = c(fits, fit)
      dd = rbind(dd, cbind(seq(from=0, to=max.pt, by = max.pt/(N-1)), fit, rep(lineage, length(fit))))
    }
    ymax = max(fits)
    colnames(dd) <- c("pseudotime", "fit", "lineage")
    dd = as.data.frame(dd)
    dd$pseudotime <- as.numeric(dd$pseudotime)
    dd$fit <- as.numeric(dd$fit)
    dd$lineage <- factor(dd$lineage, levels = lineages)
  }
  q <- ggplot(data = dd)
  if(points == T){
    for(M in 1:length(lineages)){
      loop_input1 = paste0("geom_point(aes_string(x='pseudotime',y = '", paste0('exp_', lineages[M]), "',color='pseudotime'), size=I(", point_size, "))")
      loop_input2 = paste0("scale_color_gradient2(lineages[M],low='grey', ", "high='",colors[M],"')")
      loop_input3 = "new_scale_color()"
      loop_input4 = paste0("geom_line(aes_string(x='pseudotime', y = '", paste0('fit_', lineages[M]), "',size = I(", line_size, ")), color = '", colors[M],"')")
      q <- q + eval(parse(text=loop_input1)) + eval(parse(text=loop_input2)) + eval(parse(text=loop_input3)) + eval(parse(text=loop_input4))
    }
  }
  else{
    q <- q + geom_line(aes(x = pseudotime, y = fit, color = lineage), size = I(line_size)) + scale_color_manual(values = colors)
  }
  q <- q + scale_y_log10()
  if(age.scale == T){
    if(length(scale.lineage) == 1){
    input = paste0(cds_name,"@lineages$", scale.lineage)
    cells = eval(parse(text = input))
    age = meta[cells,c("age_num", "age_range")]
    }
    else{
    age = meta[,c("age_num", "age_range")]
    }
    age = age[order(age$age_num),]
    window = nrow(age)/N
    step = ((nrow(age)-window)/N)
    age.comp = SlidingWindow("mean", age$age_num, window, step)
    d = seq(from=0, to=max.pt, by = max.pt/(N-1))
    d = cbind(as.data.frame(d), age.comp) 
    breaks.list = c(0)
    for(age.point in age.points){
      age.break = quantile(age[age$age_range == age.point,]$age_num, 0.95)
      age.break = d[which.min(abs(d[,2]-age.break)),1]
      breaks.list = append(breaks.list, age.break)
    }
    q <- q + scale_x_continuous(breaks = breaks.list, labels = breaks.labels)
  }
  q <- q + ylim(y = c(0,ymax))
  q <- q + monocle_theme_opts() + ylab("Expression") + xlab("Pseudotime") + ggtitle(gene) + theme(legend.key.size = unit(legend.key.size, 'cm'), plot.title = element_text(size = plot.title.size, face="bold", hjust = 0.5), axis.text=element_text(size=text.size), axis.text.x=element_text(angle = 60, hjust=1), axis.title=element_blank(), legend.text=element_text(size=legend.text.size), legend.title=element_text(size=text.size, face = "bold"), legend.position = legend_position)
  q
}


plot_all <- function(cds, gene, overlay = FALSE, custom_lineage_colors = NULL) {
  library(ggplot2)
  library(RColorBrewer)
  library(colorspace)  # for darken()
  
  lineages <- names(cds@lineages)
  n_lin <- length(lineages)
  
  # Order lineages numerically if they have numeric suffixes
  numeric_suffix <- suppressWarnings(as.numeric(sub("^[^0-9]*", "", lineages)))
  if (!all(is.na(numeric_suffix))) {
    ordered_lineages <- lineages[order(numeric_suffix)]
  } else {
    ordered_lineages <- lineages
  }
  
  # Colors
  base_colors <- brewer.pal(9, "Blues")[-c(1,2,3)]
  lineage_cols <- colorRampPalette(base_colors)(n_lin)
  names(lineage_cols) <- ordered_lineages
  
  if (!is.null(custom_lineage_colors)) {
    for (lin in names(custom_lineage_colors)) {
      if (lin %in% lineages) lineage_cols[lin] <- custom_lineage_colors[lin]
    }
  }
  
  lineage_cols_lines <- darken(lineage_cols, amount = 0.3)
  lineage_cols_points <- rev(lineage_cols)
  names(lineage_cols_points) <- ordered_lineages
  
  # Grid for fitted lines
  first_lin <- lineages[1]
  N <- if (!is.null(cds@pseudotime[[first_lin]]$scaled)) length(cds@pseudotime[[first_lin]]$scaled) else 100
  pt_grid <- seq(0, 1, length.out = N)
  
  df_list <- list()
  
  for (lin in lineages) {
    
    # Skip if gene missing
    if (is.null(cds@expectation[[lin]][[gene]])) next
    
    # Fitted values
    fit <- cds@expectation[[lin]][[gene]][["prediction"]]
    if (is.null(fit) || all(is.na(fit))) next
    fit <- log(fit + 1)
    
    # Raw data
    expr <- cds@expression[[lin]]$mean[[gene]]
    expr <- log(expr + 1)
    
    pt <- cds@pseudotime[[lin]]$scaled
    keep <- !is.na(expr) & !is.na(pt)
    if (sum(keep) == 0) next
    
    df_points <- data.frame(
      gene = gene,
      lineage = lin,
      pt = pt[keep],
      expr = expr[keep],
      type = "raw",
      stringsAsFactors = FALSE
    )
    
    df_fit <- data.frame(
      gene = gene,
      lineage = lin,
      pt = pt_grid,
      expr = fit,
      type = "fit",
      stringsAsFactors = FALSE
    )
    
    df_list[[lin]] <- rbind(df_points, df_fit)
  }
  
  if (length(df_list) == 0) stop("No valid data found for gene: ", gene)
  df_all <- do.call(rbind, df_list)
  
  df_all$lineage <- factor(df_all$lineage, levels = ordered_lineages)
  
  # Plot
  p <- ggplot(df_all, aes(x = pt, y = expr))
  df_points <- subset(df_all, type == "raw")
  df_points$color_mirror <- lineage_cols_points[df_points$lineage]
  
  p <- p +
    geom_point(
      data = df_points,
      aes(color = lineage),
      color = df_points$color_mirror,
      alpha = 0.8,
      size = 1
    ) +
    geom_line(
      data = subset(df_all, type == "fit"),
      aes(color = lineage),
      linewidth = 2
    )
  
  if (!overlay) {
    p <- p + facet_wrap(~ lineage, scales = "free_y")
  }
  
  p <- p +
    scale_color_manual(values = lineage_cols) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = gene,
      x = "Pseudotime",
      y = "log(Expression + 1)",
      color = "Lineage"
    )
  
  return(p)
}
