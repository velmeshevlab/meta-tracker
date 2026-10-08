densify_lineage_graph <- function(cds, lineage, spacing = NULL, factor = 1.5,
                                  reduction_method = "UMAP", reselect_cells = FALSE,
                                  update_principal_graph = FALSE,
                                  N = 5, cl = 1, sel_clusters = NULL,
                                  start_regions = FALSE, starting_clusters = FALSE) {
  
  if (is.null(cds@graphs[[lineage]]))
    stop("No graph for lineage '", lineage, "' in cds@graphs.", call. = FALSE)
  
  g_sub <- cds@graphs[[lineage]]
  Y     <- cds@principal_graph_aux[[reduction_method]]$dp_mst      # dims x nodes
  vs    <- V(g_sub)$name
  miss  <- setdiff(vs, colnames(Y))
  if (length(miss) > 0)
    stop("Lineage nodes missing from dp_mst: ", paste(miss, collapse = ", "), call. = FALSE)
  
  coords <- t(Y[, vs, drop = FALSE])                              # nodes x dims
  el     <- get.edgelist(g_sub)                                   # E x 2 (names)
  if (nrow(el) == 0) stop("Lineage graph '", lineage, "' has no edges.", call. = FALSE)
  
  edge_len <- sqrt((coords[el[, 1], 1] - coords[el[, 2], 1])^2 +
                     (coords[el[, 1], 2] - coords[el[, 2], 2])^2)
  if (is.null(spacing)) spacing <- stats::median(edge_len)
  thresh <- factor * spacing
  
  # next free Y_<k> index across the WHOLE dp_mst (keeps new names globally unique)
  all_idx <- suppressWarnings(as.integer(sub("^Y_", "", colnames(Y))))
  next_id <- max(all_idx[is.finite(all_idx)], 0L) + 1L
  
  new_names   <- character(0)
  new_coords  <- matrix(numeric(0), nrow = 0, ncol = 2)
  new_edge_df <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
  
  for (e in seq_len(nrow(el))) {
    a <- el[e, 1]; b <- el[e, 2]; L <- edge_len[e]
    n_seg <- if (L > thresh) max(2L, as.integer(round(L / spacing))) else 1L
    if (n_seg == 1L) {
      new_edge_df <- rbind(new_edge_df, data.frame(from = a, to = b, stringsAsFactors = FALSE))
      next
    }
    A <- coords[a, ]; B <- coords[b, ]
    fracs <- seq_len(n_seg - 1L) / n_seg
    ins   <- paste0("Y_", next_id + seq_len(n_seg - 1L) - 1L)
    next_id <- next_id + (n_seg - 1L)
    xy <- cbind(A[1] + fracs * (B[1] - A[1]),
                A[2] + fracs * (B[2] - A[2]))
    rownames(xy) <- ins
    new_names   <- c(new_names, ins)
    new_coords  <- rbind(new_coords, xy)
    chain       <- c(a, ins, b)
    new_edge_df <- rbind(new_edge_df,
                         data.frame(from = utils::head(chain, -1),
                                    to   = utils::tail(chain, -1),
                                    stringsAsFactors = FALSE))
  }
  
  n_subdiv <- sum(edge_len > thresh)
  if (length(new_names) == 0) {
    message(sprintf("No edges exceeded %.3g (= %.3g x spacing %.3g); nothing to densify.",
                    thresh, factor, spacing))
    return(cds)
  }
  message(sprintf("Densified lineage '%s': added %d node(s) across %d edge(s) (spacing %.3g).",
                  lineage, length(new_names), n_subdiv, spacing))
  
  # ---- add new node coordinates to the shared dp_mst ----
  add_cols <- t(new_coords)                                       # dims x new
  colnames(add_cols) <- new_names
  rownames(add_cols) <- rownames(Y)
  Y_new <- cbind(Y, add_cols)
  cds@principal_graph_aux[[reduction_method]]$dp_mst <- Y_new
  
  # ---- rebuild the densified lineage subgraph ----
  all_nodes <- c(vs, new_names)
  node_df   <- data.frame(name = all_nodes,
                          x = Y_new[1, all_nodes], y = Y_new[2, all_nodes],
                          stringsAsFactors = FALSE)
  cds@graphs[[lineage]] <- graph_from_data_frame(new_edge_df, vertices = node_df,
                                                 directed = FALSE)
  
  # ---- optionally keep the global principal graph consistent ----
  if (isTRUE(update_principal_graph)) {
    G     <- cds@principal_graph[[reduction_method]]
    G_el  <- get.edgelist(G)
    subp  <- el[edge_len > thresh, , drop = FALSE]                # original long edges
    is_sub <- rep(FALSE, nrow(G_el))
    for (i in seq_len(nrow(subp))) {
      a <- subp[i, 1]; b <- subp[i, 2]
      is_sub <- is_sub | (G_el[, 1] == a & G_el[, 2] == b) | (G_el[, 1] == b & G_el[, 2] == a)
    }
    kept        <- G_el[!is_sub, , drop = FALSE]
    chain_edges <- as.matrix(new_edge_df[new_edge_df$from %in% new_names |
                                           new_edge_df$to   %in% new_names, c("from", "to")])
    G_el_new <- rbind(kept, chain_edges)
    all_v    <- union(V(G)$name, new_names)
    vdf      <- data.frame(name = all_v, x = Y_new[1, all_v], y = Y_new[2, all_v],
                           stringsAsFactors = FALSE)
    cds@principal_graph[[reduction_method]] <-
      graph_from_data_frame(as.data.frame(G_el_new, stringsAsFactors = FALSE),
                            vertices = vdf, directed = FALSE)
  }
  
  # ---- re-select cells along the densified graph (fills the gap) ----
  if (isTRUE(reselect_cells)) {
    cds <- isolate_lineage(cds, lineage, sel_clusters = sel_clusters,
                           start_regions = start_regions,
                           starting_clusters = starting_clusters,
                           subset = FALSE, N = N, cl = cl)
  }
  cds
}


isolate_lineage <- function(cds, lineage, sel_clusters = NULL, start_regions = F, starting_clusters = F, subset = FALSE, N = 5, cl = 1, r = r){
  sel.cells = isolate_lineage_sub(cds, lineage, sel_clusters = sel_clusters, start_regions = start_regions, starting_clusters = starting_clusters, subset = subset, N = N, cl = cl, r = r)
  cds@lineages[[lineage]] <- sel.cells
  return(cds)
}


isolate_lineage_sub <- function(cds, lineage, sel_clusters = NULL, start_regions = NULL, starting_clusters = NULL,
                                excluded_ages = NULL, excluded_ages_cluster = NULL, subset = FALSE, N = 5, cl = 1, r){
  sub.graph = cds@graphs[[lineage]]
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  if(subset == F){
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names(V(sub.graph))]))
  }
  else{
    g = principal_graph(cds)[["UMAP"]]
    dd = degree(g)
    names1 = names(dd[dd > 2 | dd == 1])
    names2 = names(dd[dd == 2])
    names2 = sample(names2, length(names2)/subset, replace = F)
    names = c(names1, names2)
    names = intersect(names(V(sub.graph)), names)
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names]))
  }
  #select cells along the graph
  #mean.dist = path.distance(nodes_UMAP.sub)
  r = r
  cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
  colnames(cells_UMAP) <- toupper(colnames(cells_UMAP))
  cells_UMAP = cells_UMAP[,c("UMAP_1", "UMAP_2")]
  sel.cells = cell.selector(nodes_UMAP.sub, cells_UMAP, r, cl = cl)
  #only keep cells in the progenitor and lineage-specific clusters
  sel.cells1 = c()
  sel.cells2 = sel.cells
  if(length(starting_clusters) > 0){
    sel.cells1 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% starting_clusters])
  }
  if(length(start_regions) > 0){
    sel.cells1 = sel.cells1[sel.cells1 %in% rownames(cds@colData[cds@colData$region %in% start_regions,])]
  }
  if(length(sel_clusters) > 0){
    sel.cells2 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% sel_clusters])
  }
  
  # Optional: exclude cells in all selected clusters with specific age labels, if no cluster specified, then exclude the cells
  # in all clusters with the specific age labels
  if (length(excluded_ages) > 0) {
    cd <- colData(cds)
    
    # decide which age column the labels belong to (same check as before)
    if (all(excluded_ages %in% cd$age_reorder)) {
      age_col <- "age_reorder"
    } else if (all(excluded_ages %in% cd$age_details)) {
      age_col <- "age_details"
    } else {
      stop("excluded_ages not all found in a single column. Missing from age_reorder: ",
           paste(setdiff(excluded_ages, cd$age_reorder), collapse = ", "),
           " | missing from age_details: ",
           paste(setdiff(excluded_ages, cd$age_details), collapse = ", "))
    }
    
    cluster_labels <- cds@"clusters"[["UMAP"]]$clusters
    
    in_cluster <- if (length(excluded_ages_cluster) > 0) {
      as.character(cluster_labels) %in% as.character(excluded_ages_cluster)
    } else {
      rep(TRUE, length(cluster_labels))     # no cluster given: apply to all cells
    }
    
    cells_to_exclude <- names(cluster_labels)[
      which(in_cluster & cd[names(cluster_labels), age_col] %in% excluded_ages)
    ]
    sel.cells2 <- sel.cells2[!(sel.cells2 %in% cells_to_exclude)]
  }
  cells = unique(c(sel.cells1, sel.cells2))
  sel.cells = sel.cells[sel.cells %in% cells]
  return(sel.cells)
}


plot_combined_graph <- function(cds, reduction_method = "UMAP", color_cells_by = NULL, alpha = 0.4,
                                label_clusters = FALSE, highlight_clusters = NULL,
                                other_color = "grey85", show_others = TRUE,
                                highlight_lineage = NULL, highlight_color = "orange",
                                highlight_size = 0.3, highlight_alpha = 0.6,
                                highlight_on_top = TRUE) {
  library(igraph)
  library(ggplot2)
  library(dplyr)
  
  graph_list <- cds@graphs
  combined_graph <- Reduce(igraph::union, graph_list)
  graph_nodes <- V(combined_graph)$name
  
  node_coords <- as.data.frame(t(cds@principal_graph_aux[[reduction_method]]$dp_mst))
  node_coords$node <- rownames(node_coords)
  colnames(node_coords)[1:2] <- c("x", "y")
  node_coords <- node_coords %>% filter(node %in% graph_nodes)
  
  el <- igraph::get.edgelist(combined_graph)
  edges_df <- data.frame(
    from = el[,1], to = el[,2],
    x    = node_coords$x[match(el[,1], node_coords$node)],
    y    = node_coords$y[match(el[,1], node_coords$node)],
    xend = node_coords$x[match(el[,2], node_coords$node)],
    yend = node_coords$y[match(el[,2], node_coords$node)]
  )
  
  cell_coords <- as.data.frame(reducedDims(cds)[[reduction_method]], stringsAsFactors = FALSE)
  colnames(cell_coords) <- c("x", "y")
  cell_coords$cell_id <- rownames(cell_coords)
  has_color <- !is.null(color_cells_by) && color_cells_by %in% colnames(colData(cds))
  if (has_color) {
    cell_coords$color <- as.character(colData(cds)[cell_coords$cell_id, color_cells_by])
  } else {
    cell_coords$color <- NA_character_
  }
  
  if (!is.null(highlight_clusters) && !has_color)
    stop("highlight_clusters needs color_cells_by (the column that holds the cluster labels)")
  
  highlight_coords <- NULL
  if (!is.null(highlight_lineage)) {
    if (!highlight_lineage %in% names(cds@lineages))
      stop("lineage '", highlight_lineage, "' not found in cds@lineages")
    lin <- cds@lineages[[highlight_lineage]]
    lin_cells <- if (is.list(lin)) lin$name else lin
    highlight_coords <- cell_coords[cell_coords$cell_id %in% lin_cells, ]
    if (nrow(highlight_coords) == 0) stop("none of the lineage cells were found in the UMAP")
  }
  
  hl_layer <- NULL
  if (!is.null(highlight_coords)) {
    hl_layer <- geom_point(data = highlight_coords, aes(x = x, y = y),
                           color = highlight_color, size = highlight_size, alpha = highlight_alpha)
  }
  
  p <- ggplot()
  if (!highlight_on_top) p <- p + hl_layer
  
  if (!has_color) {
    # plain grey cells, no color aesthetic, so there is no "NA" legend
    p <- p + geom_point(data = cell_coords, aes(x = x, y = y),
                        color = other_color, size = 0.01, alpha = alpha)
    label_df <- NULL
  } else if (is.null(highlight_clusters)) {
    p <- p + geom_point(data = cell_coords, aes(x = x, y = y, color = color), size = 0.01, alpha = alpha)
    label_df <- cell_coords
  } else {
    is_hl <- cell_coords$color %in% as.character(highlight_clusters)
    if (!any(is_hl)) stop("none of highlight_clusters found in ", color_cells_by)
    if (show_others) {
      p <- p + geom_point(data = cell_coords[!is_hl, ], aes(x = x, y = y),
                          color = other_color, size = 0.01, alpha = alpha)
    }
    p <- p + geom_point(data = cell_coords[is_hl, ], aes(x = x, y = y, color = color),
                        size = 0.01, alpha = alpha)
    label_df <- cell_coords[is_hl, ]
  }
  
  if (highlight_on_top) p <- p + hl_layer          # above the cells, below the trajectory
  
  p <- p +
    geom_segment(data = edges_df, aes(x = x, y = y, xend = xend, yend = yend), color = "black", size = 1) +
    geom_point(data = node_coords, aes(x = x, y = y), color = "blue", size = 0.2) +
    theme_minimal() +
    coord_fixed()
  
  if (label_clusters && !is.null(label_df)) {
    centroids <- label_df %>%
      group_by(color) %>%
      summarise(x = median(x), y = median(y), .groups = "drop")
    p <- p + ggrepel::geom_text_repel(data = centroids, aes(x = x, y = y, label = color),
                                      size = 5, color = "black", fontface = "bold",
                                      bg.color = "white", bg.r = 0.15, max.overlaps = Inf)
  }
  
  p
}

get_lineage_object <- function(cds, lineage = FALSE, N = FALSE, recalculate_pt = TRUE){
  start = find_start_node(cds)
  if (lineage != FALSE) {
    sub.graph <- cds@graphs[[lineage]]
    if (is.list(cds@lineages[[lineage]])) {
      sel.cells <- cds@lineages[[lineage]]$name
    } else {
      sel.cells <- cds@lineages[[lineage]]
    }
    if (!is.character(sel.cells)) {
      print("sel cells are not string")
    }
  }
  else{
    sel.cells = colnames(cds)
  }
  sel.cells = sel.cells[sel.cells %in% colnames(cds)]
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  if(N != FALSE){
    if(N < length(sel.cells)){
      sel.cells = sample(sel.cells, N)
    }
  }
  #subset the moncole object
  cds_subset = cds[,sel.cells]
  #set the graph, node and cell UMAP coordinates
  if(lineage == FALSE){
    sub.graph = principal_graph(cds_subset)[["UMAP"]]
  }
  nodes_UMAP <- nodes_UMAP[,names(V(sub.graph))]
  #Reorder the vertices
  degrees <- igraph::degree(sub.graph)
  endpoints <- names(degrees[degrees == 1])
  path_result <- igraph::shortest_paths(sub.graph, from = start, to = endpoints[endpoints != start])
  path_names <- names(path_result$vpath[[1]])   # ordered root -> tip
  # Build the new sequential, zero-padded names
  n <- length(path_names)
  width <- nchar(as.character(n))
  new_names <- paste0("Y_", formatC(seq_len(n), width = width, flag = "0"))
  # Create the old-name -> new-name mapping
  rename_map <- setNames(new_names, path_names)
  igraph::V(sub.graph)$old_name <- igraph::V(sub.graph)$name
  igraph::V(sub.graph)$name <- rename_map[igraph::V(sub.graph)$name]
  #colnames(nodes_UMAP) <- rename_map[colnames(nodes_UMAP)]
  new_colnames <- unname(rename_map[colnames(nodes_UMAP)])
  colnames(nodes_UMAP) <- new_colnames
  cds_subset@principal_graph[["UMAP"]] <- sub.graph
  cds_subset@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP
  cds_subset@clusters[["UMAP"]]$partitions <- cds_subset@clusters[["UMAP"]]$partitions[colnames(cds_subset)]
  #recalculate closest vertex and pseudotime for the selected cells
  if(recalculate_pt == TRUE){
    source_url("https://raw.githubusercontent.com/cole-trapnell-lab/monocle3/master/R/learn_graph.R")
    cds_subset <- project2MST(cds_subset, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP)
    cds_subset <- order_cells(cds_subset, root_pr_nodes = unname(rename_map[start]))
  }
  return(cds_subset)
}


compress_lineage <- function(cds, lineage, n_metacells = 1000,
                             pt_name = "updated_pt") {
  
  lin <- cds@lineages[[lineage]]
  if (!is.list(lin) || is.null(lin$name) || is.null(lin[[pt_name]]))
    stop("lineage '", lineage, "' needs a list with $name and $", pt_name)
  
  pt_vec <- lin[[pt_name]]
  if (is.null(names(pt_vec))) {
    if (length(pt_vec) != length(lin$name)) stop("unnamed pseudotime does not match $name in length")
    names(pt_vec) <- lin$name
  }
  
  # cells that exist in the object and have a pseudotime
  cells <- lin$name[lin$name %in% colnames(cds)]
  cells <- cells[cells %in% names(pt_vec) & !is.na(pt_vec[cells])]
  if (length(cells) < n_metacells) warning(lineage, ": fewer cells than metacells (", length(cells), ")")
  
  cds_subset = cds[,cells]
  #preprare raw count matrix
  exp = as.data.frame(as.matrix(exprs(cds_subset)))
  exp_sum <- t(exp)
  exp_sum = exp_sum[,rownames(cds)]
  #prepare size factor
  size_factor <- (pData(cds_subset)[, 'Size_Factor'])
  size_factor <- size_factor[rownames(exp_sum)]
  #prepare pseudotime
  #pt <- cds@lineages[[lineage]][['updated_pt']]
  #pt <- cds@lineages[[lineage]]$gap_reduced_pt
  pt <- as.data.frame(pt_vec)
  pt <- pt[rownames(exp_sum), , drop = FALSE]
  colnames(pt) <- c("pseudotime")
  #prepare umap
  UMAP <- reducedDims(cds_subset)[["UMAP"]]
  UMAP <- UMAP[rownames(exp_sum),]
  exp_sum = cbind(pt, UMAP, size_factor, exp_sum)
  exp_sum = exp_sum[order(exp_sum$pseudotime),]
  exp_sum$meta_cell <- cut(rank(exp_sum$pseudotime), breaks = n_metacells, labels = FALSE)
  gene_cols <- setdiff(
    colnames(exp_sum),
    c("pseudotime", "umap_1", "umap_2", "size_factor", "meta_cell")
  )
  print(paste0("Compressing lineage ", lineage, " with sum"))
  meta_sum <- exp_sum %>%
    group_by(meta_cell) %>%
    summarise(
      n_cells = n(),                     # number of cells in this meta-cell
      pseudotime = mean(pseudotime),     # mean pseudotime
      umap_1 = mean(umap_1),             # mean UMAP_x
      umap_2 = mean(umap_2),             # mean UMAP_y
      size_factor = sum(size_factor),    # sum size factors
      across(all_of(gene_cols), sum)   # sum gene counts
    ) %>%
    ungroup()
  meta_sum_ordered <- meta_sum[order(meta_sum$pseudotime), ]
  
  exp_mean <- exp
  exp_mean = (t(exp_mean)) /  (pData(cds_subset)[, 'Size_Factor'])
  exp_mean = exp_mean[,rownames(cds)]
  pt <- pt[rownames(exp_mean), , drop = FALSE]
  UMAP <- UMAP[rownames(exp_mean),]
  exp_mean = cbind(pt, UMAP, exp_mean)
  exp_mean = exp_mean[order(exp_mean$pseudotime),]
  exp_mean$meta_cell <- cut(rank(exp_mean$pseudotime), breaks = n_metacells, labels = FALSE)
  gene_cols <- setdiff(
    colnames(exp_mean),
    c("pseudotime", "umap_1", "umap_2", "meta_cell")
  )
  print(paste0("Compressing lineage ", lineage, " with mean"))
  meta_mean <- exp_mean %>%
    group_by(meta_cell) %>%
    summarise(
      n_cells = n(),                     # number of cells in this meta-cell
      pseudotime = mean(pseudotime),     # mean pseudotime
      umap_1 = mean(umap_1),             # mean UMAP_x
      umap_2 = mean(umap_2),             # mean UMAP_y
      across(all_of(gene_cols), mean)   # mean gene counts
    ) %>%
    ungroup()
  meta_mean_ordered <- meta_mean[order(meta_mean$pseudotime), ]
  list(sum = meta_sum_ordered, mean = meta_mean_ordered)
}
compress_all_lineages <- function(cds, lineages = names(cds@lineages), save_dir = NULL, ...) {
  if (!is.null(save_dir)) dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (lineage in lineages) {
    message("Compressing lineage ", lineage)
    res <- tryCatch(compress_lineage(cds, lineage, ...),
                    error = function(e) { warning("skipping ", lineage, ": ", conditionMessage(e)); NULL })
    if (is.null(res)) next
    
    if (is.null(cds@expression[[lineage]])) cds@expression[[lineage]] <- list()
    cds@expression[[lineage]]$sum  <- res$sum
    cds@expression[[lineage]]$mean <- res$mean
    
    if (!is.null(save_dir))
      saveRDS(res, file.path(save_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", lineage), "_meta.rds")))
    rm(res); gc()
  }
  cds
}
