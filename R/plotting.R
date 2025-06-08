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
