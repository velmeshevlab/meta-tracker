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
combined_graph <- graph.empty(directed = FALSE)
for (g in graph_list) {
  combined_graph <- igraph::union(combined_graph, g)
}
combined_graph <- simplify(combined_graph)
distances <- sapply(branch_points, function(bp) {
  sp <- suppressWarnings(shortest.paths(combined_graph, v = start, to = bp))
  return(sp[1, 1])
})
branch_points <- branch_points[order(distances)]
branch_list_full = list()
for(j in 1:length(branch_points)){
  branch_vertex = branch_points[j]
  graph_list_f = list()
  i <- 1
  for(graph in graph_list){
    if(branch_vertex %in% V(graph)$name)
      graph_list_f[[names(graph_list)[i]]] <- graph
      i <- i + 1
  }
  graph_list_trunc = lapply(graph_list_f, get_subgraph_opposite_to_start, start = start, branch_vertex = branch_vertex)
  branch_list = group_graphs_by_vertex_overlap(graph_list_trunc)
  names(branch_list) <- c(paste0("BP", j, "_B1"), paste0("BP", j, "_B2"))
  branch_list_full[[branch_vertex]] <- branch_list
}
names = c()
for(i in 1:length(branch_list_full)){
  names = c(names, paste0("BP_", i))
  i <- i + 1
}
names(branch_list_full) <- names
branch_list_full
}
                        
format_branch_specific_genes <- function(branch_point, cds, branch_number = 1, p_cutoff = 0.05, FC_cutoff = 0.2, dynamic_p_cutoff = 0.05, dynamic_FC_cutoff = 0.1, p_adjust = "BH", dynamic_test = "Moran"){
  lineages = names(cds@lineages)
  branches_1 = branch_point[[branch_number]]
  if(branch_number == 1){
    branches_2 = branch_point[[2]]
  }
  else{
    branches_2 = branch_point[[1]]
  }
  branch_genes = list()
  branch_gene_names = list()
  for(lineage in branches_1){
    lineage_genes = cds@lineage_genes[[lineage]]
    FC_names = c()
    for(lin in branches_2){
      FC_names = c(FC_names, paste0("FC_", lineage, "vs", lin))
    }
    p_names = c()
    for(lin in branches_2){
      p_names = c(p_names, paste0("pvalue_", lineage, "vs", lin))
    }
    if(dynamic_FC_cutoff != F){
    dynamic = cds@dynamic_genes[[lineage]]
    if(dynamic_test == "Moran"){
      dynamic_genes = rownames(dynamic[dynamic$I >= dynamic_FC_cutoff & dynamic$padj < dynamic_p_cutoff, ])
    }
    else{
      dynamic_genes = rownames(dynamic[dynamic$scaled_FC >= dynamic_FC_cutoff, ])
    }
    lineage_genes = lineage_genes[dynamic_genes,]
    }
    lineage_genes = lineage_genes[,c(p_names, FC_names, "meta_p", "average_FC")]
    FCs = lineage_genes[,(length(branches_2)+1):(ncol(lineage_genes)-2)]
    if(is.vector(FCs)){
      lineage_spec_genes = lineage_genes[FCs >= FC_cutoff & lineage_genes$meta_p < p_cutoff, ]
    }
    else{
      lineage_spec_genes = lineage_genes[rowSums(FCs >= FC_cutoff) == ncol(FCs) & lineage_genes$meta_p < p_cutoff, ]
    }
    branch_genes[[lineage]] <- lineage_spec_genes
    branch_gene_names[[lineage]] <- rownames(lineage_spec_genes)
  }
  branch_gene_names = Reduce(intersect, branch_gene_names)
  meta_p_matrix = matrix(0,nrow = length(branch_gene_names), ncol = 0)
  median_FC_matrix = matrix(0,nrow = length(branch_gene_names), ncol = 0)
  for(lineage in branches_1){
    lineage_genes = branch_genes[[lineage]]
    FC_names = c()
    for(lin in branches_2){
      FC_names = c(FC_names, paste0("FC_", lineage, "vs", lin))
    }
    FCs = lineage_genes[branch_gene_names, FC_names]
    meta_p = lineage_genes[branch_gene_names,]$meta_p
    meta_p_matrix = cbind(meta_p_matrix, meta_p)
    if(is.vector(FCs)){
      median_FC = FCs
    }
    else{
      median_FC = apply(FCs, 1, median)
    }
    median_FC_matrix = cbind(median_FC_matrix, median_FC)
  }
  if(is.vector(FCs)){
    meta_p = meta_p_matrix
    median_FC = median_FC_matrix
  }
  else{
    meta_p = apply(meta_p_matrix, 1, get_meta_p)
    median_FC = apply(median_FC_matrix, 1, median)
  }
  if(p_adjust != FALSE){
    meta_p <- p.adjust(meta_p, method = p_adjust)
  }
  final_out = cbind(meta_p, median_FC)
  rownames(final_out) <- branch_gene_names
  final_out = as.data.frame(final_out)
  final_out$branch <- rep(names(branch_point)[branch_number], nrow(final_out))
  colnames(final_out) <- c("meta_p", "median_FC", "branch")
  final_out = final_out[with(final_out, order(-abs(median_FC), meta_p)), ]
  final_out
}

format_lineage_specific_genes <- function(lineage, cds, p_cutoff = 0.05, FC_cutoff = 0.2, dynamic_FC_cutoff = 0.1, dynamic_p_cutoff = 0.05, p_adjust = "BH", specificity = "high", dynamic_test = "Moran"){
  lineages = names(cds@lineages)
  lineage_genes = cds@lineage_genes[[lineage]]
  if(dynamic_FC_cutoff != F){
    dynamic = cds@dynamic_genes[[lineage]]
    if(dynamic_test == "Moran"){
      dynamic_genes = rownames(dynamic[dynamic$I >= dynamic_FC_cutoff & dynamic$padj < dynamic_p_cutoff, ])
    }
    else{
      dynamic_genes = rownames(dynamic[dynamic$scaled_FC >= dynamic_FC_cutoff, ])
    }
    lineage_genes = lineage_genes[dynamic_genes,]
  }
  p_values = lineage_genes[,1:length(lineages)]
  if(p_adjust != FALSE){
    p_values_adj = apply(p_values, 2, p.adjust, method = p_adjust)
  }
  if(p_adjust != FALSE){
    lineage_genes$meta_p <- p.adjust(lineage_genes$meta_p, method = p_adjust)
  }
  FCs = lineage_genes[,(length(lineages)):(ncol(lineage_genes)-2)]
  median_FC = apply(FCs, 1, median)
  lineage_genes$median_FC <- median_FC
  if(specificity == "high"){
    lineage_spec_genes = lineage_genes[rowSums(FCs >= FC_cutoff) == ncol(FCs) & lineage_genes$meta_p < p_cutoff, ]
  }
  else{
    lineage_spec_genes = lineage_genes[lineage_genes$meta_p < p_cutoff & lineage_genes$median_FC >= FC_cutoff, ]
  }
  lineage_spec_genes = lineage_spec_genes[with(lineage_spec_genes, order(-abs(median_FC), meta_p)), ]
  lineage_spec_genes$lineage <- rep(lineage, nrow(lineage_spec_genes))
  lineage_spec_genes
}

lineage_specific_genes_v2 <- function(test_lineage, cds, pattern, diffend, genes = NULL, lineages){
  print(paste0("Testing ", test_lineage))
  if(length(lineages) > 2){
    index = which(test_lineage == lineages)
    p_list = c()
    p_names = c()
    wd_list = c()
    wd_names = c()
    fc_names = c()
    for(lineage in lineages){
      index2 = which(lineage == lineages)
      if(index != index2){
        if(index<index2){
          p_name = paste0("pvalue_", index, "vs", index2)
          wd_name = paste0("waldStat_", index, "vs", index2)
        }
        else{
          p_name = paste0("pvalue_", index2, "vs", index)
          wd_name = paste0("waldStat_", index2, "vs", index)
        }
        p_name_new = paste0("pvalue_", test_lineage, "vs", lineage)
        wd_name_new = paste0("waldStat_", test_lineage, "vs", lineage)
        fc_name_new = paste0("log2FC_", test_lineage, "vs", lineage)
        p_list <- c(p_list, p_name)
        wd_list = c(wd_list, wd_name)
        p_names <- c(p_names, p_name_new)
        wd_names <- c(wd_names, wd_name_new)
        fc_names <- c(fc_names, fc_name_new)
      }
    }
    #pattern test result
    p_values_p = pattern[,p_list]
    wd_values = pattern[,wd_list]
    ref <- rownames(pattern)
    
    p_values_p <- p_values_p[ref, , drop = FALSE]
    wd_values  <- wd_values[ref, , drop = FALSE]
    
    stopifnot(
      identical(ref, rownames(p_values_p)),
      identical(ref, rownames(wd_values))
    )
    
    pattern_sel <- cbind(
      pattern[, c(1, 3), drop = FALSE],
      p_values_p,
      wd_values
    )
    colnames(pattern_sel) <- c("waldStat_combined", "pvalue_combined", p_names, wd_names)
    colnames(pattern_sel) <- paste0(colnames(pattern_sel), "_pattern")
    #diffend tes result
    p_values_d = diffend[,p_list]
    wd_values = diffend[,wd_list]
    ref <- rownames(diffend)
    
    p_values_d <- p_values_d[ref, , drop = FALSE]
    wd_values  <- wd_values[ref, , drop = FALSE]
    
    stopifnot(
      identical(ref, rownames(p_values_d)),
      identical(ref, rownames(wd_values))
    )
    
    diffend_sel <- cbind(
      diffend[, c(1, 3), drop = FALSE],
      p_values_d,
      wd_values
    )
    colnames(diffend_sel) <- c("waldStat_combined", "pvalue_combined", p_names, wd_names)
    colnames(diffend_sel) <- paste0(colnames(diffend_sel), "_diffend")
    pattern_sel <- pattern_sel[rownames(diffend_sel), ]
    res <- cbind(pattern_sel, diffend_sel)
    res_sel = as.matrix(res)
    if(nrow(res_sel) == 0){
      return(NULL)
    }
    gene_names = rownames(res_sel)
    FCs = calculate_dynamic_FC(cds, test_lineage, gene_names, lineages)
    FCs = t(FCs)
    FCs_sel = FCs
    colnames(FCs_sel) <- fc_names
    res_sel = res_sel[rownames(FCs_sel),]
    combined_pvalue_pattern <- apply(p_values_p, 1, get_meta_p)
    combined_pvalue_diffend <- apply(p_values_d, 1, get_meta_p)
    average_FC = apply(FCs_sel, 1, get_average_FC)
    #Reorder the rows
    ref <- rownames(res_sel)
    combined_pvalue_pattern <- combined_pvalue_pattern[ref]
    combined_pvalue_diffend <- combined_pvalue_diffend[ref]
    average_FC <- average_FC[ref]
    if (!identical(rownames(FCs_sel), names(average_FC))) {
      stop("Error: rownames(FCs_sel) and names(average_FC) do not match.")
    }
    final_res = cbind(res_sel, FCs_sel, combined_pvalue_pattern, combined_pvalue_diffend, average_FC)
    colnames(final_res) <- c(colnames(res_sel), colnames(FCs_sel), c("meta_p_pattern", "meta_p_diffend", "average_FC"))
    final_res = as.data.frame(final_res)
    final_res
  }
  else{
    pattern_sel = pattern[,c("waldStat", "pvalue")]
    pattern_sel = as.matrix(pattern_sel)
    diffend_sel = diffend[,c("waldStat", "pvalue")]
    diffend_sel = as.matrix(diffend_sel)
    gene_names = rownames(pattern_sel)
    FCs = calculate_dynamic_FC_single(cds, test_lineage, gene_names, lineages[lineages != test_lineage])
    FCs_sel = FCs
    pattern_sel = pattern_sel[names(FCs_sel),]
    diffend_sel = diffend_sel[names(FCs_sel),]
    final_res = as.data.frame(cbind(pattern_sel, diffend_sel, FCs_sel))
    colnames(final_res) <- c("waldStat_pattern", "meta_p_pattern", "waldStat_diffend", "meta_p_diffend", "average_FC")
    lineage_genes = as.data.frame(final_res)
    if(length(lineage_genes) > 0){
      lineage_genes
    }
    else{return(NULL)}
  }
}

calculate_dynamic_FC_single <- function(cds, test_lineage, genes, comp_lineage){
  genes = genes[genes %in% colnames(cds@expectation[[test_lineage]])]
  FCs = sapply(genes, calculate_dynamic_FC_single_gene, cds = cds, test_lineage = test_lineage, comp_lineage = comp_lineage)
  names(FCs) <- genes
  FCs
}

calculate_dynamic_FC_single_gene <- function(gene, cds, test_lineage, comp_lineage){
  exp1 = cds@expectation[[test_lineage]]
  exp2 = cds@expectation[[comp_lineage]]
  grid = seq(0, 1, length.out = nrow(exp1))
  auc_dataset1 <- trapz(grid, exp1[,gene])
  auc_dataset2 <- trapz(grid, exp2[,gene])
  auc_difference <- log2(auc_dataset1/auc_dataset2)
  auc_difference 
}

calculate_dynamic_FC <- function(cds, test_lineage, genes, lineages){
  genes = genes[genes %in% colnames(cds@expectation[[test_lineage]])]
  FC_matrix = sapply(genes, calculate_dynamic_FC_gene, cds = cds, test_lineage = test_lineage, lineages = lineages)
  rownames(FC_matrix) <- names(cds@lineages)[names(cds@lineages) != test_lineage]
  FC_matrix
}

calculate_dynamic_FC_gene <- function(gene, cds, test_lineage, lineages){
  exp1 = cds@expectation[[test_lineage]]
  grid = seq(0, 1, length.out = nrow(exp1))
  FCs = c()
  for(lin in lineages){
    if(lin != test_lineage){
      exp2 = cds@expectation[[lin]]
      auc_dataset1 <- trapz(grid, exp1[,gene])
      auc_dataset2 <- trapz(grid, exp2[,gene])
      auc_difference <- log2(auc_dataset1/auc_dataset2)
      FCs <- c(FCs, auc_difference)
    }
  }
  FCs
}

lineage_specific_genes_par <- function(cds, lineages = names(cds@lineages), parallel = F, BPPARAM = F, filter_by_dyn = TRUE){
  counts = matrix(,nrow = nrow(cds),ncol = 0)
  all_metacells = c()
  all_size_factor = c()
  lineage_list = list()
  pt_list = list()
  i = 1
  for(lineage in lineages){
    metacells = paste0(lineage, "_", c(1:nrow(cds@expression[[lineage]][['sum']])))
    d = cds@expression[[lineage]][['sum']]
    d = t(as.matrix(sapply(d[,7:ncol(d)], as.numeric)))
    d <- d[rownames(cds),]
    colnames(d) <- metacells
    counts <- cbind(counts, d)
    all_metacells <- c(all_metacells, metacells)
    size_factor = cds@expression[[lineage]][['sum']]$size_factor
    all_size_factor <- c(all_size_factor, size_factor)
    lineage_list[[i]] <- metacells
    #pt = cds@pseudotime[[lineage]]
    pt = cds@expression[[lineage]][["sum"]]$pseudotime
    names(pt) <- metacells
    pt_list[[i]] <- pt
    i <- i + 1
  }
  names(lineage_list) <- lineages
  names(pt_list) <- lineages
  cellWeights <- matrix(0, nrow = length(all_metacells), ncol = length(lineages))
  rownames(cellWeights) <- all_metacells
  colnames(cellWeights) <- lineages
  for (list_name in names(lineage_list)) {
    cellWeights[, list_name] <- as.numeric(all_metacells %in% lineage_list[[list_name]])
  }
  pseudotime <- matrix(0, nrow = length(all_metacells), ncol = length(lineages))
  rownames(pseudotime) <- all_metacells
  colnames(pseudotime) <- lineages
  for (list_name in names(pt_list)) {
    pseudotime[names(pt_list[[list_name]]), list_name] <- pt_list[[list_name]]
  }
  colnames(pseudotime) <- lineages
  rownames(pseudotime) <- all_metacells
  counts = counts[rownames(cds),]
  print(paste0("Testing ", nrow(counts), " genes"))
  gamlist = tradeSeq::fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, U = NULL, nknots = 6, offset = log(all_size_factor), parallel = parallel, BPPARAM = BPPARAM)
  pattern = tradeSeq::patternTest(models = gamlist, global = T, pairwise = T)
  diffend = tradeSeq::diffEndTest(models = gamlist, global = T, pairwise = T)
  out <- sapply(lineages, lineage_specific_genes_v2, cds = cds, pattern = pattern, diffend = diffend, lineages = lineages, simplify = FALSE)
  cds@lineage_genes <- out
  cds 
}

lineage_specific_genes <- function(test_lineage, cds, U = NULL, nknots = 6, parallel = F, BPPARAM = BPPARAM, p_cutoff = 0.05, FC_cutoff = 0, dyn_FC_cutoff = 0){
    library(monocle3)
    library(igraph)
    library(ggplot2)
    library(pbapply)
    library(devtools)
    library(dplyr)
    library(plotly)
    library(parallel)
    library(evobiR)
    library(shiny)
    library(colorspace)
    library(BiocParallel)
    library(tradeSeq)
    library(stringr)
    library(ggnewscale)
    library(pracma)
    library("metap")
    source_url("https://raw.githubusercontent.com/velmeshevlab/meta-tracker/dev/R/between_lineage_DE.R")
    source_url("https://raw.githubusercontent.com/velmeshevlab/meta-tracker/dev/R/compress_lineages.R")
    source_url("https://raw.githubusercontent.com/velmeshevlab/meta-tracker/dev/R/differential_expression.R")
    source_url("https://raw.githubusercontent.com/velmeshevlab/meta-tracker/dev/R/plotting.R")
    source_url("https://raw.githubusercontent.com/velmeshevlab/meta-tracker/dev/R/select_lineages.R")
    lineages = names(cds@lineages)
    lineage_genes = get_lineage_genes(cds, test_lineage, U = U, nknots = nknots, parallel = parallel, BPPARAM = BPPARAM, p_cutoff = p_cutoff, FC_cutoff = FC_cutoff, lineages = lineages, dyn_FC_cutoff = dyn_FC_cutoff)
    if(length(lineage_genes) > 0){
    lineage_genes = lineage_genes[with(lineage_genes, order(meta_p, -abs(average_FC))), ]
    lineage_genes
    }
    else{return(NULL)}
}

branch_specific_genes <- function(cds, test_lineages, name, U = NULL, nknots = 6, parallel = F, p_cutoff = 0.05, FC_cutoff = 0, dyn_FC_cutoff = 0){
  lineages = names(cds@lineages)
  genes = c()
  for(test_lineage in test_lineages){
    dynamic = cds@dynamic_genes[[test_lineage]]
    dynamic_genes = rownames(dynamic[dynamic$scaled_FC >= dyn_FC_cutoff,])
    genes = c(genes, dynamic_genes)
    }
  genes = unique(genes)
  for(lineage in lineages){
    genes = genes[genes %in% colnames(cds@expectation[[lineage]])]
    }
  print(paste0("Testing ", length(genes), " genes"))
  all_Ps <- matrix(0, nrow = length(genes), ncol = 0)
  all_FCs <- matrix(0, nrow = length(genes), ncol = 0)
  for(test_lineage in test_lineages){
    print(paste0("Testing lineage ", test_lineage))
    lineage_genes = get_lineage_genes(cds, test_lineage, genes = genes, U = U, nknots = nknots, parallel = parallel, p_cutoff = F, FC_cutoff = F, lineages = c(test_lineage, lineages[!(lineages %in% test_lineages)]))
    if(ncol(lineage_genes) == 2){
      Ps = as.matrix(lineage_genes[,1])
      rownames(Ps) <- rownames(lineage_genes)
      colnames(Ps) <- paste0("pvalue_",  test_lineage)
      FCs = as.matrix(lineage_genes[,2])
      rownames(FCs) <- rownames(lineage_genes)
      colnames(FCs) <- paste0("FC_",  test_lineage)
      }
    else{
      lineage_genes = lineage_genes[,1:(ncol(lineage_genes)-2)]
      Ps = lineage_genes[,1:(ncol(lineage_genes)/2)]
      FCs = lineage_genes[,((ncol(lineage_genes)/2)+1):ncol(lineage_genes)]
      }
    all_Ps = cbind(all_Ps, Ps)
    all_FCs = cbind(all_FCs, FCs)
    }
  #all_FCs_sel = all_FCs[rowSums(FCs >= FC_cutoff) == ncol(FCs) | rowSums(FCs <= -FC_cutoff) == ncol(FCs), ]
  all_FCs_sel = all_FCs[rowSums(FCs >= FC_cutoff) == ncol(FCs), ]
  all_Ps_sel = all_Ps[rownames(all_FCs_sel),]
  all_Ps_sel = all_Ps_sel[rowSums(all_Ps_sel < p_cutoff) == ncol(all_Ps_sel), ]
  all_FCs_sel = all_FCs_sel[rownames(all_Ps_sel),]
  meta_p <- apply(all_Ps_sel, 1, get_meta_p)
  average_FC = apply(all_FCs_sel, 1, get_average_FC)
  res = as.data.frame(cbind(all_Ps_sel, all_FCs_sel, meta_p, average_FC))
  res = res[with(res, order(meta_p, -abs(average_FC))), ]
  cds@lineage_genes[[name]] <- res
  cds
}

get_lineage_genes <- function(cds, test_lineage, genes = NULL, U = NULL, nknots = 6, parallel = F, BPPARAM = F, p_cutoff = F, FC_cutoff = 0, lineages = NULL, dyn_FC_cutoff = 0){
  counts = matrix(,nrow = ncol(cds@expression[[lineages[1]]])-3,ncol = 0)
  all_metacells = c()
  lineage_list = list()
  pt_list = list()
  i = 1
  for(lineage in lineages){
    metacells = paste0(lineage, "_", c(1:nrow(cds@expression[[lineage]])))
    d = cds@expression[[lineage]]
    d = t(as.matrix(sapply(d[,4:ncol(d)], as.numeric)))
    colnames(d) <- metacells
    counts <- cbind(counts, d)
    all_metacells <- c(all_metacells, metacells)
    lineage_list[[i]] <- metacells
    pt = cds@pseudotime[[lineage]][,1]
    names(pt) <- metacells
    pt_list[[i]] <- pt
    i <- i + 1
    }
  dynamic = cds@dynamic_genes[[test_lineage]]
  dynamic_genes = rownames(dynamic[dynamic$scaled_FC >= dyn_FC_cutoff,])
  if(length(genes) == 0){
    counts = counts[dynamic_genes,]
  }
  else{
    counts = counts[genes,]
  }
  print(paste0("Testing ", nrow(counts), " genes"))
  names(lineage_list) <- lineages
  names(pt_list) <- lineages
  cellWeights <- matrix(0, nrow = length(all_metacells), ncol = length(lineages))
  rownames(cellWeights) <- all_metacells
  colnames(cellWeights) <- lineages
  for (list_name in names(lineage_list)) {
    cellWeights[, list_name] <- as.numeric(all_metacells %in% lineage_list[[list_name]])
    }
  pseudotime <- matrix(0, nrow = length(all_metacells), ncol = length(lineages))
  rownames(pseudotime) <- all_metacells
  colnames(pseudotime) <- lineages
  for (list_name in names(pt_list)) {
    pseudotime[names(pt_list[[list_name]]), list_name] <- pt_list[[list_name]]
    }
  colnames(pseudotime) <- lineages
  rownames(pseudotime) <- all_metacells
  gamlist = tradeSeq::fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, U = U, nknots = nknots, parallel = parallel, BPPARAM = BPPARAM)
  res = tradeSeq::patternTest(models = gamlist, global = T, pairwise = T)
  if(length(lineages) > 2){
    index = which(test_lineage == lineages)
    p_list = c()
    p_names = c()
    for(linege in lineages){
      index2 = which(linege == lineages)
      if(index != index2){
        if(index<index2){
          p_name = paste0("pvalue_", index, "vs", index2)
          }
          else{
          p_name = paste0("pvalue_", index2, "vs", index)
          }
        p_name_new = paste0("pvalue_", test_lineage, "vs", linege)
        p_list <- c(p_list, p_name)
        p_names <- c(p_names, p_name_new)
        }
    }
    p_values = res[,p_list]
    colnames(p_values) <- p_names
    if(p_cutoff != F){
      p_values_sel = as.matrix(p_values[rowSums(p_values < p_cutoff) == ncol(p_values), ])
      }
    else{
      p_values_sel = as.matrix(p_values)
      }
    if(nrow(p_values_sel) == 0){
      return(NULL)
    }
    gene_names = rownames(p_values_sel)
    FCs = calculate_dynamic_FC(cds, test_lineage, gene_names, lineages)
    FCs = t(FCs)
    if(FC_cutoff != F){
        #FCs_sel = FCs[rowSums(FCs >= FC_cutoff) == ncol(FCs) | rowSums(FCs <= -FC_cutoff) == ncol(FCs), ]
        FCs_sel = FCs[rowSums(FCs >= FC_cutoff) == ncol(FCs), ]
      }
    else{
      FCs_sel = FCs
      }
    p_values_sel = p_values_sel[rownames(FCs_sel),]
    combined_pvalue <- apply(p_values_sel, 1, get_meta_p)
    average_FC = apply(FCs_sel, 1, get_average_FC)
    colnames_old = c(colnames(p_values_sel), colnames(FCs_sel))
    final_res = cbind(p_values_sel, FCs_sel, combined_pvalue, average_FC)
    colnames(final_res) <- c(colnames_old, c("meta_p", "average_FC"))
    final_res = as.data.frame(final_res)
    final_res
  }
  else{
    if(p_cutoff != F){
      res_sel = res[res$pvalue < p_cutoff,]
      }
    else{
      res_sel = res
      }
    p_values_sel = as.matrix(res_sel[,"pvalue"])
    rownames(p_values_sel) <- rownames(res_sel)
    colnames(p_values_sel) <- paste0("pvalue_", test_lineage, "vs", lineages[lineages != test_lineage])
    gene_names = rownames(p_values_sel)
    FCs = calculate_dynamic_FC_single(cds, test_lineage, gene_names, lineages[lineages != test_lineage])
    if(FC_cutoff != F){
        #FCs_sel = FCs[FCs >= FC_cutoff | FCs <= -FC_cutoff]
        FCs_sel = FCs[FCs >= FC_cutoff]
      }
    else
      {
      FCs_sel = FCs
      }
    p_values_sel = p_values_sel[names(FCs_sel),]
    final_res = as.data.frame(cbind(p_values_sel, FCs_sel))
    colnames(final_res) <- c("meta_p", "average_FC")
    final_res = as.data.frame(final_res)
    final_res
  }
  }

get_average_FC <- function(FC){
  linear_values <- 2^FC
  mean_linear <- mean(linear_values)
  mean_log2 <- log2(mean_linear)
  mean_log2
}

get_meta_p <- function(Ps) {
  Ps <- Ps[!is.na(Ps)]
  
  # Replace zeros with a tiny positive number
  Ps[Ps == 0] <- .Machine$double.xmin
  
  if (length(Ps) >= 2) {
    sumlog(Ps)$p
  } else {
    NA_real_  # Not enough p-values to combine
  }
}

between_lineage_DE <- function(counts, # A matrix with genes in rows and cells in columns, cells should be aligned in the order of separate lineages
                               pseudotime, # A matrix of pseudotime values, each row represents a cell and each column represents a lineage, the order of lineages should correspond to the order of cells in "counts"
                               cellWeights, # A matrix of cell weights defining the probability that a cell belongs to a particular lineage.
                               U = NULL, # design matrix for covariates
                               nknots = 6, # number of knots to fit the GAM
                               npoints = 2*nknots, # number of points to be compared between lineages
                               lineage_labels, # vector of lineage names, should correspond to the order of cells in "counts"
                               padjust_method = "BH" # for other methods, see ?p.adjust()
) {
  # Run tradeSeq test
  set.seed(10000)
  gamlist = tradeSeq::fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, U = U, nknots = nknots)
  res = tradeSeq::patternTest(models = gamlist, global = TRUE, pairwise = TRUE, nPoints = npoints)
  
  # Extract p values
  n_lineage = length(lineage_labels)
  
  pvalue = res[,grepl("pvalue_.*",colnames(res))]
  for (i in 1:n_lineage){
    colnames(pvalue)=gsub(as.character(i),lineage_labels[i],colnames(pvalue))
  }
  
  # Adjust p values
  padj = data.frame(apply(pvalue,2,function(x) p.adjust(x, method=padjust_method)))
  
  # Extract DE genes
  DE_list = list()
  
  for (i in 1:n_lineage){
    baseline_label = lineage_labels[i]
    baseline_padj = padj[,grepl(baseline_label,colnames(padj))]
    baseline_DE <- baseline_padj %>% 
      mutate(highest=apply(baseline_padj,1,max)) %>%
      filter(highest<0.05)
    DE_list[[i]] = baseline_DE
  }
  
  names(DE_list)=lineage_labels
  
  return(DE_list)
}
