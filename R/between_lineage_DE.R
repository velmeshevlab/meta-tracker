format_lineage_genes <- function(cds){
lineages <- names(cds@lineages)
for (lineage in lineages) {
  # 1. Generate the filtered genes
  lin_genes <- format_lineage_specific_genes(lineage, cds)
  
  # 2. Get existing data
  lineage_genes_orig <- cds_F@lineage_genes[[lineage]]
  
  # 3. Update the global object DIRECTLY
  cds@lineage_genes[[lineage]] <- list(
    "lineage_genes" = lineage_genes_orig, 
    "filtered" = lin_genes
  )
}
  cds
}

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
    lineage_genes = cds@lineage_genes[[lineage]][["lineage_genes"]][["pattern_test"]]
    #First filter out genes that express in fewer than 100 cells from lineages in both branches
    expressed_genes <- filter_by_expression_lineage(cds=cds, lineage = lineage, mode = "number", N = 100, ratio = 0.01)
    common_genes <- intersect(expressed_genes, rownames(lineage_genes))
    lineage_genes <- lineage_genes[common_genes, , drop = FALSE]
    #Second filter out genes that has p value >= 0.05 and moran's I statistics above 0.1
    if(dynamic_FC_cutoff != F){
      dynamic = cds@dynamic_genes[[lineage]]
      if(dynamic_test == "Moran"){
        dynamic_genes = rownames(dynamic[dynamic$I >= dynamic_FC_cutoff & dynamic$padj <= dynamic_p_cutoff, ])
      }
      else{
        dynamic_genes = rownames(dynamic[dynamic$scaled_FC >= dynamic_FC_cutoff, ])
      }
      common_genes <- intersect(dynamic_genes, rownames(lineage_genes))
      lineage_genes <- lineage_genes[common_genes, , drop = FALSE]
    }
    FC_names = c()
    for(lin in branches_2){
      FC_names = c(FC_names, paste0("log2FC_", lineage, "vs", lin, "_pattern"))
    }
    p_names = c()
    for(lin in branches_2){
      p_names = c(p_names, paste0("pvalue_", lineage, "vs", lin, "_pattern"))
    }
    lineage_genes_filtered = lineage_genes[,c(p_names, FC_names)]
    FCs = lineage_genes_filtered[,grepl("log2FC", colnames(lineage_genes_filtered))]
    p_values <- lineage_genes_filtered[,grepl("pvalue", colnames(lineage_genes_filtered))]
    if(is.vector(FCs)){
      #Third filter out genes with NA values in averageFC or p values
      colnames(lineage_genes_filtered) <- c("meta_p", "averageFC")
      lineage_genes_filtered <- lineage_genes_filtered[(!is.na(lineage_genes_filtered$meta_p) & !is.na(lineage_genes_filtered$averageFC)),]
      lineage_genes_filtered$p_adjusted <- p.adjust(lineage_genes_filtered$meta_p, method = p_adjust)
      lineage_spec_genes = lineage_genes_filtered[lineage_genes_filtered$averageFC >= FC_cutoff & lineage_genes_filtered$p_adjusted <= p_cutoff, ]
    }
    else{
      average_FC = apply(FCs, 1, get_average_FC)
      meta_p <- apply(p_values, 1, get_meta_p)
      ref <- rownames(lineage_genes_filtered)
      col <- colnames(lineage_genes_filtered)
      average_FC <- average_FC[ref]
      meta_p <- meta_p[ref]
      lineage_genes_filtered <- cbind(lineage_genes_filtered, average_FC, meta_p)
      colnames(lineage_genes_filtered) <- c(col, "average_FC", "meta_p")
      #Third filter out genes with NA values in averageFC or p values
      lineage_genes_filtered <- lineage_genes_filtered[(!is.na(lineage_genes_filtered$meta_p) & !is.na(lineage_genes_filtered$average_FC)),]
      lineage_genes_filtered$p_adjusted <- p.adjust(lineage_genes_filtered$meta_p, method = p_adjust)
      FCs <- FCs[rownames(lineage_genes_filtered),]
      lineage_spec_genes = lineage_genes_filtered[rowSums(FCs >= FC_cutoff) == ncol(FCs) & lineage_genes_filtered$p_adjusted <= p_cutoff, ]
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
      FC_names = c(FC_names, paste0("log2FC_", lineage, "vs", lin, "_pattern"))
    }
    p_names = c()
    for(lin in branches_2){
      p_names = c(p_names, paste0("pvalue_", lineage, "vs", lin, "_pattern"))
    }
    if(ncol(lineage_genes) <= 3){
      FCs = lineage_genes[branch_gene_names, "averageFC", drop = FALSE]
      p_values <- lineage_genes[branch_gene_names, "meta_p", drop = FALSE]
    }else{
      FCs = lineage_genes[branch_gene_names, FC_names]
      p_values <- lineage_genes[branch_gene_names, p_names]
    }
    meta_p_matrix = cbind(meta_p_matrix, p_values)
    median_FC_matrix = cbind(median_FC_matrix, FCs)
  }
  if(ncol(meta_p_matrix)==1){
    meta_p <- setNames(meta_p_matrix[,1], rownames(meta_p_matrix))
    median_FC <- setNames(median_FC_matrix[,1], rownames(median_FC_matrix))
  }else{
    meta_p = apply(meta_p_matrix, 1, get_meta_p)
    median_FC = apply(median_FC_matrix, 1, median)
  }
  if(p_adjust != FALSE){
    meta_p <- p.adjust(meta_p, method = p_adjust)
  }
  meta_p <- meta_p[names(median_FC)]
  final_out = cbind(meta_p, median_FC)
  rownames(final_out) <- names(median_FC)
  final_out = as.data.frame(final_out)
  final_out$branch <- rep(names(branch_point)[branch_number], nrow(final_out))
  colnames(final_out) <- c("meta_p", "median_FC", "branch")
  final_out = final_out[with(final_out, order(-abs(median_FC), meta_p)), ]
  final_out
}
                        
format_lineage_specific_genes <- function(lineage, cds, p_cutoff = 0.05, FC_pattern_cutoff = 0.2, FC_diffend_cutoff = 0.2, dynamic_I_cutoff = 0.1, dynamic_p_cutoff = 0.05, p_adjust = "BH", specificity = "high", dynamic_test = "Moran"){
  lineages = names(cds@lineages)
  pattern_genes = cds@lineage_genes[[lineage]]$pattern_test
  diffend_genes = cds@lineage_genes[[lineage]]$diffend_test
  #First filter out genes that express in fewer than 100 cells
  expressed_genes <- filter_by_expression_lineage(cds=cds, lineage=lineage, mode = "number", N = 100, ratio = 0.01)
  common_genes <- intersect(expressed_genes, rownames(pattern_genes))
  pattern_genes <- pattern_genes[common_genes, , drop = FALSE]
  common_genes <- intersect(expressed_genes, rownames(diffend_genes))
  diffend_genes <- diffend_genes[common_genes, , drop = FALSE]
  #Second filter out genes that has p value >= 0.05 and moran's I statistics above 0.1
  if(dynamic_I_cutoff != F){
    dynamic = cds@dynamic_genes[[lineage]]
    if(dynamic_test == "Moran"){
      dynamic_genes = rownames(dynamic[dynamic$I >= dynamic_I_cutoff & dynamic$padj <= dynamic_p_cutoff, ])
    }
    else{
      dynamic_genes = rownames(dynamic[dynamic$scaled_FC >= dynamic_I_cutoff, ])
    }
    common_genes <- intersect(dynamic_genes, rownames(pattern_genes))
    pattern_genes <- pattern_genes[common_genes, , drop = FALSE]
    common_genes <- intersect(dynamic_genes, rownames(diffend_genes))
    diffend_genes <- diffend_genes[common_genes, , drop = FALSE]
  }
  #Third filter out genes with NA values in averageFC 
  pattern_genes <- pattern_genes[(!is.na(pattern_genes$average_FC_pattern) & !is.na(pattern_genes$meta_p_pattern)),]
  diffend_genes <- diffend_genes[(!is.na(diffend_genes$average_FC_diffend) & !is.na(diffend_genes$meta_p_diffend)),]
  if(p_adjust != FALSE){
    pattern_genes$p_adjusted <- p.adjust(pattern_genes$meta_p_pattern, method = p_adjust)
    diffend_genes$p_adjusted <- p.adjust(diffend_genes$meta_p_diffend, method = p_adjust)
  }
  if(length(lineages) > 2){
    FCs_p = pattern_genes[, grepl("log2FC", colnames(pattern_genes))]
    median_FC_p = apply(FCs_p, 1, median)
    median_FC_p = median_FC_p[rownames(pattern_genes)]
    pattern_genes$median_FC <- median_FC_p
    FCs_d = diffend_genes[, grepl("log2FC", colnames(diffend_genes))]
    median_FC_d = apply(FCs_d, 1, median)
    median_FC_d = median_FC_d[rownames(diffend_genes)]
    diffend_genes$median_FC <- median_FC_d
    if(specificity == "high"){
      #lineage_genes_p = pattern_genes[(rowSums(FCs_p >= FC_cutoff) == ncol(FCs_p) | rowSums(FCs <= -FC_cutoff) == ncol(FCs)) & pattern_filtered$p_adjusted <= p_cutoff, ]
      lineage_genes_p = pattern_genes[rowSums(FCs_p >= FC_pattern_cutoff) == ncol(FCs_p) & pattern_genes$p_adjusted <= p_cutoff, ]
      lineage_genes_d = diffend_genes[rowSums(FCs_d >= FC_diffend_cutoff) == ncol(FCs_d) & diffend_genes$p_adjusted <= p_cutoff, ]
    }
    else{
      lineage_genes_p = pattern_genes[pattern_genes$p_adjusted <= p_cutoff & pattern_genes$median_FC >= FC_pattern_cutoff, ]
      lineage_genes_d = diffend_genes[diffend_genes$p_adjusted <= p_cutoff & diffend_genes$median_FC >= FC_diffend_cutoff, ]
    }
    gene = union(rownames(lineage_genes_p), rownames(lineage_genes_d))
    p_df <- pattern_genes[gene, c("waldStat_combined_pattern"), drop = FALSE]
    d_df <- diffend_genes[gene, c("waldStat_combined_diffend"), drop = FALSE]
    p_df <- p_df[rownames(d_df), , drop = FALSE]
    merged <- cbind(p_df, d_df)
    merged$transientScore <- 
      rank(-merged$waldStat_combined_pattern, ties.method = "min")^2 + rank(-merged$waldStat_combined_diffend, ties.method = "min")^2
    lineage_spec_genes <- merged[order(merged$transientScore), ]
    lineage_spec_genes$lineage <- rep(lineage, nrow(lineage_spec_genes))
    lineage_spec <- list("pattern_filtered" = lineage_genes_p, "diffend_filtered" = lineage_genes_d, "combined" = lineage_spec_genes)
    lineage_spec
  }
  else{
    pattern_genes['median_FC'] <- pattern_genes['average_FC_pattern']
    diffend_genes['median_FC'] <- diffend_genes['average_FC_diffend']
    if(specificity == "high"){
      lineage_genes_p = pattern_genes[(pattern_genes$median_FC >= FC_pattern_cutoff) & pattern_genes$p_adjusted <= p_cutoff, ]
      lineage_genes_d = diffend_genes[(diffend_genes$median_FC >= FC_diffend_cutoff) & diffend_genes$p_adjusted <= p_cutoff, ]
    }
    else{
      lineage_genes_p = pattern_genes[(pattern_genes$median_FC >= FC_pattern_cutoff) & pattern_genes$p_adjusted <= p_cutoff, ]
      lineage_genes_d = diffend_genes[(diffend_genes$median_FC >= FC_diffend_cutoff) & diffend_genes$p_adjusted <= p_cutoff, ]
      message("Same as high specificity test for 2 lineages")
    }
    gene = union(rownames(lineage_genes_p), rownames(lineage_genes_d))
    p_df <- pattern_genes[gene, c("waldStat_pattern"), drop = FALSE]
    d_df <- diffend_genes[gene, c("waldStat_diffend"), drop = FALSE]
    p_df <- p_df[rownames(d_df), , drop = FALSE]
    merged <- cbind(p_df, d_df)
    merged$transientScore <- 
      rank(-merged$waldStat_pattern, ties.method = "min")^2 + rank(-merged$waldStat_diffend, ties.method = "min")^2
    lineage_spec_genes <- merged[order(merged$transientScore), ]
    lineage_spec_genes$lineage <- rep(lineage, nrow(lineage_spec_genes))
    lineage_spec <- list("pattern_filtered" = lineage_genes_p, "diffend_filtered" = lineage_genes_d, "combined" = lineage_spec_genes)
    lineage_spec
  }
}

lineage_specific_genes_v2 <- function(test_lineage, cds, pattern, diffend, genes = NULL, lineages){
  print(paste0("Testing ", test_lineage))
  if(length(lineages) > 2){
    index = which(test_lineage == lineages)
    p_list = c()
    p_names = c()
    wd_list = c()
    fc_list = c()
    wd_names = c()
    fc_names = c()
    flip_fc <- c()
    for(lineage in lineages){
      index2 = which(lineage == lineages)
      if(index != index2){
        if(index<index2){
          p_name = paste0("pvalue_", index, "vs", index2)
          wd_name = paste0("waldStat_", index, "vs", index2)
          fc_name = paste0("logFC", index, "_", index2)
          flip_fc <- c(flip_fc, FALSE) 
        }
        else{
          p_name = paste0("pvalue_", index2, "vs", index)
          wd_name = paste0("waldStat_", index2, "vs", index)
          fc_name = paste0("logFC", index2, "_", index)
          flip_fc <- c(flip_fc, TRUE)
        }
        p_name_new = paste0("pvalue_", test_lineage, "vs", lineage)
        wd_name_new = paste0("waldStat_", test_lineage, "vs", lineage)
        fc_name_new = paste0("log2FC_", test_lineage, "vs", lineage)
        p_list <- c(p_list, p_name)
        wd_list = c(wd_list, wd_name)
        fc_list = c(fc_list, fc_name)
        p_names <- c(p_names, p_name_new)
        wd_names <- c(wd_names, wd_name_new)
        fc_names <- c(fc_names, fc_name_new)
      }
    }
    #pattern test result
    p_values_p = pattern[,p_list]
    wd_values = pattern[,wd_list]
    combined_pvalue_pattern <- apply(p_values_p, 1, get_meta_p)
    ref <- rownames(pattern)
    
    p_values_p <- p_values_p[ref, , drop = FALSE]
    wd_values  <- wd_values[ref, , drop = FALSE]
    combined_pvalue_pattern <- combined_pvalue_pattern[ref]
    
    stopifnot(
      identical(ref, rownames(p_values_p)),
      identical(ref, rownames(wd_values)),
      identical(ref, names(combined_pvalue_pattern))
    )
    
    pattern_sel <- cbind(
      pattern[, c(1, 3), drop = FALSE],
      p_values_p,
      wd_values,
      combined_pvalue_pattern
    )
    colnames(pattern_sel) <- c("waldStat_combined", "pvalue_combined", p_names, wd_names, "meta_p")
    colnames(pattern_sel) <- paste0(colnames(pattern_sel), "_pattern")
    #diffend tes result
    p_values_d = diffend[,p_list]
    wd_values = diffend[,wd_list]
    fc_values = diffend[,fc_list]
    fc_values[, flip_fc] <- -fc_values[, flip_fc]
    combined_pvalue_diffend <- apply(p_values_d, 1, get_meta_p)
    average_FC = apply(fc_values, 1, get_average_FC)
    
    ref <- rownames(diffend)
    p_values_d <- p_values_d[ref, , drop = FALSE]
    wd_values  <- wd_values[ref, , drop = FALSE]
    fc_values  <- fc_values[ref, , drop = FALSE]
    combined_pvalue_diffend <- combined_pvalue_diffend[ref]
    average_FC <- average_FC[ref]
   
     stopifnot(
      identical(ref, rownames(p_values_d)),
      identical(ref, rownames(wd_values)),
      identical(ref, rownames(fc_values)),
      identical(ref, names(combined_pvalue_diffend)),
      identical(ref, names(average_FC))
    )
    
    diffend_sel <- cbind(
      diffend[, c(1, 3), drop = FALSE],
      p_values_d,
      wd_values,
      combined_pvalue_diffend,
      fc_values,
      average_FC
    )
    colnames(diffend_sel) <- c("waldStat_combined", "pvalue_combined", p_names, wd_names, "meta_p", fc_names, "average_FC")
    colnames(diffend_sel) <- paste0(colnames(diffend_sel), "_diffend")
    pattern_sel <- pattern_sel[rownames(diffend_sel), ]
    #Calculate the FC for pattern_test
    gene_names = rownames(pattern_sel)
    FCs = calculate_dynamic_FC(cds, test_lineage, gene_names, lineages)
    FCs = t(FCs)
    FCs_sel = FCs
    colnames(FCs_sel) <- paste0(fc_names, "_pattern")
    pattern_sel = pattern_sel[rownames(FCs_sel),]
    average_FC = apply(FCs_sel, 1, get_average_FC)
    average_FC <- average_FC[rownames(FCs_sel)]
    if (!identical(rownames(FCs_sel), names(average_FC))) {
      stop("Error: rownames(FCs_sel) and names(average_FC) do not match.")
    }
    pattern_fin = cbind(pattern_sel, FCs_sel, average_FC)
    colnames(pattern_fin) <- c(colnames(pattern_sel), colnames(FCs_sel), "average_FC_pattern")
    pattern_fin = as.data.frame(pattern_fin)
    diffend_fin = as.data.frame(diffend_sel)
    final_res = list("pattern_test" = pattern_fin, "diffend_test" = diffend_fin)
    final_res
  }
  else{
    pattern_sel = pattern[,c("waldStat", "pvalue")]
    diffend_sel = diffend[,c("waldStat", "pvalue", "logFC1_2")]
    if (test_lineage == lineages[2]) {
      diffend_sel$logFC1_2 <- -diffend_sel$logFC1_2
    }
    gene_names = rownames(pattern_sel)
    FCs = calculate_dynamic_FC_single(cds, test_lineage, gene_names, lineages[lineages != test_lineage])
    FCs_sel = FCs
    pattern_sel = pattern_sel[names(FCs_sel),]
    diffend_sel = diffend_sel[names(FCs_sel),]
    pattern_fin = as.data.frame(cbind(pattern_sel, FCs_sel))
    diffend_fin = as.data.frame(diffend_sel)
    colnames(pattern_fin) <- c("waldStat_pattern", "meta_p_pattern", "average_FC_pattern")
    colnames(diffend_fin) <- c("waldStat_diffend", "meta_p_diffend", "average_FC_diffend")
    lineage_genes = list("pattern_test" = pattern_fin, "diffend_test" = diffend_fin)
    lineage_genes
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

lineage_specific_genes_par <- function(cds, lineages = names(cds@lineages)){
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
    pt <-  (pt - min(pt))/(max(pt) - min(pt))
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
  gamlist = tradeSeq::fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, U = NULL, nknots = 6, offset = log(all_size_factor))
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
  Ps[Ps == 0] <- .Machine$double.xmin
  if (length(Ps) == 0) return(NA_real_)
  
  stat <- -2 * sum(log(Ps))
  df <- 2 * length(Ps)
  pchisq(stat, df=df, lower.tail=FALSE)
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

                          
prepare_pt <- function(cds, lineages){
  all_pt = c()
  i = 1
  for(lineage in lineages){
    metacells = paste0(lineage, "_", c(1:nrow(cds@expression[[lineage]][['sum']])))
    #pt = cds@pseudotime[[lineage]]
    pt = cds@expression[[lineage]][["sum"]]$pseudotime
    pt <-  (pt - min(pt))/(max(pt) - min(pt))
    names(pt) <- metacells
    all_pt <- c(all_pt, pt)
    i <- i + 1
  }
  return(all_pt)
}

quasi_refit <- function(cds, cores = 1, N = 1000, nknot = 8) {
  lineage_names <- names(cds@lineages)
  #prepare interior knots and boundary knots
  pseudotime <- prepare_pt(cds, lineage_names)
  interior_knots <- quantile(pseudotime, probs = seq(0, 1, length.out = nknot)[-c(1, nknot)])
  for (lineage in lineage_names) {
    message(paste("Processing lineage:", lineage))
    meta_sum_ordered <- cds@expression[[lineage]][["sum"]]
    mat <- meta_sum_ordered[,7:(ncol(meta_sum_ordered))]
    size_factor = meta_sum_ordered$size_factor
    d <-  (meta_sum_ordered$pseudotime - min(meta_sum_ordered$pseudotime))/(max(meta_sum_ordered$pseudotime)-min(meta_sum_ordered$pseudotime))
    predict_pt <- seq(0, 1, length.out = N)
    mat_m <- as.matrix(mat)
    genes <- colnames(mat_m)
    model <- as.formula(
      substitute(
        expression ~ splines::ns(pseudotime, knots = k, Boundary.knots = c(0, 1)) + offset(log(size_factor)),
        list(k = interior_knots)
      )
    )
    print(paste("Fitting curves for", lineage, "using scaled pseudotime"))
    fit_list_2 <- pbapply::pblapply(setNames(seq_along(genes), genes),function(i) {
      fit.m3_3_wald(exp.sel = mat_m[, i], pt = d, size_factor = size_factor, predict_pt = predict_pt, model = model, N = N)
    },
    cl = cores
    )
    cds@expectation[[lineage]] <- fit_list_2
  }
  return(cds)
}

calculate_separate_wald <- function(fit_A, fit_B, l2fc= 0.4, eigenThresh = 1e-8) {
  # 1. Extract betas and vcovs
  betaA <- fit_A$beta
  betaB <- fit_B$beta
  VA <- fit_A$vcov
  VB <- fit_B$vcov
  L <- fit_A$X 
  
  # 2. Calculate the difference (Contrast)
  estFC <- L %*% betaA - L %*% betaB
  logFCCutoff <- log(2^l2fc) # log2 to log scale
  est <- sign(estFC)*pmax(0, abs(estFC) - logFCCutoff) # zero or remainder
  
  # 3. Calculate Combined Sigma (Independence assumption)
  # Sigma = L*VA*L' + L*VB*L'
  # Optimized as: L %*% (VA + VB) %*% t(L)
  sigma <- L %*% (VA + VB) %*% t(L)
  
  # 4. Eigen-decomposition
  eSigma <- eigen(sigma, symmetric = TRUE)
  
  # 5. Determine rank using the PatternTest threshold
  r <- try(sum(eSigma$values / eSigma$values[1] > eigenThresh), silent = TRUE)
  
  if (inherits(r, "try-error") || is.na(r) || r < 1) {
    return(c(waldStat = NA, df = NA, p_val = NA, log2FC = NA, predictA = NA, predictB = NA))
  }
  
  if (r == 1) {
    # If rank is 1, halfCovInv is just a column vector
    # sigma = v * lambda * v' => sigma^-1 = v * (1/lambda) * v'
    # halfCovInv = v * 1/sqrt(lambda)
    halfCovInv <- eSigma$vectors[, 1, drop = FALSE] * (1 / sqrt(eSigma$values[1]))
  } else {
    # Standard matrix logic for r > 1
    inv_sqrt_vals <- 1 / sqrt(eSigma$values[1:r])
    halfCovInv <- eSigma$vectors[, 1:r, drop = FALSE] %*% diag(inv_sqrt_vals, nrow = r)
  }
  # 6. Project the difference into the eigen-space
  # t(halfCovInv) is (r x N), estFC is (N x 1) -> halfStat is (r x 1)
  halfStat <- t(halfCovInv) %*% est
  stat <- sum(halfStat^2)
  p_val   <- stats::pchisq(as.numeric(stat), df = r, lower.tail = FALSE)
  # 7. Calculate the AUC between fit_A and fit_B
  exp1 = fit_A$prediction
  exp2 = fit_B$prediction
  grid = seq(0, 1, length.out = length(exp1))
  auc_dataset1 <- trapz(grid, exp1)
  auc_dataset2 <- trapz(grid, exp2)
  auc_difference <- log2(auc_dataset1/auc_dataset2)
  pred_val_A <- auc_dataset1
  pred_val_B <- auc_dataset2
  return(c(waldStat = as.numeric(stat), df = r, p_val = p_val, log2FC = auc_difference, predictA = pred_val_A, predictB = pred_val_B))
}

run_pairwise_comparison <- function(list_A, list_B) {
  # 1. Identify common genes
  common_genes <- intersect(names(list_A), names(list_B))
  if (length(common_genes) == 0) stop("No common genes found.")
  
  list_A <- list_A[common_genes]
  list_B <- list_B[common_genes]
  message(paste("Processing", length(common_genes), "common genes..."))
  
  # 2. Process genes in parallel
  results_list <- pbapply::pblapply(common_genes, function(gene) {
    fit_A <- list_A[[gene]]
    fit_B <- list_B[[gene]]
    
    # Template for failure (ensuring column count consistency)
    empty_res_p <- c(waldStat = NA, df = NA, p_val = NA, log2FC = NA, predictA = NA, predictB = NA)
    empty_res_d <- c(waldStat = NA, df = NA, p_val = NA, log2FC = NA, predictA = NA, predictB = NA)
    
    # --- Convergence and Null Checks ---
    # Check if beta or vcov exists
    total_failure <- is.null(fit_A$beta)  || is.null(fit_B$beta) || 
      is.null(fit_A$vcov)  || is.null(fit_B$vcov) ||
      isFALSE(fit_A$converged) || isFALSE(fit_B$converged)
    
    if (total_failure) {
      return(list(pattern = empty_res_p, diffend = empty_res_d))
    }
    
    # --- Statistical Tests ---
    res_1 <- tryCatch({ 
      calculate_separate_wald(fit_A, fit_B) 
    }, error = function(e) empty_res_p)
    
    res_2 <- tryCatch({ 
      calculate_separate_diff(fit_A, fit_B) 
    }, error = function(e) empty_res_d)
    
    return(list(pattern = res_1, diffend = res_2))
  })
  
  # 3. Extract and format the Pattern Test Data Frame
  pattern_df <- as.data.frame(do.call(rbind, lapply(results_list, `[[`, "pattern")))
  rownames(pattern_df) <- common_genes
  colnames(pattern_df) <- c("waldStat", "df", "pvalue", "log2FC", "predictA", "predictB")
  
  # 4. Extract and format the DiffEnd Test Data Frame
  diffend_df <- as.data.frame(do.call(rbind, lapply(results_list, `[[`, "diffend")))
  rownames(diffend_df) <- common_genes
  colnames(diffend_df) <- c("waldStat", "df", "pvalue", "log2FC", "predictA", "predictB")
  
  # 5. Return as a list containing the two data frames
  return(list(
    pattern = pattern_df,
    diffend = diffend_df
  ))
}


run_all_pairwise <- function(cds, lineages) {
  # 1. Validation
  lineage_names <- lineages
  if (length(lineage_names) < 2) {
    stop("You need at least 2 lineages to perform pairwise comparisons.")
  }
  
  # 2. Generate all unique pairs
  pairs <- combn(lineage_names, 2, simplify = FALSE)
  
  # 3. Storage for results
  # We will store results in a nested list: results$pattern and results$diffend
  final_pattern_list <- list()
  final_diffend_list <- list()
  
  for (pair in pairs) {
    l1 <- pair[1]
    l2 <- pair[2]
    comp_name <- paste(l1, "vs", l2, sep = "_")
    
    message(sprintf("\n>>> Comparing: %s vs %s", l1, l2))
    
    # 4. Run your existing comparison logic
    # Assuming run_pairwise_comparison returns list(pattern = df, diffend = df)
    res <- run_pairwise_comparison(cds@expectation[[l1]], cds@expectation[[l2]])
    
    # 5. Store the dataframes
    final_pattern_list[[comp_name]] <- res$pattern
    final_diffend_list[[comp_name]] <- res$diffend
  }
  
  # 6. Return as a structured list
  return(list(
    pattern = final_pattern_list,
    diffend = final_diffend_list
  ))
}


calculate_global_wald <- function(fit_list, l2fc = 0.4, eigenThresh = 1e-2) {
  
  # 1. Filter for valid fits (Not NULL, have beta/vcov, and converged)
  is_valid <- sapply(fit_list, function(x) {
    !is.null(x) && 
      !is.null(x$beta) && 
      !is.null(x$vcov) && 
      !isFALSE(x$converged)
  })
  
  valid_fits <- fit_list[is_valid]
  k_valid <- length(valid_fits)
  
  # 2. Need at least two valid lineages to perform a comparison
  if (k_valid < 2) {
    return(c(waldStat = NA, df = NA, pvalue = NA))
  }
  
  # 3. Setup using only valid lineages
  # Assuming X is consistent across all lineages, take it from the first valid one
  L_single <- as.matrix(valid_fits[[1]]$X) 
  n_points <- nrow(L_single)
  n_params <- ncol(L_single)
  
  # 3. Stack Betas and create Block-Diagonal Sigma
  beta_stack <- do.call(rbind, lapply(valid_fits, function(x) as.matrix(x$beta)))
  V_block <- as.matrix(Matrix::bdiag(lapply(valid_fits, `[[`, "vcov")))
  
  # 4. Create the Omnibus Contrast Matrix (L_omni) for k_valid lineages
  # This compares Valid Lineage 1 against Lineages 2...k_valid
  L_omni <- matrix(0, nrow = (k_valid - 1) * n_points, ncol = k_valid * n_params)
  
  for (i in 2:k_valid) {
    row_idx <- ((i - 2) * n_points + 1):((i - 1) * n_points)
    
    # Reference (The first valid lineage)
    L_omni[row_idx, 1:n_params] <- L_single
    
    # Comparison (The i-th valid lineage)
    col_idx <- ((i - 1) * n_params + 1):(i * n_params)
    L_omni[row_idx, col_idx] <- -L_single
  }
  
  # 5. Calculate Global Difference and Sigma
  estFC <- L_omni %*% beta_stack
  logFCCutoff <- log(2^l2fc) 
  est <- sign(estFC) * pmax(0, abs(estFC) - logFCCutoff) 
  
  sigma <- L_omni %*% V_block %*% t(L_omni)
  
  # 6. Eigen-decomposition for Rank-Deficient Inverse
  eSigma <- eigen(sigma, symmetric = TRUE)
  
  # 7. Avoid division by zero if max eigenvalue is non-positive
  max_ev <- eSigma$values[1]
  if (is.na(max_ev) || max_ev <= 0) {
    return(c(waldStat = NA, df = NA, pvalue = NA))
  }
  
  r <- try(sum(eSigma$values / max_ev > eigenThresh), silent = TRUE)
  
  if (inherits(r, "try-error") || is.na(r) || r < 1) {
    return(c(waldStat = NA, df = NA, pvalue = NA))
  }
  
  # 8. Calculate the Wald Statistic
  inv_sqrt_vals <- 1 / sqrt(eSigma$values[1:r])
  if (r == 1) {
    halfCovInv <- eSigma$vectors[, 1, drop = FALSE] * inv_sqrt_vals
  } else {
    halfCovInv <- eSigma$vectors[, 1:r, drop = FALSE] %*% diag(inv_sqrt_vals, nrow = r)
  }
  
  halfStat <- t(halfCovInv) %*% est
  stat <- sum(halfStat^2)
  p_val <- stats::pchisq(as.numeric(stat), df = r, lower.tail = FALSE)
  
  return(c(waldStat = as.numeric(stat), df = r, pvalue = p_val))
}

calculate_separate_diff <- function(fit_A, fit_B, l2fc = 0.1, eigenThresh = 1e-2) {
  # 1. Extract parameters
  betaA <- fit_A$beta
  betaB <- fit_B$beta
  VA <- fit_A$vcov
  VB <- fit_B$vcov
  
  # 2. L is the contrast
  L <- matrix(as.numeric(fit_A$X_end), ncol = 1)
  
  # 3. Combined Variance-Covariance Matrix
  Sigma <- VA + VB
  
  # 4. Project Difference and Variance into Contrast Space
  # estFC is the actual log-natural difference
  beta_diff <- betaA - betaB
  estFC <- t(L) %*% beta_diff
  
  # 5. sigma is the variance of that difference
  sigma <- t(L) %*% Sigma %*% L
  
  # 6. Apply the LFC Threshold (The tradeSeq "FC" logic)
  logFCCutoff <- log(2^l2fc)
  est <- sign(estFC) * pmax(0, abs(estFC) - logFCCutoff)
  
  # 7. Eigendecomposition Method (tradeSeq Stability Logic)
  eSigma <- eigen(sigma, symmetric = TRUE)
  
  # 8. Determine rank based on eigenvalue threshold
  # We use the ratio of eigenvalues to the maximum eigenvalue
  r <- sum(eSigma$values / eSigma$values[1] > eigenThresh)
  
  # 9. Handle cases where the matrix is numerically zero or invalid
  if (r == 0 || is.na(r)) {
    return(c(waldStat = NA, df = NA, p_val = NA, log2FC = NA, predictA = NA, predictB = NA))
  }
  
  # 10. Compute the pseudo-inverse component (halfCovInv)
  # This follows: halfCovInv = Vectors * diag(1/sqrt(Values))
  halfCovInv <- eSigma$vectors[, seq_len(r), drop = FALSE] %*% 
    diag(1 / sqrt(eSigma$values[seq_len(r)]), nrow = r)
  
  # 11. Final Wald Statistic
  # stat = t(est) %*% Sigma_Inv %*% est
  halfStat <- t(est) %*% halfCovInv
  wald_stat <- as.numeric(crossprod(t(halfStat)))
  
  # 12. Degrees of Freedom and P-value
  df_wald <- r
  p_val <- stats::pchisq(wald_stat, df = df_wald, lower.tail = FALSE)
  
  # 13. Convert natural log FC to Log2FC for reporting
  log2FC_val <- as.numeric(estFC / log(2))
  
  pred_val_A <- fit_A$prediction[length(fit_A$prediction)]
  pred_val_B <- fit_B$prediction[length(fit_B$prediction)]
  
  return(c(waldStat = wald_stat, df = df_wald, p_val = p_val, log2FC = log2FC_val, predictA = pred_val_A, predictB = pred_val_B))
}

calculate_global_diffend <- function(list_of_fits, l2fc = 0.1, eigenThresh = 1e-2) {
  
  # 1. Filter for valid fits (Not NULL, have beta/vcov, and converged)
  is_valid <- sapply(list_of_fits, function(x) {
    !is.null(x) && 
      !is.null(x$beta) && 
      !is.null(x$vcov) && 
      !isFALSE(x$converged)
  })
  
  valid_fits <- list_of_fits[is_valid]
  n_valid <- length(valid_fits)
  
  # 2. Need at least two lineages to compare differences
  if (n_valid < 2) {
    return(c(waldStat = NA, df = NA, pvalue = NA))
  }
  
  # 3. Setup Dimensions using valid fits only
  K <- length(valid_fits[[1]]$beta)
  beta_all  <- do.call(c, lapply(valid_fits, `[[`, "beta"))
  Sigma_all <- as.matrix(Matrix::bdiag(lapply(valid_fits, `[[`, "vcov")))
  
  # 3. Build the Global Contrast Matrix (L)
  # We compare valid lineages 2, 3, ... N against valid lineage 1
  L_cols <- list()
  X_end  <- as.numeric(valid_fits[[1]]$X_end) 
  
  for (i in 2:n_valid) {
    l_vec <- rep(0, length(beta_all))
    l_vec[1:K] <- X_end                # Reference (First valid lineage)
    
    start_idx <- ((i - 1) * K) + 1
    end_idx   <- i * K
    l_vec[start_idx:end_idx] <- -X_end # Comparison (Lineage i)
    
    L_cols[[i-1]] <- l_vec
  }
  L_global <- do.call(cbind, L_cols)
  
  # 4. Project Differences and Variance
  estFC_vec <- t(L_global) %*% beta_all
  sigma_mat <- t(L_global) %*% Sigma_all %*% L_global
  
  # 5. Apply the LFC Threshold
  logFCCutoff <- log(2^l2fc)
  est_shrunk <- sign(estFC_vec) * pmax(0, abs(estFC_vec) - logFCCutoff)
  
  # 6. Eigendecomposition Method
  eSigma <- eigen(sigma_mat, symmetric = TRUE)
  max_ev <- eSigma$values[1]
  
  if (is.na(max_ev) || max_ev <= 0) {
    return(c(waldStat = NA, df = NA, pvalue = NA))
  }
  
  r <- sum(eSigma$values / max_ev > eigenThresh)
  
  if (is.na(r) || r < 1) {
    return(c(waldStat = NA, df = NA, pvalue = NA))
  }
  
  # 7. Compute Pseudo-inverse component
  halfCovInv <- eSigma$vectors[, seq_len(r), drop = FALSE] %*% 
    diag(1 / sqrt(eSigma$values[seq_len(r)]), nrow = r)
  
  # 8. Final Wald Statistic
  halfStat <- t(est_shrunk) %*% halfCovInv
  wald_stat <- as.numeric(crossprod(t(halfStat)))
  
  # 9. Result
  return(c(waldStat = wald_stat, df = r, pvalue = stats::pchisq(wald_stat, df = r, lower.tail = FALSE)))
}
                                      
                                      
run_global_comparison <- function(cds) {
  lineage_data_list <- cds@expectation
  common_genes <- Reduce(intersect, lapply(lineage_data_list, names))
  
  results_list <- pbapply::pblapply(common_genes, function(gene) {
    # Extract the fit for this gene from every lineage
    fits_for_gene <- lapply(lineage_data_list, function(lin) lin[[gene]])
    
    # Run the global Pattern test (using your calculate_global_wald)
    res_pattern <- tryCatch({
      calculate_global_wald(fits_for_gene)
    }, error = function(e) {
      return(c(waldStat = NA, df = NA, pvalue = NA))
    })
    
    # Run the global DiffEnd test (using your calculate_global_diffend)
    res_diffend <- tryCatch({
      calculate_global_diffend(fits_for_gene)
    }, error = function(e) {
      return(c(waldStat = NA, df = NA, pvalue = NA))
    })
    
    return(list(pattern = res_pattern, diffend = res_diffend))
  })
  
  # 1. Format Pattern results
  # We use lapply to extract the "pattern" element from each list item
  pattern_df <- as.data.frame(do.call(rbind, lapply(results_list, `[[`, "pattern")))
  rownames(pattern_df) <- common_genes
  
  # 2. Format DiffEnd results
  # We use lapply to extract the "diffend" element from each list item
  diffend_df <- as.data.frame(do.call(rbind, lapply(results_list, `[[`, "diffend")))
  rownames(diffend_df) <- common_genes
  
  # 3. Return as a named list
  return(list(
    pattern = pattern_df,
    diffend = diffend_df
  ))
}

extract_lineage_df <- function(res, lineage) {
  nm <- names(res)
  keep <- sapply(nm, function(x) {
    parts <- strsplit(x, "_vs_")[[1]]
    lineage %in% parts
  })
  sub_res <- res[keep]
  if (length(sub_res) == 0) return(NULL)
  
  out_list <- list()
  for (i in seq_along(sub_res)) {
    name <- names(sub_res)[i]
    df   <- sub_res[[i]]
    
    parts <- strsplit(name, "_vs_")[[1]]
    a <- parts[1] # The original 'Reference' (predictA)
    b <- parts[2] # The original 'Comparison' (predictB)
    
    # If the lineage we are interested in was 'b', it was the 'Comparison'.
    # To make it the 'Reference', we must flip FC and swap predictA/B values.
    if (lineage == b) {
      # 1. Flip log2FC
      if ("log2FC" %in% colnames(df)) {
        df$log2FC <- -df$log2FC
      }
      
      # 2. Swap PredictA and PredictB values
      if ("predictA" %in% colnames(df) && "predictB" %in% colnames(df)) {
        temp_A <- df$predictA
        df$predictA <- df$predictB
        df$predictB <- temp_A
      }
      
      other <- a
    } else {
      other <- b
    }
    
    new_name <- paste0(lineage, "_vs_", other)
    # Rename columns to include the comparison name for clarity
    colnames(df) <- paste0(colnames(df), "_", new_name)
    out_list[[new_name]] <- df
  }
  
  names(out_list) <- NULL
  # Combine all pairwise comparisons for this lineage into one wide data frame
  out <- do.call(cbind, out_list)
  return(out)
}

reorder_metrics_grouped <- function(df, type) {
  cn <- colnames(df)
  
  # 1. Identify global columns
  global_cols <- c("waldStat", "pvalue", "df")
  global_present <- global_cols[global_cols %in% cn]
  
  # Identify pairwise columns (everything else)
  pairwise <- setdiff(cn, global_present)
  
  # 2. Extract base comparison names
  # Logic: Remove the metric prefix (waldStat_|pvalue_|df_|log2FC_) 
  # AND remove the suffix (_diffend)
  # Result should be "ExN12_vs_ExN2", etc.
  base_names <- sub("^(waldStat_|pvalue_|df_|log2FC_|predictA_|predictB_)", "", pairwise)
  unique_comps <- unique(base_names)
  
  # 3. NATURAL SORT (ExN2 before ExN10)
  if (requireNamespace("gtools", quietly = TRUE)) {
    comp_order <- gtools::mixedsort(unique_comps)
  } else {
    comp_order <- sort(unique_comps)
  }
  
  # 4. Reconstruct columns in grouped order
  # Structure: metric_COMP_diffend
  wald_cols   <- paste0("waldStat_", comp_order)
  pval_cols   <- paste0("pvalue_",   comp_order)
  df_cols     <- paste0("df_",       comp_order)
  l2fc_cols   <- paste0("log2FC_",   comp_order)
  predict_A_cols   <- paste0("predictA_",   comp_order)
  predict_B_cols   <- paste0("predictB_",   comp_order)
  
  # Filter for columns that actually exist
  ordered_pairwise <- c(wald_cols, pval_cols, df_cols, l2fc_cols, predict_A_cols, predict_B_cols)
  ordered_pairwise <- ordered_pairwise[ordered_pairwise %in% cn]
  
  # 5. Final Order
  final_order <- c(global_present, ordered_pairwise)
  df <- df[, final_order, drop = FALSE]
  cn_fin <- colnames(df)
  # Check if name is exactly waldStat, pvalue, or df
  is_global <- cn_fin %in% c("waldStat", "pvalue", "df")
  if(type == "pattern"){
    colnames(df) <- ifelse(is_global, 
                           paste0(cn_fin, "_combined_pattern"), 
                           paste0(cn_fin, "_pattern"))
    colnames(df) <- gsub("_vs_", "vs", colnames(df))
  }
  if(type == "diffend"){
    colnames(df) <- ifelse(is_global, 
                           paste0(cn_fin, "_combined_diffend"), 
                           paste0(cn_fin, "_diffend"))
    
  }
  colnames(df) <- gsub("_vs_", "vs", colnames(df))
  return(df)
}

quasi_test <- function(cds){
  lineages <- names(cds@lineages)
  
  # Both 'res' and 'waldResOmnibus' are now lists with names 'pattern' and 'diffend'
  res <- run_all_pairwise(cds, lineages)
  waldResOmnibus <- run_global_comparison(cds)
  
  for (lin in lineages) {
    # 1. Extract pairwise dataframes for this lineage
    # Assuming 'extract_lineage_df' handles one test type at a time
    out_pattern <- extract_lineage_df(res$pattern, lin)
    out_diffend <- extract_lineage_df(res$diffend, lin)
    
    # 2. Match rows with the Global (Omnibus) results for each test
    # This ensures genes align across Global and Pairwise results
    global_p <- waldResOmnibus$pattern[rownames(out_pattern), ]
    global_d <- waldResOmnibus$diffend[rownames(out_diffend), ]
    
    df_pattern <- cbind(global_p, out_pattern)
    df_diffend <- cbind(global_d, out_diffend)
    
    cds@lineage_genes[[lin]] <- list() 
    cds@lineage_genes[[lin]]$lineage_genes <- list()
    
    # 3. Reorder metrics and store as a list of 2 dataframes
    # Reusing 'reorder_metrics_grouped' which now handles the 'global_' prefix
    pattern_test = reorder_metrics_grouped(df_pattern, type = "pattern")
    fc_cols <- grep("^log2FC_.*_pattern$", colnames(pattern_test), value = TRUE)
    if (length(fc_cols) == 1) {
      # Direct assignment is much faster than apply
      pattern_test$average_FC_pattern <- pattern_test[[fc_cols]]
    } else {
      # Use apply for multiple columns
      pattern_test$average_FC_pattern <- apply(pattern_test[, fc_cols], 1, get_average_FC)
    }
    diffend_test = reorder_metrics_grouped(df_diffend, type = "diffend")
    fc_cols <- grep("^log2FC_.*_diffend$", colnames(diffend_test), value = TRUE)
    if (length(fc_cols) == 1) {
      # Direct assignment is much faster than apply
      diffend_test$average_FC_diffend <- diffend_test[[fc_cols]]
    } else {
      # Use apply for multiple columns
      diffend_test$average_FC_diffend <- apply(diffend_test[, fc_cols], 1, get_average_FC)
    }
    
    cds@lineage_genes[[lin]]$lineage_genes$quasipoisson <- list(
      pattern_test = pattern_test,
      diffend_test = diffend_test
    )
  }
  return(cds)
}

fit.m3_3_wald <- function(exp.sel, pt, size_factor, predict_pt, model, N) {
  
  if (!requireNamespace("speedglm", quietly = TRUE)) {
    stop("speedglm not installed")
  }
  
  exp_data.sel <- data.frame(
    pseudotime  = as.numeric(pt),
    size_factor = as.numeric(size_factor),
    expression  = as.numeric(exp.sel)
  )
  
  # Try to fit the model and extract wald components
  result <- tryCatch({
    
    # 1. Fit the model
    fit_model <- speedglm::speedglm(
      model,
      data = exp_data.sel,
      family = quasipoisson(),
      acc = 1e-3,
      model = TRUE, # We need this temporarily to extract terms
      y = FALSE
    )
    
    # 2. Get Predicted Values (y-hat) for the 100 points
    grid_data <- data.frame(
      pseudotime  = as.numeric(predict_pt),
      size_factor = 1
    )
    pred <- predict(fit_model, newdata = grid_data, type = "response")
    
    # Determine number of points for Wald stats: 2 * nknots (number of coefficients)
    n_coefs <- length(coef(fit_model))
    # nPoints_wald <- 2 * (n_coefs - 1) # If nknots means just the spline bases
    nPoints_wald <- 2 * n_coefs # If nknots means total parameters as per the user logic
    
    # 2. Wald Grid (100 points for the statistics)
    # This is what you use for the p-value math
    grid_wald <- data.frame(
      pseudotime = seq(0, 1, length.out = nPoints_wald), 
      size_factor = 1
    )
    
    end_pt_value <- 1
    df_end <- data.frame(pseudotime = end_pt_value, size_factor = 1)
    
    # 3. Extract Basis Matrix (X) for the 100 points
    # This ensures the Wald test uses the same basis as the fit
    model_terms <- terms(fit_model)
    clean_terms <- delete.response(model_terms)
    X_matrix <- model.matrix(clean_terms, data = grid_wald)
    X_end <- model.matrix(clean_terms, data = df_end)
    
    # 4. Return all pieces needed for Wald
    list(
      prediction = as.numeric(pred),
      beta       = coef(fit_model),
      vcov       = vcov(fit_model),
      X          = X_matrix,
      X_end      = X_end,
      converged  = fit_model$convergence
    )
    
  }, error = function(e) {
    # Return NAs/NULLs if the gene fails to fit
    list(
      prediction = rep(NA_real_, N),
      beta       = NULL,
      vcov       = NULL,
      X          = NULL,
      X_end      = NULL,
      converged  = FALSE
    )
  })
  
  return(result)
}

                          
