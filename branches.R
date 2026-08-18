# Branch-point detection and lineage-specific gene formatting.

#' @export
format_lineage_genes <- function(cds, type = NULL){
  lineages <- names(cds@lineages)
  for (lineage in lineages) {
    # 1. Generate the filtered genes
    lin_genes <- .format_lineage_specific_genes(lineage, cds, type = type)

    # 2. Get the existing per-lineage node from the SAME object we are updating
    #    (was `cds_F` before -- a stray global reference).
    lineage_genes_orig <- cds@lineage_genes[[lineage]]

    # 3. Update the object
    cds@lineage_genes[[lineage]] <- list(
      "lineage_genes" = lineage_genes_orig,
      "filtered"      = lin_genes
    )
  }
  cds
}

.format_lineage_specific_genes <- function(lineage, cds, p_cutoff = 0.05,
                                          FC_pattern_cutoff = 0.2, FC_diffend_cutoff = 0.4,
                                          dynamic_I_cutoff = 0.1, dynamic_p_cutoff = 0.05,
                                          threshold = 0.1, p_adjust = "BH",
                                          specificity = "high", dynamic_test = "Moran",
                                          type = NULL) {
  lineages = names(cds@lineages)

  tabs          <- .get_lineage_test_tables(cds, lineage, type = type)
  pattern_genes <- tabs$pattern
  diffend_genes <- tabs$diffend
  type          <- tabs$type   # resolved model type (may still be NULL for nb/par output)

  is_quasi <- !is.null(type) && type == "quasipoisson"

  #First filter out genes that express in fewer than 100 cells
  expressed_genes <- .filter_by_expression_lineage(cds=cds, lineage=lineage, mode = "number", N = 100, ratio = 0.01)
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
  #Third filter out genes with NA values in averageFC and global pvalue
  pattern_genes_adjusted <- pattern_genes[(!is.na(pattern_genes$average_FC_pattern) & !is.na(pattern_genes$pvalue_combined_pattern)),]
  diffend_genes_adjusted <- diffend_genes[(!is.na(diffend_genes$average_FC_diffend) & !is.na(diffend_genes$pvalue_combined_diffend)),]
  #post-hoc p values
  if(p_adjust != FALSE){
    pattern_genes_adjusted$pvalue_combined_pattern <- p.adjust(pattern_genes_adjusted$pvalue_combined_pattern, method = "BH")
    diffend_genes_adjusted$pvalue_combined_diffend <- p.adjust(diffend_genes_adjusted$pvalue_combined_diffend, method = "BH")
  }
  if(length(lineages) > 2){
    if(p_adjust != FALSE){
      pairwise_columns_p <- grep("^pvalue_.*vs.*_pattern$", colnames(pattern_genes_adjusted), value = TRUE)
      pattern_genes_adjusted[, pairwise_columns_p] <- lapply(pattern_genes_adjusted[, pairwise_columns_p, drop = FALSE], p.adjust, method = "BH")
      p_matrix_adj <- pattern_genes_adjusted[, pairwise_columns_p, drop = FALSE]
      pairwise_columns_d <- grep("^pvalue_.*vs.*_diffend$", colnames(diffend_genes_adjusted), value = TRUE)
      diffend_genes_adjusted[, pairwise_columns_d] <- lapply(diffend_genes_adjusted[, pairwise_columns_d, drop = FALSE], p.adjust, method = "BH")
      d_matrix_adj <- diffend_genes_adjusted[, pairwise_columns_d, drop = FALSE]
    }
    FCs_p = pattern_genes_adjusted[, grepl("log2FC", colnames(pattern_genes_adjusted))]
    median_FC_p <- apply(FCs_p, 1, median, na.rm = TRUE)
    median_FC_p = median_FC_p[rownames(pattern_genes_adjusted)]
    pattern_genes_adjusted$median_FC_pattern <- median_FC_p
    FCs_d = diffend_genes_adjusted[, grepl("log2FC", colnames(diffend_genes_adjusted))]
    median_FC_d = apply(FCs_d, 1, median, na.rm = TRUE)
    median_FC_d = median_FC_d[rownames(diffend_genes_adjusted)]
    diffend_genes_adjusted$median_FC_diffend <- median_FC_d
    if(specificity == "high"){
      lineage_genes_p <- pattern_genes_adjusted[
        rowSums(FCs_p >= FC_pattern_cutoff, na.rm = TRUE) == rowSums(!is.na(FCs_p)) &
          rowSums(p_matrix_adj <= p_cutoff, na.rm = TRUE) == rowSums(!is.na(p_matrix_adj)),
      ]
      if(is_quasi){
        final_df <- lineage_genes_p[4:ncol(lineage_genes_p)] %>%
          select(
            matches("^waldStat|^pvalue|^log2FC"),
            predictA = matches("^predictA") %>% head(1),
            starts_with("predictB"),
            matches("^average_FC|^median_FC")
          )
        lineage_genes_p <- cbind(lineage_genes_p[,1:2], final_df)
      }
      lineage_genes_d <- diffend_genes_adjusted[
        rowSums(FCs_d >= FC_diffend_cutoff, na.rm = TRUE) == rowSums(!is.na(FCs_d)) &
          rowSums(d_matrix_adj <= p_cutoff, na.rm = TRUE) == rowSums(!is.na(d_matrix_adj)),
      ]
      if(is_quasi){
        final_df <- lineage_genes_d[4:ncol(lineage_genes_d)] %>%
          select(
            matches("^waldStat|^pvalue|^log2FC"),
            predictA = matches("^predictA") %>% head(1),
            starts_with("predictB"),
            matches("^average_FC|^median_FC")
          )
        lineage_genes_d <- cbind(lineage_genes_d[,1:2], final_df)
      }
    }
  }
  else{
    pattern_genes_adjusted['median_FC_pattern'] <- pattern_genes_adjusted['average_FC_pattern']
    diffend_genes_adjusted['median_FC_diffend'] <- diffend_genes_adjusted['average_FC_diffend']
    if(is_quasi){
      if (p_adjust != FALSE){
        pairwise_columns_p <- grep("^pvalue_.*vs.*_pattern$", colnames(pattern_genes_adjusted), value = TRUE)
        pattern_genes_adjusted[, pairwise_columns_p] <- lapply(pattern_genes_adjusted[, pairwise_columns_p, drop = FALSE], p.adjust, method = "BH")
        p_matrix_adj <- pattern_genes_adjusted[, pairwise_columns_p, drop = FALSE]
        pairwise_columns_d <- grep("^pvalue_.*vs.*_diffend$", colnames(diffend_genes_adjusted), value = TRUE)
        diffend_genes_adjusted[, pairwise_columns_d] <- lapply(diffend_genes_adjusted[, pairwise_columns_d, drop = FALSE], p.adjust, method = "BH")
        d_matrix_adj <- diffend_genes_adjusted[, pairwise_columns_d, drop = FALSE]
      }
      lineage_genes_p = pattern_genes_adjusted[pattern_genes_adjusted[, pairwise_columns_p] <= p_cutoff & pattern_genes_adjusted$median_FC_pattern >= FC_pattern_cutoff, ]
      lineage_genes_d = diffend_genes_adjusted[diffend_genes_adjusted[, pairwise_columns_d] <= p_cutoff & diffend_genes_adjusted$median_FC_diffend >= FC_diffend_cutoff, ]
    }else{
      lineage_genes_p = pattern_genes_adjusted[pattern_genes_adjusted$pvalue_combined_pattern <= p_cutoff & pattern_genes_adjusted$median_FC_pattern >= FC_pattern_cutoff, ]
      lineage_genes_d = diffend_genes_adjusted[diffend_genes_adjusted$pvalue_combined_diffend <= p_cutoff & diffend_genes_adjusted$median_FC_diffend >= FC_diffend_cutoff, ]
    }
  }
  gene = union(rownames(lineage_genes_p), rownames(lineage_genes_d))
  p_df <- pattern_genes_adjusted[gene, c("median_FC_pattern"), drop = FALSE]
  d_df <- diffend_genes_adjusted[gene, c("median_FC_diffend"), drop = FALSE]
  p_df <- p_df[rownames(d_df), , drop = FALSE]
  merged <- cbind(p_df, d_df)
  merged$transientScore <-
    rank(-merged$median_FC_pattern, ties.method = "min")^2 + rank(-merged$median_FC_diffend, ties.method = "min")^2
  lineage_spec_genes <- merged[order(merged$transientScore), ]
  lineage_spec_genes$lineage <- rep(lineage, nrow(lineage_spec_genes))
  lineage_spec <- list("pattern_filtered" = lineage_genes_p, "pattern_prefiltered" = pattern_genes_adjusted, "diffend_filtered" = lineage_genes_d, "diffend_prefiltered" = diffend_genes_adjusted, "combined" = lineage_spec_genes)
  lineage_spec
}

# Pull $pattern_test / $diffend_test out of cds@lineage_genes[[lineage]],
# regardless of which pipeline populated it. Returns list(pattern, diffend, type).
.get_lineage_test_tables <- function(cds, lineage, type = NULL) {
  node <- cds@lineage_genes[[lineage]]

  # Descend through an optional "lineage_genes" wrapper.
  if (is.list(node) && !is.null(node[["lineage_genes"]])) {
    node <- node[["lineage_genes"]]
  }

  # Descend through an optional model-type level (e.g. "quasipoisson", "nb").
  if (is.list(node) && is.null(node[["pattern_test"]])) {
    if (!is.null(type) && !is.null(node[[type]])) {
      node <- node[[type]]
    } else if (length(node) >= 1 && is.list(node[[1]]) &&
               !is.null(node[[1]][["pattern_test"]])) {
      # Fall back to the first (usually only) model type present.
      if (is.null(type)) type <- names(node)[1]
      node <- node[[1]]
    }
  }

  if (is.null(node[["pattern_test"]]) || is.null(node[["diffend_test"]])) {
    stop(sprintf(
      "Could not locate pattern_test/diffend_test for lineage '%s'. Inspect with: str(cds@lineage_genes[['%s']], max.level = 3)",
      lineage, lineage))
  }

  list(pattern = node[["pattern_test"]],
       diffend = node[["diffend_test"]],
       type    = type)
}

.find_branches <- function(cds, drop_non_binary = TRUE){
  graph_list = cds@graphs
  start = .find_start_point(graph_list)
  branch_points = .get_branch_points(graph_list)

  combined_graph <- graph.empty(directed = FALSE)
  for (g in graph_list) combined_graph <- igraph::union(combined_graph, g)
  combined_graph <- simplify(combined_graph)

  distances <- sapply(branch_points, function(bp) {
    sp <- suppressWarnings(shortest.paths(combined_graph, v = start, to = bp))
    sp[1, 1]
  })
  branch_points <- branch_points[order(distances)]

  branch_list_full = list()
  kept_vertices    = character(0)
  j_kept <- 0

  for (j in seq_along(branch_points)) {
    branch_vertex <- branch_points[j]

    # Lineage graphs that actually pass through this vertex.
    keep_idx <- vapply(graph_list, function(g) branch_vertex %in% V(g)$name, logical(1))
    graph_list_f <- graph_list[keep_idx]

    graph_list_trunc <- lapply(graph_list_f, .get_subgraph_opposite_to_start,
                               start = start, branch_vertex = branch_vertex)
    branch_list <- .group_graphs_by_vertex_overlap(graph_list_trunc)
    k <- length(branch_list)

    if (k != 2) {
      msg <- sprintf(
        "Branch point '%s' produced %d downstream branch(es) (expected 2); %d lineage(s) pass through it.",
        branch_vertex, k, length(graph_list_f))
      if (drop_non_binary) {
        warning(msg, " Skipping it.", call. = FALSE)
        next
      } else {
        warning(msg, " Keeping it (downstream binary-branch code may fail).", call. = FALSE)
      }
    }

    j_kept <- j_kept + 1
    names(branch_list) <- paste0("BP", j_kept, "_B", seq_along(branch_list))
    branch_list_full[[branch_vertex]] <- branch_list
    kept_vertices <- c(kept_vertices, branch_vertex)
  }

  if (length(branch_list_full) == 0) {
    warning("No branch points split cleanly into two branches. Returning empty list.",
            call. = FALSE)
    return(branch_list_full)
  }

  names(branch_list_full) <- paste0("BP_", seq_along(branch_list_full))
  branch_list_full
}

.find_start_point <-function(graph_list){
  start_ends = c()
  for(graph_name in names(graph_list)){
    graph = graph_list[[graph_name]]
    start_end = V(graph)[degree(graph) == 1]$name
    start_ends = c(start_ends, start_end)
  }
  start = names(sort(table(start_ends),decreasing=TRUE)[1])
  start
}

.get_subgraph_opposite_to_start <- function(graph, branch_vertex, start) {
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

.get_branch_points <- function(graph_list){
  combined_graph <- graph.empty(directed = FALSE)
  for (g in graph_list) {
    combined_graph <- igraph::union(combined_graph, g)
  }
  combined_graph <- simplify(combined_graph)
  branch_points <- V(combined_graph)[degree(combined_graph) >= 3]
  branch_points = branch_points$name
  branch_points
}

.group_graphs_by_vertex_overlap <- function(graph_list) {
  graph_names <- names(graph_list)
  if (is.null(graph_names)) {
    graph_names <- paste0("graph_", seq_along(graph_list))
    names(graph_list) <- graph_names
  }

  n <- length(graph_list)

  # Degenerate sizes: 0 or 1 graph can't form an overlap matrix loop.
  if (n == 0) return(list())
  if (n == 1) return(setNames(list(graph_names), "1"))

  for (i in seq_along(graph_list)) {
    g <- graph_list[[i]]
    if (is.null(V(g)$name)) V(g)$name <- as.character(V(g))
  }

  vertex_sets <- lapply(graph_list, function(g) V(g)$name)

  overlap_matrix <- matrix(0, nrow = n, ncol = n,
                           dimnames = list(graph_names, graph_names))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      if (length(intersect(vertex_sets[[i]], vertex_sets[[j]])) > 0) {
        overlap_matrix[i, j] <- 1
        overlap_matrix[j, i] <- 1
      }
    }
  }

  overlap_graph <- graph_from_adjacency_matrix(overlap_matrix,
                                               mode = "undirected", diag = FALSE)
  comps <- components(overlap_graph)
  split(names(comps$membership), comps$membership)
}
