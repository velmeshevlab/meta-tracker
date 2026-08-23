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

#' Find binary branch points across a set of lineage trajectories
#'
#' Detects branch points where lineages diverge and groups the downstream arms.
#' Branch points that don't split cleanly into exactly two arms are skipped
#' (with a warning) unless \code{drop_non_binary = FALSE}.
#'
#' @param cds A \code{metatracker_data_set} with \code{@graphs} populated.
#' @param drop_non_binary Skip non-binary branch points (default TRUE).
#' @return A named list of branch points; each element is a list of two arms,
#'   each arm a character vector of lineage names.
#' @export
find_branches <- function(cds, drop_non_binary = TRUE) {
  .find_branches(cds, drop_non_binary = drop_non_binary)
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


# --- Branch-specific genes ----------------------------------------------------

#' Branch-specific genes for one arm of a branch point
#'
#' Genes specifically up in the lineages of one arm (\code{branch_number}) of a
#' branch point relative to the other arm, using the pattern-test fold-changes
#' from \code{lineage_specific_genes_par()} and (optionally) the Moran dynamic
#' genes from \code{within_lineage_DE_par()}.
#'
#' @param branch_point One element of \code{find_branches(cds)} (a list of two
#'   arms, each a character vector of lineage names).
#' @param cds A \code{metatracker_data_set} with \code{@lineage_genes} populated.
#' @param branch_number Which arm (1 or 2) to test.
#' @param p_cutoff,FC_cutoff Adjusted-p and log2FC thresholds for specificity.
#' @param dynamic_p_cutoff,dynamic_FC_cutoff Dynamic-gene thresholds; set
#'   \code{dynamic_FC_cutoff = FALSE} to skip the dynamic filter.
#' @param p_adjust Method for \code{p.adjust} (or FALSE to skip).
#' @param dynamic_test "Moran" (uses \code{$I}/\code{$padj}) or otherwise
#'   \code{$scaled_FC}.
#' @return A data.frame (meta_p, median_FC, branch), ordered by |median_FC|.
#' @export
format_branch_specific_genes <- function(branch_point, cds, branch_number = 1,
                                         p_cutoff = 0.05, FC_cutoff = 0.2,
                                         dynamic_p_cutoff = 0.05, dynamic_FC_cutoff = 0.1,
                                         p_adjust = "BH", dynamic_test = "Moran") {
  branches_1 <- branch_point[[branch_number]]
  branches_2 <- if (branch_number == 1) branch_point[[2]] else branch_point[[1]]

  branch_genes      <- list()
  branch_gene_names <- list()
  for (lineage in branches_1) {
    lineage_genes <- .get_lineage_test_tables(cds, lineage)$pattern
    # 1. drop genes expressed in fewer than 100 cells
    expressed_genes <- .filter_by_expression_lineage(cds = cds, lineage = lineage,
                                                     mode = "number", N = 100, ratio = 0.01)
    common_genes  <- intersect(expressed_genes, rownames(lineage_genes))
    lineage_genes <- lineage_genes[common_genes, , drop = FALSE]
    # 2. keep only dynamic genes (Moran I / padj), if requested and available
    if (!isFALSE(dynamic_FC_cutoff)) {
      dynamic <- cds@dynamic_genes[[lineage]]
      if (is.null(dynamic)) {
        warning(sprintf(paste0("No @dynamic_genes for lineage '%s'; skipping the dynamic ",
                               "filter. Run within_lineage_DE_par() first, or set ",
                               "dynamic_FC_cutoff = FALSE."), lineage), call. = FALSE)
      } else {
        dynamic_genes <- if (dynamic_test == "Moran")
          rownames(dynamic[dynamic$I >= dynamic_FC_cutoff & dynamic$padj <= dynamic_p_cutoff, ])
        else
          rownames(dynamic[dynamic$scaled_FC >= dynamic_FC_cutoff, ])
        common_genes  <- intersect(dynamic_genes, rownames(lineage_genes))
        lineage_genes <- lineage_genes[common_genes, , drop = FALSE]
      }
    }
    FC_names <- paste0("log2FC_", lineage, "vs", branches_2, "_pattern")
    p_names  <- paste0("pvalue_", lineage, "vs", branches_2, "_pattern")
    lineage_genes_filtered <- lineage_genes[, c(p_names, FC_names)]
    FCs      <- lineage_genes_filtered[, grepl("log2FC", colnames(lineage_genes_filtered))]
    p_values <- lineage_genes_filtered[, grepl("pvalue", colnames(lineage_genes_filtered))]
    if (is.vector(FCs)) {
      # single "other" lineage: columns are (meta_p, averageFC)
      colnames(lineage_genes_filtered) <- c("meta_p", "averageFC")
      lineage_genes_filtered <- lineage_genes_filtered[(!is.na(lineage_genes_filtered$meta_p) &
                                                         !is.na(lineage_genes_filtered$averageFC)), ]
      lineage_genes_filtered$p_adjusted <- p.adjust(lineage_genes_filtered$meta_p, method = p_adjust)
      lineage_spec_genes <- lineage_genes_filtered[lineage_genes_filtered$averageFC >= FC_cutoff &
                                                     lineage_genes_filtered$p_adjusted <= p_cutoff, ]
    } else {
      average_FC <- apply(FCs, 1, .get_average_FC)
      meta_p     <- apply(p_values, 1, .get_meta_p)
      ref <- rownames(lineage_genes_filtered); col <- colnames(lineage_genes_filtered)
      average_FC <- average_FC[ref]; meta_p <- meta_p[ref]
      lineage_genes_filtered <- cbind(lineage_genes_filtered, average_FC, meta_p)
      colnames(lineage_genes_filtered) <- c(col, "average_FC", "meta_p")
      lineage_genes_filtered <- lineage_genes_filtered[(!is.na(lineage_genes_filtered$meta_p) &
                                                         !is.na(lineage_genes_filtered$average_FC)), ]
      lineage_genes_filtered$p_adjusted <- p.adjust(lineage_genes_filtered$meta_p, method = p_adjust)
      FCs <- FCs[rownames(lineage_genes_filtered), ]
      lineage_spec_genes <- lineage_genes_filtered[rowSums(FCs >= FC_cutoff) == ncol(FCs) &
                                                     lineage_genes_filtered$p_adjusted <= p_cutoff, ]
    }
    branch_genes[[lineage]]      <- lineage_spec_genes
    branch_gene_names[[lineage]] <- rownames(lineage_spec_genes)
  }

  # genes specific in EVERY lineage of the arm
  branch_gene_names <- Reduce(intersect, branch_gene_names)
  meta_p_matrix    <- matrix(0, nrow = length(branch_gene_names), ncol = 0)
  median_FC_matrix <- matrix(0, nrow = length(branch_gene_names), ncol = 0)
  for (lineage in branches_1) {
    lineage_genes <- branch_genes[[lineage]]
    FC_names <- paste0("log2FC_", lineage, "vs", branches_2, "_pattern")
    p_names  <- paste0("pvalue_", lineage, "vs", branches_2, "_pattern")
    if (ncol(lineage_genes) <= 3) {
      FCs      <- lineage_genes[branch_gene_names, "averageFC", drop = FALSE]
      p_values <- lineage_genes[branch_gene_names, "meta_p",   drop = FALSE]
    } else {
      FCs      <- lineage_genes[branch_gene_names, FC_names]
      p_values <- lineage_genes[branch_gene_names, p_names]
    }
    meta_p_matrix    <- cbind(meta_p_matrix, p_values)
    median_FC_matrix <- cbind(median_FC_matrix, FCs)
  }
  if (ncol(meta_p_matrix) == 1) {
    meta_p    <- setNames(meta_p_matrix[, 1], rownames(meta_p_matrix))
    median_FC <- setNames(median_FC_matrix[, 1], rownames(median_FC_matrix))
  } else {
    meta_p    <- apply(meta_p_matrix, 1, .get_meta_p)
    median_FC <- apply(median_FC_matrix, 1, median)
  }
  if (!isFALSE(p_adjust)) meta_p <- p.adjust(meta_p, method = p_adjust)
  meta_p <- meta_p[names(median_FC)]

  final_out <- as.data.frame(cbind(meta_p, median_FC))
  rownames(final_out) <- names(median_FC)
  final_out$branch <- rep(names(branch_point)[branch_number], nrow(final_out))
  colnames(final_out) <- c("meta_p", "median_FC", "branch")
  final_out[with(final_out, order(-abs(median_FC), meta_p)), ]
}

#' Compute branch-specific genes for every branch point and store them on the cds
#'
#' One-call wrapper: detects branch points with \code{find_branches()}, computes
#' branch-specific genes for each arm, and stores the tables in the cds under
#' \code{metadata(cds)$branch_genes} (a named list keyed \code{<BP>_B<arm>}).
#' Retrieve with \code{metadata(cds)$branch_genes}. Optionally also writes .tsv.
#'
#' @param cds A \code{metatracker_data_set} with \code{@lineage_genes} populated
#'   (run \code{lineage_specific_genes_par()}; and \code{within_lineage_DE_par()}
#'   if using the Moran dynamic filter).
#' @param branch_numbers Which arms to compute per branch point (default both).
#' @param drop_non_binary Passed to \code{find_branches()}.
#' @param write If TRUE, also write \code{<BP>_B<arm>_branch_genes.tsv}.
#' @param out_dir Directory for the .tsv files when \code{write = TRUE}.
#' @param ... Thresholds forwarded to \code{format_branch_specific_genes()}.
#' @return \code{cds} with \code{metadata(cds)$branch_genes} populated.
#' @export
find_branch_genes <- function(cds, branch_numbers = c(1, 2), drop_non_binary = TRUE,
                              write = FALSE, out_dir = ".", ...) {
  branches <- .find_branches(cds, drop_non_binary = drop_non_binary)
  if (length(branches) == 0) {
    warning("No branch points found; nothing to store.", call. = FALSE)
    return(cds)
  }
  if (write && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  out <- list()
  for (index in seq_along(branches)) {
    branch_point <- branches[[index]]
    for (bn in branch_numbers) {
      if (bn > length(branch_point)) next
      key <- paste0(names(branches)[index], "_B", bn)
      tab <- format_branch_specific_genes(branch_point, cds, branch_number = bn, ...)
      out[[key]] <- tab
      if (write) {
        utils::write.table(tab, file = file.path(out_dir, paste0(key, "_branch_genes.tsv")),
                           sep = "\t", quote = FALSE)
      }
    }
  }

  md <- metadata(cds)
  md[["branch_genes"]] <- out
  metadata(cds) <- md
  cds
}# Branch-point detection and lineage-specific gene formatting.

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
