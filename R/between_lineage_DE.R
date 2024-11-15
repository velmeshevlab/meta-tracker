calculate_weight_matrix <- function(cds)
  {
  lineages = names(cds@lineages)
  all_cells = c()
  lineage_list = list()
  i = 1
  for(lineage in lineages){
    all_cells = c(all_cells, cds@lineages[[lineage]])
    lineage_list[[i]] <- cds@lineages[[lineage]]
    i <- i + 1
    }
  all_cells = unique(all_cells)
  weight_matrix <- matrix(0, nrow = length(all_cells), ncol = length(lineages))
  rownames(weight_matrix) <- all_cells
  colnames(weight_matrix) <- names(cds@lineages)
  for (i in seq_along(all_cells)) {
    member <- all_cells[i]
    # Check membership in each list
    membership <- sapply(lineage_list, function(lst) member %in% lst)
    # Determine the number of lists the member belongs to
    num_lists <- sum(membership)
    # Assign weights based on membership
    weight_matrix[i, membership] <- 1 / num_lists
    }
  return(weight_matrix)
  }

lineage_specific_genes <- function(cds, cellWeights){
  lineages = names(cds@lineages)
  counts = matrix(,nrow = ncol(cds@expression[[lineages[1]]])-7,ncol = 0)
  for(lineage in lineages){
    d = cds@expression[[lineage]]
    d = t(as.matrix(sapply(d[,8:ncol(d)], as.numeric)))
    counts <- cbind(counts, d)
    }
  gamlist = tradeSeq::fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, U = U, nknots = nknots)
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
