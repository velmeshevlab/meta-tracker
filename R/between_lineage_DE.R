calculate_dynamic_FC <- function(cds, lineage, genes, lineages){
  genes = genes[genes %in% colnames(cds@expectation[[lineage]])]
  FC_matrix = sapply(genes, calculate_dynamic_FC_gene, cds = cds, test_lineage = lineage, lineages = lineages)
  rownames(FC_matrix) <- names(cds@lineages)[names(cds@lineages) != lineage]
  FC_matrix
  }

calculate_dynamic_FC_gene <- function(gene, cds, test_lineage, lineages){
  exp1 = cds@expectation[[lineage]]
  pt = c(1:nrow(exp1))
  FCs = c()
  for(lin in names(lineages)){
    if(lin != lineage){
      exp2 = cds@expectation[[lin]]
      auc_dataset1 <- trapz(pt, exp1[,gene])
      auc_dataset2 <- trapz(pt, exp2[,gene])
      auc_difference <- log2(auc_dataset1/auc_dataset2)
      FCs <- c(FCs, auc_difference)
      }
    }
  FCs
}

lineage_specific_genes <- function(cds, test_lineage, U = NULL, nknots = 6, parallel = F, p_cutoff = 0.05, FC_cutoff = 1){
    lineages = names(cds@lineages)
    lineage_genes = get_lineage_genes(cds, test_lineage, U = U, nknots = nknots, parallel = parallel, p_cutoff = p_cutoff, FC_cutoff = FC_cutoff, lineages = lineages)
    cds@lineage_genes[[test_lineage]] <- lineage_genes
    cds
}

branch_specific_genes <- function(cds, test_lineages, name, U = NULL, nknots = 6, parallel = F, p_cutoff = 0.05, FC_cutoff = 1){
  lineages = names(cds@lineages)
  genes = c()
  for(test_lineage in test_lineages){
    dynamic_genes = rownames(cds@dynamic_genes[[test_lineage]])
    genes = c(genes, dynamic_genes)
    }
  genes = unique(genes)
  for(lineage in lineages){
    genes = genes[genes %in% colnames(cds@expectation[[lineage]])]
    }
  all_Ps <- matrix(0, nrow = length(genes), ncol = 0)
  all_FCs <- matrix(0, nrow = length(genes), ncol = 0)
  for(test_lineage in test_lineages){
    lineage_genes = get_lineage_genes(cds, test_lineage, genes = genes, U = U, nknots = nknots, parallel = parallel, p_cutoff = p_cutoff, FC_cutoff = FC_cutoff, lineages = c(test_lineage, lineages[!(lineages %in% test_lineages)]))
    lineage_genes = lineage_genes[,1:(ncol(lineage_genes)-2)]
    Ps = lineage_genes[,1:(ncol(lineage_genes)/2)]
    FCs = lineage_genes[,((ncol(lineage_genes)/2)+1):ncol(lineage_genes)]
    all_Ps = cbind(all_Ps, Ps)
    all_FCs = cbind(all_FCs, FCs)
    }
  res = cbind(all_Ps, all_FCs)
  rownames(res) <- genes
  cds@lineage_genes[[name]] <- res
  cds
}

get_lineage_genes <- function(cds, test_lineage, genes = NULL, U = NULL, nknots = 6, parallel = F, p_cutoff = 0.05, FC_cutoff = 1, lineages = NULL){
  counts = matrix(,nrow = ncol(cds@expression[[lineages[1]]])-7,ncol = 0)
  all_metacells = c()
  lineage_list = list()
  pt_list = list()
  i = 1
  for(lineage in lineages){
    metacells = paste0(lineage, "_", c(1:nrow(cds@expression[[lineage]])))
    d = cds@expression[[lineage]]
    d = t(as.matrix(sapply(d[,8:ncol(d)], as.numeric)))
    colnames(d) <- metacells
    counts <- cbind(counts, d)
    all_metacells <- c(all_metacells, metacells)
    lineage_list[[i]] <- metacells
    pt = cds@pseudotime[[lineage]][,1]
    names(pt) <- metacells
    pt_list[[i]] <- pt
    i <- i + 1
    }
  dynamic_genes = rownames(cds@dynamic_genes[[test_lineage]])
  if(length(genes) > 0){
    counts = counts[dynamic_genes,]
  }
  else{
    counts = counts[genes,]
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
  gamlist = tradeSeq::fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, U = U, nknots = nknots, parallel = parallel)
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
  }
  else{
    p_values = as.matrix(res[,"pvalue"])
    colnames(p_values) <- paste0("pvalue_", test_lineage, "vs", lineages[lineages != test_lineage])
  }
  p_values_sel = p_values[rowSums(p_values < p_cutoff) == ncol(p_values), ]
  gene_names = rownames(p_values_sel)
  FCs = calculate_dynamic_FC(cds, test_lineage, gene_names, lineages)
  FCs = t(FCs)
  FCs_sel = FCs[rowSums(FCs >= FC_cutoff) == ncol(FCs), ]
  p_values_sel = p_values_sel[rownames(FCs_sel),]
  combined_pvalue <- apply(p_values_sel, 1, get_meta_p)
  average_FC = apply(FCs_sel, 1, average_FC)
  colnames_old = c(colnames(p_values_sel), colnames(FCs_sel))
  final_res = cbind(p_values_sel, FCs_sel, combined_pvalue, average_FC)
  colnames(final_res) <- c(colnames_old, c("meta_p", "average_FC"))
  final_res = final_res[with(final_res, order(meta_p, -average_FC)), ]
  final_res
}

average_FC <- function(FC){
  linear_values <- 2^FC
  mean_linear <- mean(linear_values)
  mean_log2 <- log2(mean_linear)
  mean_log2
}

get_meta_p <- function(Ps){
  nonZero = length(which(Ps!=0))
  if(nonZero >= 2){
    sumlog(Ps)$p
  }
  else{
    0
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
