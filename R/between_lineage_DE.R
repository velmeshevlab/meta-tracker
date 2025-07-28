lineage_specific_genes <- function(test_lineage, cds, U = NULL, nknots = 6, dyn_FC_cutoff = 0){
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
    lineage_genes = get_lineage_genes_v2(cds, test_lineage, U = U, nknots = nknots, lineages = lineages, dyn_FC_cutoff = dyn_FC_cutoff)
    if(length(lineage_genes) > 0){
    lineage_genes = lineage_genes[with(lineage_genes, order(meta_p, -abs(average_FC))), ]
    lineage_genes
    }
    else{return(NULL)}
}

get_lineage_genes_v2 <- function(cds, test_lineage, genes = NULL, U = NULL, nknots = 6, lineages = NULL, dyn_FC_cutoff = 0){
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
    counts = counts[all_dynamic_genes ,]
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
  gamlist = tradeSeq::fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, U = U, nknots = nknots)
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
    p_values_sel = as.matrix(p_values)
    if(nrow(p_values_sel) == 0){
      return(NULL)
    }
    gene_names = rownames(p_values_sel)
    FCs = calculate_dynamic_FC(cds, test_lineage, gene_names, lineages)
    FCs = t(FCs)
    FCs_sel = FCs
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
    res_sel = res
    p_values_sel = as.matrix(res_sel[,"pvalue"])
    rownames(p_values_sel) <- rownames(res_sel)
    colnames(p_values_sel) <- paste0("pvalue_", test_lineage, "vs", lineages[lineages != test_lineage])
    gene_names = rownames(p_values_sel)
    FCs = calculate_dynamic_FC_single(cds, test_lineage, gene_names, lineages[lineages != test_lineage])
    FCs_sel = FCs
    p_values_sel = p_values_sel[names(FCs_sel),]
    final_res = as.data.frame(cbind(p_values_sel, FCs_sel))
    colnames(final_res) <- c("meta_p", "average_FC")
    final_res = as.data.frame(final_res)
    final_res
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
  pt = c(1:nrow(exp1))
  exp2 = cds@expectation[[comp_lineage]]
  auc_dataset1 <- trapz(pt, exp1[,gene])
  auc_dataset2 <- trapz(pt, exp2[,gene])
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
  pt = c(1:nrow(exp1))
  FCs = c()
  for(lin in lineages){
    if(lin != test_lineage){
      exp2 = cds@expectation[[lin]]
      auc_dataset1 <- trapz(pt, exp1[,gene])
      auc_dataset2 <- trapz(pt, exp2[,gene])
      auc_difference <- log2(auc_dataset1/auc_dataset2)
      FCs <- c(FCs, auc_difference)
      }
    }
  FCs
}

lineage_specific_genes_par <- function(cds, cores, dyn_FC_cutoff = 0){
  lineages = names(cds@lineages)
  n.cores <- cores
  clust <- makeCluster(n.cores)
  out <- parSapply(clust, names(cds@lineages), lineage_specific_genes, cds = cds_new, dyn_FC_cutoff = dyn_FC_cutoff, simplify = FALSE)
  stopCluster(clust)
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
