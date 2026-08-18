# Expression-based gene filtering.

#' @export
filter_by_expression <- function(cds, mode = "number", N = 100, ratio = 0.01){
  data = counts(cds)
  lineages = names(cds@lineages)
  all.expressed_genes = c()
  for(lineage in lineages){
    cells = as.character(cds@lineages[[lineage]])
    data.sub = data[,cells]
    if(mode == "ratio"){
     cutoff = ratio*ncol(data.sub)
    }
    else{
     cutoff = N
    }
    expressed_genes = rownames(data.sub)[rowSums(data.sub> 0) > cutoff]
    all.expressed_genes = append(all.expressed_genes, expressed_genes)
  }
  all.expressed_genes = unique(all.expressed_genes)
  all.expressed_genes
}

.filter_by_expression_lineage <- function(cds, lineage, mode = "number", N = 100, ratio = 0.01){
  data = counts(cds)
  cells = as.character(cds@lineages[[lineage]][['name']])
  data.sub = data[,cells]
  if(mode == "ratio"){
    cutoff = ratio*ncol(data.sub)
  }
  else{
    cutoff = N
  }
  expressed_genes = rownames(data.sub)[rowSums(data.sub> 0) > cutoff]
  expressed_genes
}
