within_lineage_DE <- function(cds, #metatracker object
                          lineage, #name of the lineage to analyze
                          conditions = NULL, #a vector of condition information
                          nknots = 3, #number of knots used to fit the GAM
                          pairwise=TRUE, #pairwise comparison between different conditions
                          contrast_type = "end",
                          p = 0.05, #p value threshold
                          parallel = F
                          ){
  d =cds@expression[[lineage]]
  #conditions = d[,conditions]
  counts = t(as.matrix(sapply(d[,8:ncol(d)], as.numeric))) #a matrix of expression values, with genes in rows and cells in columns
  pseudotime = as.matrix(cds@pseudotime[[lineage]]) #a matrix of pseudotime values, each row for a cell
  cell_wt<-as.matrix(rep(1,ncol(counts)),ncol=1)
  #gamList<-tradeSeq::fitGAM(counts=counts,conditions=conditions,nknots=nknots,pseudotime=pseudotime,cellWeights=cell_wt,sce=FALSE)
  #res<-tradeSeq::conditionTest(gamList,global=TRUE,pairwise=pairwise,lineages=FALSE)
  #res_full<-res[!is.na(res$pvalue),]
  #res_sig<-res_full[res_full$pvalue<p,]
  #p_value <- getSmootherPvalues(gamList)
  #p_value = p_value[p_value<p,]
  #statLineage <- getSmootherTestStats(gamList)
  #statLineage = statLineage[names(p_value),]
  #final_res = cbind(statLineage, p_value)
  #final_res= final_res[order(final_res[,1], decreasing = T),]
  gamList<-tradeSeq::fitGAM(counts=counts,conditions=conditions,nknots=nknots,pseudotime=pseudotime,cellWeights=cell_wt, parallel = parallel)
  res = associationTest(gamList, contrastType = contrast_type)
  exp_95 = matrix(,nrow = nrow(res),ncol = 0)
  for(lin in names(cds@lineages)){
    exp_lin = cds@expectation[[lin]]
    exp_lin = exp_lin[,rownames(res)]
    exp_95_lin = apply(exp_lin, 2, function(x) as.numeric(quantile(x, 0.95)))
    exp_95 = cbind(exp_95, exp_95_lin)
    }
  rownames(exp_95) <- rownames(res)
  exp_95_max = rowMax(exp_95)
  FC_factor = apply(cds@expectation[[lineage]], 2, function(x) as.numeric(quantile(x, 0.95)))/exp_95_max
  scaled_FC = res$meanLogFC*FC_factor
  res$scaled_FC <- scaled_FC                      
  res_sig = res[res$pvalue<p,]
  res_sig = res_sig[with(res_sig, order(pvalue, -meanLogFC)), ]
  cds@dynamic_genes[[lineage]] <- res_sig
  return(cds)
}
