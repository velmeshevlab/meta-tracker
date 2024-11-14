#Load dependent packages
library(tradeSeq)

within_lineage_DE <- function(cds, #metatracker object
                          lineage, #name of the lineage to analyze
                          conditions = NULL, #a vector of condition information
                          nknots = 3, #number of knots used to fit the GAM
                          pairwise=TRUE, #pairwise comparison between different conditions
                          contrast_type = "end",
                          p = 0.05 #p value threshold
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
  gamList<-tradeSeq::fitGAM(counts=counts,conditions=conditions,nknots=nknots,pseudotime=pseudotime,cellWeights=cell_wt)
  res = associationTest(gamList, contrastType = contrast_type)
  res_sig = res[res$pvalue<p,]
  res_sig = res_sig[order(res_sig$waldStat, decreasing = T),]
  return(res_sig)
}
