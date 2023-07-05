#Load dependent packages
library(tradeSeq)

within_lineage_DE <- function(counts = expr_matrix,  #a matrix of expression values, with genes in rows and cells in columns
                          conditions = condition_vector, #a vector of condition information
                          nknots = 3, #number of knots used to fit the GAM
                          pseudotime = pseudotime_matrix, #a matrix of pseudotime values, each row for a cell
                          pairwise=TRUE #pairwise comparison between different conditions
                          ){
  cell_wt<-as.matrix(rep(1,ncol(expr_matrix)),ncol=1)
  gamList<-tradeSeq::fitGAM(counts=counts,conditions=conditions,nknots=nknots,pseudotime=pseudotime,cellWeights=cell_wt)
  res<-tradeSeq::conditionTest(gamList,global=TRUE,pairwise=pairwise,lineages=FALSE)
  res_full<-res[!is.na(res$pvalue),]
  res_sig<-res_full[res_full$pvalue<0.05,]
  return(res_sig)
}