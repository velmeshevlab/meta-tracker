#Load dependent packages
library(tradeSeq)

within_lineage_DE <- function(cds, #metatracker object
                          lineage, #name of the lineage to analyze
                          conditions = NULL, #a vector of condition information
                          nknots = 3, #number of knots used to fit the GAM
                          pairwise=TRUE, #pairwise comparison between different conditions
                          p = 0.05 #p value threshold
                          ){
  counts = cds_new@expression[lineage] #a matrix of expression values, with genes in rows and cells in columns
  pseudotime = cds_new@pseudotime[lineage] #a matrix of pseudotime values, each row for a cell
  cell_wt<-as.matrix(rep(1,ncol(expr_matrix)),ncol=1)
  gamList<-tradeSeq::fitGAM(counts=counts,conditions=conditions,nknots=nknots,pseudotime=pseudotime,cellWeights=cell_wt)
  res<-tradeSeq::conditionTest(gamList,global=TRUE,pairwise=pairwise,lineages=FALSE)
  res_full<-res[!is.na(res$pvalue),]
  res_sig<-res_full[res_full$pvalue<,]
  return(res_sig)
}
