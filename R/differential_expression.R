#scaled_FC is fold change of dynamic gene expression scaled based on 95 percentile of expression of the entire dataset.
#Helps to filter out genes that are dynamically expressed in the given lineage but are expressed at much lower level than in other lineages.

within_lineage_DE_par <- function(cds, cores){
lineages = names(cds@lineages)
out <- mclapply(lineages, within_lineage_DE, cds = cds, mc.cores = cores)
names(out) <- lineages
cds@dynamic_genes <- out
cds
}
        
within_lineage_DE <- function(lineage, #name of the lineage to analyze
                          cds, #metatracker object
                          conditions = NULL, #a vector of condition information
                          nknots = 3, #number of knots used to fit the GAM
                          pairwise=TRUE, #pairwise comparison between different conditions
                          contrast_type = "end",
                          p = 0.05, #p value threshold
                          FC = 0.5, #fold change value threshold
                          parallel = F
                          ){
  print(paste0("Testing lineage ", lineage))
  d =cds@expression[[lineage]]
  #conditions = d[,conditions]
  counts = t(as.matrix(sapply(d[,4:ncol(d)], as.numeric))) #a matrix of expression values, with genes in rows and cells in columns
  counts = counts[rownames(cds),]
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
  expectation = cds@expectation[[lineage]]
  expectation = expectation[rownames(cds),]
  FC_factor = apply(expectation, 2, function(x) as.numeric(quantile(x, 0.95)))/exp_95_max
  scaled_FC = res$meanLogFC*FC_factor
  res$scaled_FC <- scaled_FC                      
  res_sig = res[res$pvalue<p & res$scaled_FC >= FC,]
  res_sig = res_sig[with(res_sig, order(pvalue, -scaled_FC)), ]
  res_sig
  #cds@dynamic_genes[[lineage]] <- res_sig
  #return(cds)
}
