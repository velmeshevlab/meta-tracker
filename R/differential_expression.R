#scaled_FC is fold change of dynamic gene expression scaled based on 95 percentile of expression of the entire dataset.
#Helps to filter out genes that are dynamically expressed in the given lineage but are expressed at much lower level than in other lineages.

within_lineage_DE_par <- function(cds, parallel = F, BPPARAM = BPPARAM, test = "tradeSeq"){
lineages = names(cds@lineages)
if(test == "tradeSeq"){
out <- lapply(lineages, within_lineage_DE_trade, cds = cds, parallel = T, BPPARAM = BPPARAM)
}
else(
out <- lapply(lineages, within_lineage_DE_Moran, cds = cds, parallel = T, BPPARAM = BPPARAM)     
)
names(out) <- lineages
cds@dynamic_genes <- out
cds
}

make_time_nb <- function(n, k = 5) {
  stopifnot(n >= 2, k >= 1)
  k <- min(k, n - 1)
  nb <- vector("list", n)
  for (i in seq_len(n)) {
    lo <- max(1, i - k)
    hi <- min(n, i + k)
    nb[[i]] <- setdiff(lo:hi, i)
  }
  class(nb) <- "nb"
  attr(nb, "region.id") <- as.character(seq_len(n))
  nb
}

within_lineage_DE_Moran <- function(lineage, #name of the lineage to analyze
                          cds, #metatracker object
                          k = 5
                          ){
        d =cds@expression[[lineage]]
        expr = t(as.matrix(sapply(d[,4:ncol(d)], as.numeric))) #a matrix of expression values, with genes in rows and cells in columns
        expr = expr[rownames(cds),]
        n_time <- ncol(expr)
        nb <- make_time_nb(n_time, k = k)
        lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
        res <- pblapply(seq_len(nrow(expr)), function(i) {
                  mt <- moran.test(expr[i, ], lw, zero.policy = TRUE)
                  c(I = unname(mt$estimate[["Moran I statistic"]]),
                  p = mt$p.value)
        })
        res <- as.data.frame(do.call(rbind, res))
        res$padj <- p.adjust(res$p, method = "fdr")
        rownames(res) <- rownames(expr)
        res <- res[order(res$padj, -res$I),]
        res
}

within_lineage_DE_trade <- function(lineage, #name of the lineage to analyze
                          cds, #metatracker object
                          conditions = NULL, #a vector of condition information
                          nknots = 3, #number of knots used to fit the GAM
                          pairwise=TRUE, #pairwise comparison between different conditions
                          contrast_type = "end",
                          parallel = F,
                          BPPARAM = F
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
  gamList<-tradeSeq::fitGAM(counts=counts,conditions=conditions,nknots=nknots,pseudotime=pseudotime,cellWeights=cell_wt, parallel = parallel, BPPARAM =BPPARAM)
  res = associationTest(gamList, contrastType = contrast_type)
  exp_95 = matrix(,nrow = nrow(res),ncol = 0)
  for(lin in names(cds@lineages)){
    exp_lin = cds@expectation[[lin]]
    exp_lin = exp_lin[,rownames(res)]
    exp_95_lin = apply(exp_lin, 2, function(x) as.numeric(quantile(x, 0.95, na.rm= TRUE)))
    exp_95 = cbind(exp_95, exp_95_lin)
    }
  rownames(exp_95) <- rownames(res)
  exp_95_max = apply(exp_95, 1, max, na.rm = TRUE)
  expectation = cds@expectation[[lineage]]
  expectation = expectation[,rownames(cds)]
  FC_factor = apply(expectation, 2, function(x) as.numeric(quantile(x, 0.95, na.rm= TRUE)))/exp_95_max
  scaled_FC = res$meanLogFC*FC_factor
  res$scaled_FC <- scaled_FC                      
  res = res_sig[with(res, order(pvalue, -scaled_FC)), ]
  res
}
