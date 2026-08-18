# Between-lineage differential expression (tradeSeq).

#' @export
lineage_specific_genes_par <- function(cds, lineages = names(cds@lineages)){
  counts = matrix(,nrow = nrow(cds),ncol = 0)
  all_metacells = c()
  all_size_factor = c()
  lineage_list = list()
  pt_list = list()
  i = 1
  for(lineage in lineages){
    metacells = paste0(lineage, "_", c(1:nrow(cds@expression[[lineage]][['sum']])))
    d = cds@expression[[lineage]][['sum']]
    d = t(as.matrix(sapply(d[,7:ncol(d)], as.numeric)))
    d <- d[rownames(cds),]
    colnames(d) <- metacells
    counts <- cbind(counts, d)
    all_metacells <- c(all_metacells, metacells)
    size_factor = cds@expression[[lineage]][['sum']]$size_factor
    all_size_factor <- c(all_size_factor, size_factor)
    lineage_list[[i]] <- metacells
    #pt = cds@pseudotime[[lineage]]
    pt = cds@expression[[lineage]][["sum"]]$pseudotime
    pt <-  (pt - min(pt))/(max(pt) - min(pt))
    names(pt) <- metacells
    pt_list[[i]] <- pt
    i <- i + 1
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
  counts = counts[rownames(cds),]
  print(paste0("Testing ", nrow(counts), " genes"))
  gamlist = tradeSeq::fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, U = NULL, nknots = 6, offset = log(all_size_factor))
  pattern = tradeSeq::patternTest(models = gamlist, global = T, pairwise = T)
  diffend = tradeSeq::diffEndTest(models = gamlist, global = T, pairwise = T)
  out <- sapply(lineages, .lineage_specific_genes, cds = cds, model = gamlist, pattern = pattern, diffend = diffend, lineages = lineages, simplify = FALSE)
  cds@lineage_genes <- out
  cds 
}

.lineage_specific_genes <- function(test_lineage, cds, model, pattern, diffend, genes = NULL, lineages){
  print(paste0("Testing ", test_lineage))
  if(length(lineages) > 2){
    index = which(test_lineage == lineages)
    p_list = c()
    p_names = c()
    wd_list = c()
    fc_list = c()
    wd_names = c()
    fc_names = c()
    flip_fc <- c()
    for(lineage in lineages){
      index2 = which(lineage == lineages)
      if(index != index2){
        if(index<index2){
          p_name = paste0("pvalue_", index, "vs", index2)
          wd_name = paste0("waldStat_", index, "vs", index2)
          fc_name = paste0("logFC", index, "_", index2)
          flip_fc <- c(flip_fc, FALSE) 
        }
        else{
          p_name = paste0("pvalue_", index2, "vs", index)
          wd_name = paste0("waldStat_", index2, "vs", index)
          fc_name = paste0("logFC", index2, "_", index)
          flip_fc <- c(flip_fc, TRUE)
        }
        p_name_new = paste0("pvalue_", test_lineage, "vs", lineage)
        wd_name_new = paste0("waldStat_", test_lineage, "vs", lineage)
        fc_name_new = paste0("log2FC_", test_lineage, "vs", lineage)
        p_list <- c(p_list, p_name)
        wd_list = c(wd_list, wd_name)
        fc_list = c(fc_list, fc_name)
        p_names <- c(p_names, p_name_new)
        wd_names <- c(wd_names, wd_name_new)
        fc_names <- c(fc_names, fc_name_new)
      }
    }
    #pattern test result
    p_values_p = pattern[,p_list]
    wd_values = pattern[,wd_list]
    #combined_pvalue_pattern <- apply(p_values_p, 1, .get_meta_p)
    ref <- rownames(pattern)
    
    p_values_p <- p_values_p[ref, , drop = FALSE]
    wd_values  <- wd_values[ref, , drop = FALSE]
    #combined_pvalue_pattern <- combined_pvalue_pattern[ref]
    
    stopifnot(
      identical(ref, rownames(p_values_p)),
      identical(ref, rownames(wd_values))
      #identical(ref, names(combined_pvalue_pattern))
    )
    
    pattern_sel <- cbind(
      pattern[, c(1, 3), drop = FALSE],
      p_values_p,
      wd_values
      #combined_pvalue_pattern
    )
    colnames(pattern_sel) <- c("waldStat_combined", "pvalue_combined", p_names, wd_names)
    colnames(pattern_sel) <- paste0(colnames(pattern_sel), "_pattern")
    #diffend tes result
    p_values_d = diffend[,p_list]
    wd_values = diffend[,wd_list]
    fc_values = diffend[,fc_list]
    fc_values[, flip_fc] <- -fc_values[, flip_fc]
    #combined_pvalue_diffend <- apply(p_values_d, 1, .get_meta_p)
    #transfer from logFC to log2FC
    fc_values <- as.data.frame(fc_values / log(2))
    average_FC = apply(fc_values, 1, .get_average_FC)
    
    ref <- rownames(diffend)
    p_values_d <- p_values_d[ref, , drop = FALSE]
    wd_values  <- wd_values[ref, , drop = FALSE]
    fc_values  <- fc_values[ref, , drop = FALSE]
    #combined_pvalue_diffend <- combined_pvalue_diffend[ref]
    average_FC <- average_FC[ref]
    
    stopifnot(
      identical(ref, rownames(p_values_d)),
      identical(ref, rownames(wd_values)),
      identical(ref, rownames(fc_values)),
      #identical(ref, names(combined_pvalue_diffend)),
      identical(ref, names(average_FC))
    )
    
    diffend_sel <- cbind(
      diffend[, c(1, 3), drop = FALSE],
      p_values_d,
      wd_values,
      #combined_pvalue_diffend,
      fc_values,
      average_FC
    )
    colnames(diffend_sel) <- c("waldStat_combined", "pvalue_combined", p_names, wd_names, fc_names, "average_FC")
    colnames(diffend_sel) <- paste0(colnames(diffend_sel), "_diffend")
    pattern_sel <- pattern_sel[rownames(diffend_sel), ]
    #Calculate the FC for pattern_test
    gene_names = rownames(pattern_sel)
    FCs = .calculate_FC_tradeseq(models = model, test_lineage = test_lineage, lineages = lineages, fc_names = fc_names, genes = gene_names, N = 1000)
    FCs_sel <- FCs
    pattern_sel = pattern_sel[rownames(FCs_sel),]
    average_FC = apply(FCs_sel, 1, .get_average_FC)
    average_FC <- average_FC[rownames(FCs_sel)]
    if (!identical(rownames(FCs_sel), names(average_FC))) {
      stop("Error: rownames(FCs_sel) and names(average_FC) do not match.")
    }
    pattern_fin = cbind(pattern_sel, FCs_sel, average_FC)
    colnames(pattern_fin) <- c(colnames(pattern_sel), colnames(FCs_sel), "average_FC_pattern")
    pattern_fin = as.data.frame(pattern_fin)
    diffend_fin = as.data.frame(diffend_sel)
    final_res = list("pattern_test" = pattern_fin, "diffend_test" = diffend_fin)
    final_res
  }
  else{
    pattern_sel = pattern[,c("waldStat", "pvalue")]
    diffend_sel = diffend[,c("waldStat", "pvalue", "logFC1_2")]
    diffend_sel$logFC1_2 <- diffend_sel$logFC1_2 / log(2)
    if (test_lineage == lineages[2]) {
      diffend_sel$logFC1_2 <- -diffend_sel$logFC1_2
    }
    gene_names = rownames(pattern_sel)
    FCs = .calculate_FC_tradeseq(models = model, test_lineage = test_lineage, lineages = lineages, genes = gene_names, N = 1000)
    FCs_sel = FCs
    pattern_sel = pattern_sel[rownames(FCs_sel),]
    diffend_sel = diffend_sel[rownames(FCs_sel),]
    pattern_fin = as.data.frame(cbind(pattern_sel, FCs_sel))
    diffend_fin = as.data.frame(diffend_sel)
    colnames(pattern_fin) <- c("waldStat_combined_pattern", "pvalue_combined_pattern", "average_FC_pattern")
    colnames(diffend_fin) <- c("waldStat_combined_diffend", "pvalue_combined_diffend", "average_FC_diffend")
    lineage_genes = list("pattern_test" = pattern_fin, "diffend_test" = diffend_fin)
    lineage_genes
  }
}

.calculate_FC_tradeseq <- function(models, test_lineage, lineages, fc_names = NULL, genes, knots = NULL, N){
  dm <- colData(models)$tradeSeq$dm # design matrix
  X <- colData(models)$tradeSeq$X # linear predictor
  slingshotColData <- colData(models)$crv
  pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                       pattern = "pseudotime")]
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  beta <- betaMat[genes,]
  knotPoints <- S4Vectors::metadata(models)$tradeSeq$knots
  lineage_map <- data.frame(
    lineage_id = 1:length(lineages),
    lineage_name = lineages
  )
  lookup <- setNames(lineage_map$lineage_id, lineage_map$lineage_name)
  time_points <- seq(from = 0, to = 1, length.out = N)
  #For test lineage
  lineage_id <- unname(lookup[test_lineage])
  df <- .getPredictRangeDf(dm, lineage_id, nPoints = N)
  Xdf <- predictGAM(lpmatrix = X,
                    df = df,
                    pseudotime = pseudotime)
  yhat_mat_test <- exp(Xdf %*% t(beta) + df$offset)
  auc_vector_test <- apply(yhat_mat_test, 2, function(y) trapz(time_points, y))
  if(!is.null(knots)){
    t1 <- unname(knotPoints[knots[[1]]])
    t2 <- unname(knotPoints[knots[[2]]])
    indices <- which(time_points >= t1 & time_points <= t2)
    yhat_mat_test <- yhat_mat_test[indices,]
    auc_vector_test <- apply(yhat_mat_test, 2, function(y) trapz(time_points[indices], y))
  }
  all_auc_log_diff <- c()
  for(lin in lineages){
    if (lin != test_lineage){
      lineage_id <- unname(lookup[lin])
      df <- .getPredictRangeDf(dm, lineage_id, nPoints = N)
      Xdf <- predictGAM(lpmatrix = X,
                        df = df,
                        pseudotime = pseudotime)
      yhat_mat <- exp(Xdf %*% t(beta) + df$offset)
      auc_vector <- apply(yhat_mat, 2, function(y) trapz(time_points, y))
      if(!is.null(knots)){
        t1 <- unname(knotPoints[knots[[1]]])
        t2 <- unname(knotPoints[knots[[2]]])
        indices <- which(time_points >= t1 & time_points <= t2)
        yhat_mat <- yhat_mat[indices,]
        auc_vector <- apply(yhat_mat, 2, function(y) trapz(time_points[indices], y))
      }
      auc_log_diff <- log2(auc_vector_test/auc_vector)
      all_auc_log_diff <- cbind(all_auc_log_diff, auc_log_diff)
    }
  }
  if (length(lineages)==2){
    colnames(all_auc_log_diff) <- paste0("log2FC_", test_lineage, "vs", lineages[lineages!=test_lineage])
  }
  else{
    colnames(all_auc_log_diff) <- paste0(fc_names, "_pattern")
  }
  return(all_auc_log_diff)
}

.get_average_FC <- function(FC){
  linear_values <- 2^FC
  mean_linear <- mean(linear_values, na.rm = TRUE)
  mean_log2 <- log2(mean_linear)
  mean_log2
}

.get_meta_p <- function(Ps) {
  Ps <- Ps[!is.na(Ps)]
  Ps[Ps == 0] <- .Machine$double.xmin
  if (length(Ps) == 0) return(NA_real_)
  
  stat <- -2 * sum(log(Ps))
  df <- 2 * length(Ps)
  pchisq(stat, df=df, lower.tail=FALSE)
}
