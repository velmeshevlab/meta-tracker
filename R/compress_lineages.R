compress_3_2 <- function(df, leftover, n){
 df_comp = sw(df[1:(length(df) - (n + leftover))], n, mean, n)
 df_comp = c(df_comp, mean(df[(length(df[1:(length(df) - (n + leftover))])+ 1): length(df)]))
 return(df_comp)
}

filter_by_expression_lineage <- function(cds, lineage, mode = "number", N = 100, ratio = 0.01){
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

compress_lineage_v3_2 <- function(cds, lineage, method, N, cores = 1, ID){
  exp = compress_expression_v3_3(cds, lineage = lineage, method = method, N = N, cores = cores, ID = ID)
  cds@lineages[[lineage]] <- exp$lineage
  cds@expression[[lineage]] <- exp$expression
  cds@expectation[[lineage]] <- exp$expectation
  cds@pseudotime[[lineage]] <- exp$pseudotime
  cds
}

compress_expression_v3_3 <- function(cds, lineage, N, cores = 1, method = "sum", ID = FALSE){
  print("Updating pseudotime")
  #Extract the lineage object
  cds_sub <- get_lineage_object(cds, lineage)
  #Get pseudotime from the principal graph aux slot
  updated_pt <- cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  #Store back in cds@lineages as a list
  #need to add a new slot for updated pseudotime
  cell_barcodes <- cds@lineages[[lineage]]
  cds@lineages[[lineage]] <- list(
    name = cell_barcodes,
    updated_pt = updated_pt)
  if(lineage != FALSE){
    sel.cells = cds@lineages[[lineage]][['name']]
  }
  sel.cells = sel.cells[sel.cells %in% colnames(cds)]
  cds_subset = cds[,sel.cells]
  #preprare raw count matrix
  exp = as.data.frame(as.matrix(exprs(cds_subset)))
  exp_sum <- t(exp)
  exp_sum = exp_sum[,rownames(cds)]
  #prepare size factor
  size_factor <- (pData(cds_subset)[, 'Size_Factor'])
  size_factor <- size_factor[rownames(exp_sum)]
  #prepare pseudotime
  pt <- cds@lineages[[lineage]][['updated_pt']]
  pt <- as.data.frame(pt)
  pt <- pt[rownames(exp_sum), , drop = FALSE]
  colnames(pt) <- c("pseudotime")
  #prepare umap
  UMAP <- reducedDims(cds_subset)[["UMAP"]]
  UMAP <- UMAP[rownames(exp_sum),]
  if(ID == FALSE){
    if (N >= nrow(exp_sum)){
      stop(sprintf(
        "out of boundary: N (%d) must be smaller than number of rows in exp (%d)",
        N, nrow(exp_sum)
      ))
    }
    exp_sum = cbind(pt, UMAP, size_factor, exp_sum)
    exp_sum = exp_sum[order(exp_sum$pseudotime),]
    exp_sum$meta_cell <- cut(rank(exp_sum$pseudotime), breaks = N, labels = FALSE)
    gene_cols <- setdiff(
      colnames(exp_sum),
      c("pseudotime", "umap_1", "umap_2", "size_factor", "meta_cell")
    )
    print(paste0("Compressing lineage ", lineage, " with sum"))
    meta_sum <- exp_sum %>%
      group_by(meta_cell) %>%
      summarise(
        n_cells = n(),                     # number of cells in this meta-cell
        pseudotime = mean(pseudotime),     # mean pseudotime
        umap_1 = mean(umap_1),             # mean UMAP_x
        umap_2 = mean(umap_2),             # mean UMAP_y
        size_factor = sum(size_factor),    # sum size factors
        across(all_of(gene_cols), sum)   # sum gene counts
      ) %>%
      ungroup()
    meta_sum_ordered <- meta_sum[order(meta_sum$pseudotime), ]
    mat <- meta_sum_ordered[,7:(ncol(meta_sum_ordered))]
    #fit expression
    model <- expression ~ splines::ns(pseudotime, df = 7) + offset(log(size_factor))
    size_factor = meta_sum_ordered$size_factor
    d <-  (meta_sum_ordered$pseudotime - min(meta_sum_ordered$pseudotime))/(max(meta_sum_ordered$pseudotime)-min(meta_sum_ordered$pseudotime))
    print("Fitting curves scaled pseudotime")
    predict_pt <- seq(0, 1, length.out = N)
    mat_m <- as.matrix(mat)
    genes <- colnames(mat_m)
    fit_list_2 <- pbapply::pbsapply(setNames(seq_along(genes), genes),function(i) {
        fit.m3_3(exp.sel = mat_m[, i], pt = d, size_factor = size_factor, predict_pt = predict_pt, lineage = lineage, model = model, N = N)
      },
      cl = cores
    )
    fit_list_2 = apply(fit_list_2, 2, as.numeric)
    if (method != "sum") {
      exp_mean <- exp
      exp_mean = (t(exp_mean)) /  (pData(cds_subset)[, 'Size_Factor'])
      exp_mean = exp_mean[,rownames(cds)]
      pt <- pt[rownames(exp_mean), , drop = FALSE]
      UMAP <- UMAP[rownames(exp_mean),]
      exp_mean = cbind(pt, UMAP, exp_mean)
      exp_mean = exp_mean[order(exp_mean$pseudotime),]
      exp_mean$meta_cell <- cut(rank(exp_mean$pseudotime), breaks = N, labels = FALSE)
      gene_cols <- setdiff(
        colnames(exp_mean),
        c("pseudotime", "umap_1", "umap_2", "meta_cell")
      )
      print(paste0("Compressing lineage ", lineage, " with mean"))
      meta_mean <- exp_mean %>%
        group_by(meta_cell) %>%
        summarise(
          n_cells = n(),                     # number of cells in this meta-cell
          pseudotime = mean(pseudotime),     # mean pseudotime
          umap_1 = mean(umap_1),             # mean UMAP_x
          umap_2 = mean(umap_2),             # mean UMAP_y
          across(all_of(gene_cols), mean)   # mean gene counts
        ) %>%
        ungroup()
      meta_mean_ordered <- meta_mean[order(meta_mean$pseudotime), ]
      return(list(
        "lineage" = cds@lineages[[lineage]], "expression" = list("sum" = meta_sum_ordered, "mean" = meta_mean_ordered), "expectation" = fit_list_2, "pseudotime"  = list("real" = meta_sum_ordered$pseudotime, "scaled" = d)))
    }
    return(list("lineage"= cds@lineages[[lineage]], "expression" = meta_sum_ordered, "expectation" = fit_list_2, "pseudotime"  = list("real" = meta_sum_ordered$pseudotime, "scaled" = d)))
  }
}

fit.m3_3 <- function(exp.sel, pt, size_factor, predict_pt, lineage, model, N) {
  
  if (!requireNamespace("speedglm", quietly = TRUE)) {
    stop("speedglm not installed")
  }
  
  exp_data.sel <- data.frame(
    pseudotime  = as.numeric(pt),
    size_factor = as.numeric(size_factor),
    expression  = as.numeric(exp.sel)
  )
  
  # Try to fit the model and predict
  pred <- tryCatch({
    
    fit_model <- speedglm::speedglm(
      model,
      data = exp_data.sel,
      family = quasipoisson(),
      acc = 1e-3,
      model = FALSE,
      y = FALSE
    )
    
    newdata <- data.frame(
      pseudotime  = as.numeric(predict_pt),
      size_factor = 1
    )
    
    predict(fit_model, newdata = newdata, type = "response")
    
  }, error = function(e) {
    # If error occurs, return NA vector
    rep(NA_real_, N)
  })
  
  return(pred)
}

compress_expression_v3_2 <- function(cds, lineage, N, cores = 1, ID = TRUE){
  if(lineage != FALSE){
    sel.cells = cds@lineages[[lineage]]
  }
  sel.cells = sel.cells[sel.cells %in% colnames(cds)]
  cds_subset = cds[,sel.cells]
  family = stats::quasipoisson()
  model = "expression ~ splines::ns(pseudotime, df=3)"
  #names(cds_subset) <- rowData(cds_subset)$gene_short_name
  exp = as.data.frame(as.matrix(exprs(cds_subset)))
  exp = (t(exp)) /  (pData(cds_subset)[, 'Size_Factor'])
  exp = exp[,rownames(cds)]
  pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt <- pt[rownames(exp)]
  pt <- as.data.frame(pt)
  colnames(pt) <- c("pseudotime")
  pt$cell <- rownames(exp)
  UMAP <- reducedDims(cds_subset)[["UMAP"]]
  UMAP <- UMAP[rownames(exp),]
  if(ID == FALSE){
    exp = cbind(pt, UMAP, exp)
    exp = exp[order(exp$pseudotime),]
    pt = exp[,"pseudotime"]
    if(N >= nrow(exp)){
      return(FALSE)
    }
    n <- floor(nrow(exp)/(N))
    leftover <- nrow(exp) - n * N
    x <- UMAP[,'umap_1']
    y <- UMAP[,'umap_2']
    if (leftover <= 1/4*n){
      pt.comp = compress_3_2(pt, leftover, n)
      UMAP.comp.x = compress_3_2(x, leftover, n)
      UMAP.comp.y = compress_3_2(y, leftover, n)
      max.pt = max(pt.comp)
      print(paste0("Compressing lineage ", lineage))
      mat <- as.data.frame(exp[,5:ncol(exp)])
      exp.comp = pbsapply(mat, compress_3_2, leftover, n)
      exp_data <- cbind(pt.comp, UMAP.comp.x, UMAP.comp.y, exp.comp)
      exp_data <- as.data.frame(exp_data)
      exp_data$pt.comp <- as.numeric(exp_data$pt.comp)
      exp_data_ordered <- exp_data[order(exp_data$pt.comp), ]
      mat <- exp_data_ordered[,4:(ncol(exp_data_ordered))]
      d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
      print("Fitting curves")
      fit = pbsapply(mat, fit.m3_3, pt = d, max.pt = max(d), N = N, cl = cores)
      fit = apply(fit, 2, as.numeric)
      return(list("expression" = exp_data_ordered, "expectation" = fit, "pseudotime" = d))
      exp$expression[exp$expression < 0] <- 0
      exp$expectation[exp$expectation < 0] <- 0
      }else{
      k = N/50
      n_1 <- floor(leftover/(k))
      leftover_2 <- leftover - n_1 * k
      pt.comp <- compress_3_3(pt, leftover, leftover_2, n, n_1, k)
      UMAP.comp.x <- compress_3_3(x, leftover, leftover_2, n, n_1, k)
      UMAP.comp.y <- compress_3_3(y, leftover, leftover_2, n, n_1, k)
      max.pt = max(pt.comp)
      print(paste0("Compressing lineage ", lineage))
      mat <- as.data.frame(exp[,5:ncol(exp)])
      exp.comp = pbsapply(mat, compress_3_3, leftover, leftover_2, n, n_1, k)
      exp_data <- cbind(pt.comp, UMAP.comp.x, UMAP.comp.y, exp.comp)
      exp_data <- as.data.frame(exp_data)
      exp_data$pt.comp <- as.numeric(exp_data$pt.comp)
      exp_data_ordered <- exp_data[order(exp_data$pt.comp), ]
      mat <- exp_data_ordered[,4:(ncol(exp_data_ordered))]
      d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
      print("Fitting curves")
      fit = pbsapply(mat, fit.m3_3, pt = d, max.pt = max(d), N = N, cl = cores)
      fit = apply(fit, 2, as.numeric)
      return(list("expression" = exp_data_ordered, "expectation" = fit, "pseudotime" = d))
      exp$expression[exp$expression < 0] <- 0
      exp$expectation[exp$expectation < 0] <- 0
      }
  }else{
    ID <- as.data.frame(as.factor((cds_subset@colData[rownames(exp),]['sample'][,1])))
    colnames(ID) <- c("ID")
    age <- as.data.frame(cds_subset@colData[rownames(exp),]['age'])
    sex <- as.data.frame(cds_subset@colData[rownames(exp),]['sex'])
    region_broad <- as.data.frame(cds_subset@colData[rownames(exp),]['region_broad'])
    exp = cbind(pt, ID, age, sex, region_broad, UMAP, exp)
    #use sliding window to compress pseudotime with ID information
    #add the mean values of pt if length of the ID is smaller than n
    exp = exp[order(exp$ID, exp$pseudotime),]
    pt = exp[,"pseudotime"]
    length <- length(unique(exp$ID))
    pt.comp <- vector(mode = "list", length = length)
    ID.comp <- vector(mode = "list", length = length)
    age.comp <- vector(mode = "list", length = length)
    sex.comp <- vector(mode = "list", length = length)
    region_broad.comp <- vector(mode = "list", length = length)
    UMAP.comp.x <- vector(mode = "list", length = length)
    UMAP.comp.y <- vector(mode = "list", length = length)
    len <- c()
    unique <- unique(exp$ID)
    position <- match(unique, exp$ID)
    unique_age <- exp$age[position]
    unique_sex <- exp$sex[position]
    unique_region_broad <- exp$region_broad[position]
    ID <- exp$ID
    if (length > N)
    {
      return(FALSE)
    }
    n <- round(N/length)
    abs <- (N-n*length)
    if (abs < 0){
      n <- c(rep(n,length-abs(abs)), rep(n-1,abs(abs)))
      n <- sample(n, size = length, replace = FALSE)
    }
    if (abs > 0){
      n <- c(rep(n,length-abs(abs)), rep(n+1,abs(abs)))
      n <- sample(n, size = length, replace = FALSE)
    }
    for (i in 1:(length-1)){
      pos <- match(unique[i], ID)
      pos_1 <- match(unique[i+1], ID)
      n_1 <- n[i]
      if (length(pt[pos:(pos_1-1)]) <= n_1){
        pt.comp_1 = pt[pos:(pos_1-1)]
        UMAP.comp.x_1 = UMAP[,'UMAP_1'][pos:(pos_1-1)]
        UMAP.comp.y_1 = UMAP[,'UMAP_2'][pos:(pos_1-1)]
        ID.comp_1 <- as.character(rep(unique[i], times=length(pt.comp_1)))
        age.comp_1 <- as.character(rep(unique_age[i], times=length(pt.comp_1)))
        sex.comp_1 <- as.character(rep(unique_sex[i], times=length(pt.comp_1)))
        region_broad.comp_1 <- as.character(rep(unique_region_broad[i], times=length(pt.comp_1)))
      }
      else{
        window <- (pos_1-pos)/(n_1+1)
        #step <- (pos_1-pos-window)/n_1
        step <- window
        len <- c(len, (pos_1-pos))
        pt.comp_1 = SlidingWindow("mean", pt[pos:(pos_1-1)], window, step)
        UMAP.comp.x_1 = SlidingWindow("mean", UMAP[,'UMAP_1'][pos:(pos_1-1)], window, step)
        UMAP.comp.y_1 = SlidingWindow("mean", UMAP[,'UMAP_2'][pos:(pos_1-1)], window, step)
        ID.comp_1 <- as.character(rep(unique[i], times=length(pt.comp_1)))
        age.comp_1 <- as.character(rep(unique_age[i], times=length(pt.comp_1)))
        sex.comp_1 <- as.character(rep(unique_sex[i], times=length(pt.comp_1)))
        region_broad.comp_1 <- as.character(rep(unique_region_broad[i], times=length(pt.comp_1)))
      }
      pt.comp[[i]] <- pt.comp_1
      UMAP.comp.x[[i]] = UMAP.comp.x_1
      UMAP.comp.y[[i]] = UMAP.comp.y_1
      ID.comp[[i]] <- ID.comp_1
      age.comp[[i]] <- age.comp_1
      sex.comp[[i]] <- sex.comp_1
      region_broad.comp[[i]] <- region_broad.comp_1
      if (i == (length-1)){
        pos <- match(unique[i+1], ID)
        pos_1 <- length(ID)
        n_1 <- n[(length(n))]
        window <- (pos_1-pos+1)/(n_1+1)
        #step = ((pos_1-pos+1-window)/n_1)
        step <- window
        pt.comp_1 = SlidingWindow("mean", pt[pos:pos_1], window, step)
        UMAP.comp.x_1 = SlidingWindow("mean", UMAP[,'UMAP_1'][pos:(pos_1)], window, step)
        UMAP.comp.y_1 = SlidingWindow("mean", UMAP[,'UMAP_2'][pos:(pos_1)], window, step)
        ID.comp_1 <- as.character(rep(unique[i+1], times=length(pt.comp_1)))
        age.comp_1 <- as.character(rep(unique_age[i+1], times=length(pt.comp_1)))
        sex.comp_1 <- as.character(rep(unique_sex[i+1], times=length(pt.comp_1)))
        region_broad.comp_1 <- as.character(rep(unique_region_broad[i+1], times=length(pt.comp_1)))
        if (length(pt[pos:pos_1]) <= n_1){
          pt.comp_1 = pt[pos:pos_1]
          UMAP.comp.x_1 = UMAP[,'UMAP_1'][pos:(pos_1)]
          UMAP.comp.y_1 = UMAP[,'UMAP_2'][pos:(pos_1)]
          ID.comp_1 <- as.character(rep(unique[i+1], times=length(pt.comp_1)))
          age.comp_1 <- as.character(rep(unique_age[i+1], times=length(pt.comp_1)))
          sex.comp_1 <- as.character(rep(unique_sex[i+1], times=length(pt.comp_1)))
          region_broad.comp_1 <- as.character(rep(unique_region_broad[i+1], times=length(pt.comp_1)))
        }
        pt.comp[[i+1]] <- pt.comp_1
        UMAP.comp.x[[i+1]] = UMAP.comp.x_1
        UMAP.comp.y[[i+1]] = UMAP.comp.y_1
        ID.comp[[i+1]] <- ID.comp_1
        age.comp[[i+1]] <- age.comp_1
        sex.comp[[i+1]] <- sex.comp_1
        region_broad.comp[[i+1]] <- region_broad.comp_1
      }
    }
    eta <- N-sum(lengths(pt.comp))
    if(eta > 0){
      max <- which.max(len)
      pos <- match(unique[max], ID)
      pos_1 <- match(unique[max+1], ID)
      window_new <- (pos_1-pos)/(n[max]+eta+1)
      #step_new <- (pos_1-pos-window_new)/(n[max]+eta)
      step_new <- window_new
      pt.comp[[max]] <- SlidingWindow("mean", pt[pos:(pos_1-1)], window_new, step_new)
      UMAP.comp.x[[max]] <- SlidingWindow("mean", UMAP[,'UMAP_1'][pos:(pos_1-1)], window_new, step_new)
      UMAP.comp.y[[max]] <- SlidingWindow("mean", UMAP[,'UMAP_2'][pos:(pos_1-1)], window_new, step_new)
      ID.comp[[max]] <- as.character(rep(unique[max], times=length(pt.comp[[max]])))
      age.comp[[max]] <- as.character(rep(unique_age[max], times=length(pt.comp[[max]])))
      sex.comp[[max]] <- as.character(rep(unique_sex[max], times=length(pt.comp[[max]])))
      region_broad.comp[[max]] <- as.character(rep(unique_region_broad[max], times=length(pt.comp[[max]])))
    }
    pt.comp <- unlist(pt.comp)
    UMAP.comp.x <- unlist(UMAP.comp.x)
    UMAP.comp.y <- unlist(UMAP.comp.y)
    ID.comp <- unlist(ID.comp)
    age.comp <- unlist(age.comp)
    sex.comp <- unlist(sex.comp)
    region_broad.comp <- unlist(region_broad.comp)
    max.pt = max(pt.comp)
    
    #use sliding window to compress expression with ID information
    len <- c()
    exp.comp <- vector(mode = "list", length = length)
    print(paste0("Compressing lineage ", lineage))
    mat <- as.data.frame(exp[,9:ncol(exp)])
    exp.comp <- vector(mode = "list", length = length)
    exp.comp = pbsapply(mat, compress_3, length = length, unique = unique, ID = ID, n = n, N = N, l = 108, cl = cores)
    exp_data <- cbind(pt.comp, ID.comp, age.comp, sex.comp, region_broad.comp, UMAP.comp.x, UMAP.comp.y, exp.comp)
    exp_data <- as.data.frame(exp_data)
    exp_data$pt.comp <- as.numeric(exp_data$pt.comp)
    exp_data_ordered <- exp_data[order(exp_data$pt.comp), ]
    mat <- exp_data_ordered[,8:(ncol(exp_data_ordered))]
    d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    print("Fitting curves")
    fit = pbsapply(mat, fit.m3_3, pt = d, max.pt = max(d), N = N, cl = cores)
    fit = apply(fit, 2, as.numeric)
    return(list("expression" = exp_data_ordered, "expectation" = fit, "pseudotime" = d))
    exp$expression[exp$expression < 0] <- 0
    exp$expectation[exp$expectation < 0] <- 0
  }
  return(exp)
}

compress_3_3 <- function(df, leftover, leftover_2, n, n_1, k){
  df.comp_1 = sw(df[1:((length(df) - (n*k + leftover)))], n, mean, n)
  df.comp_2 = sw(df[((length(df) - (n*k + leftover))+1):(length(df) - (n + n_1 + leftover_2))], (n+n_1), mean, (n+n_1))
  df.comp_3 = mean(df[((length(df) - (n + n_1 + leftover_2))+1):(length(df))])
  df.comp <- c(df.comp_1, df.comp_2, df.comp_3)
  return(df.comp)
}

sw <- function(vec, window_size, FUN, step = 1) {
  n <- length(vec)
  starts <- seq(1, n - window_size + 1, by = step)
  sapply(starts, function(i) FUN(vec[i:(i + window_size - 1)]))
}

compress_lineage_v3 <- function(cds, lineage, N, cores = 1){
  cds_name = deparse(substitute(cds))
  input = paste0("compress_expression_v3(",cds_name,", lineage = '", lineage, "', N = ", N, ", cores = ", cores, ")")
  exp = eval(parse(text=input))
  input = paste0(cds_name, "@expression$", lineage, " <- exp$expression")
  eval(parse(text=input))
  input = paste0(cds_name, "@expectation$", lineage, " <- exp$expectation")
  eval(parse(text=input))
  input = paste0(cds_name, "@pseudotime$", lineage, " <- exp$pseudotime")
  eval(parse(text=input))
  eval(parse(text=paste0("return(",cds_name, ")")))
}

compress_expression_v3 <- function(cds, lineage, N, cores = 1){
  cds_name = deparse(substitute(cds))
  if(lineage != FALSE){
    input = paste0("sel.cells = ",cds_name,"@lineages$", lineage)
    eval(parse(text=input))
  }
  sel.cells = sel.cells[sel.cells %in% colnames(cds)]
  cds_subset = cds[,sel.cells]
  family = stats::quasipoisson()
  model = "expression ~ splines::ns(pseudotime, df=3)"
  #names(cds_subset) <- rowData(cds_subset)$gene_short_name
  exp = as.data.frame(as.matrix(exprs(cds_subset)))
  exp = (t(exp)) /  (pData(cds_subset)[, 'Size_Factor'])
  pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt <- pt[rownames(exp)]
  pt <- as.data.frame(pt)
  colnames(pt) <- c("pseudotime")
  pt$cell <- rownames(exp)
  ID <- as.data.frame(as.factor((cds_subset@colData[rownames(exp),]['sample'][,1])))
  colnames(ID) <- c("ID")
  age <- as.data.frame(cds_subset@colData[rownames(exp),]['age'])
  sex <- as.data.frame(cds_subset@colData[rownames(exp),]['sex'])
  region_broad <- as.data.frame(cds_subset@colData[rownames(exp),]['region_broad'])
  UMAP <- reducedDims(cds_subset)[["UMAP"]]
  UMAP <- UMAP[rownames(exp),]
  exp = cbind(pt, ID, age, sex, region_broad, UMAP, exp)
  #use sliding window to compress pseudotime with ID information
  #add the mean values of pt if length of the ID is smaller than n
  exp = exp[order(exp$ID, exp$pseudotime),]
  pt = exp[,"pseudotime"]
  length <- length(unique(exp$ID))
  pt.comp <- vector(mode = "list", length = length)
  ID.comp <- vector(mode = "list", length = length)
  age.comp <- vector(mode = "list", length = length)
  sex.comp <- vector(mode = "list", length = length)
  region_broad.comp <- vector(mode = "list", length = length)
  UMAP.comp.x <- vector(mode = "list", length = length)
  UMAP.comp.y <- vector(mode = "list", length = length)
  len <- c()
  unique <- unique(exp$ID)
  position <- match(unique, exp$ID)
  unique_age <- exp$age[position]
  unique_sex <- exp$sex[position]
  unique_region_broad <- exp$region_broad[position]
  ID <- exp$ID
  if (length > N)
  {
    return(FALSE)
  }
  n <- round(N/length)
  abs <- (N-n*length)
  if (abs < 0){
    n <- c(rep(n,length-abs(abs)), rep(n-1,abs(abs)))
    n <- sample(n, size = length, replace = FALSE)
  }
  if (abs > 0){
    n <- c(rep(n,length-abs(abs)), rep(n+1,abs(abs)))
    n <- sample(n, size = length, replace = FALSE)
  }
  for (i in 1:(length-1)){
    pos <- match(unique[i], ID)
    pos_1 <- match(unique[i+1], ID)
    n_1 <- n[i]
    if (length(pt[pos:(pos_1-1)]) <= n_1){
      pt.comp_1 = pt[pos:(pos_1-1)]
      UMAP.comp.x_1 = UMAP[,'UMAP_1'][pos:(pos_1-1)]
      UMAP.comp.y_1 = UMAP[,'UMAP_2'][pos:(pos_1-1)]
      ID.comp_1 <- as.character(rep(unique[i], times=length(pt.comp_1)))
      age.comp_1 <- as.character(rep(unique_age[i], times=length(pt.comp_1)))
      sex.comp_1 <- as.character(rep(unique_sex[i], times=length(pt.comp_1)))
      region_broad.comp_1 <- as.character(rep(unique_region_broad[i], times=length(pt.comp_1)))
    }
    else{
      window <- (pos_1-pos)/(n_1+1)
      #step <- (pos_1-pos-window)/n_1
      step <- window
      len <- c(len, (pos_1-pos))
      pt.comp_1 = SlidingWindow("mean", pt[pos:(pos_1-1)], window, step)
      UMAP.comp.x_1 = SlidingWindow("mean", UMAP[,'UMAP_1'][pos:(pos_1-1)], window, step)
      UMAP.comp.y_1 = SlidingWindow("mean", UMAP[,'UMAP_2'][pos:(pos_1-1)], window, step)
      ID.comp_1 <- as.character(rep(unique[i], times=length(pt.comp_1)))
      age.comp_1 <- as.character(rep(unique_age[i], times=length(pt.comp_1)))
      sex.comp_1 <- as.character(rep(unique_sex[i], times=length(pt.comp_1)))
      region_broad.comp_1 <- as.character(rep(unique_region_broad[i], times=length(pt.comp_1)))
    }
    pt.comp[[i]] <- pt.comp_1
    UMAP.comp.x[[i]] = UMAP.comp.x_1
    UMAP.comp.y[[i]] = UMAP.comp.y_1
    ID.comp[[i]] <- ID.comp_1
    age.comp[[i]] <- age.comp_1
    sex.comp[[i]] <- sex.comp_1
    region_broad.comp[[i]] <- region_broad.comp_1
    if (i == (length-1)){
      pos <- match(unique[i+1], ID)
      pos_1 <- length(ID)
      n_1 <- n[(length(n))]
      window <- (pos_1-pos+1)/(n_1+1)
      #step = ((pos_1-pos+1-window)/n_1)
      step <- window
      pt.comp_1 = SlidingWindow("mean", pt[pos:pos_1], window, step)
      UMAP.comp.x_1 = SlidingWindow("mean", UMAP[,'UMAP_1'][pos:(pos_1)], window, step)
      UMAP.comp.y_1 = SlidingWindow("mean", UMAP[,'UMAP_2'][pos:(pos_1)], window, step)
      ID.comp_1 <- as.character(rep(unique[i+1], times=length(pt.comp_1)))
      age.comp_1 <- as.character(rep(unique_age[i+1], times=length(pt.comp_1)))
      sex.comp_1 <- as.character(rep(unique_sex[i+1], times=length(pt.comp_1)))
      region_broad.comp_1 <- as.character(rep(unique_region_broad[i+1], times=length(pt.comp_1)))
      if (length(pt[pos:pos_1]) <= n_1){
        pt.comp_1 = pt[pos:pos_1]
        UMAP.comp.x_1 = UMAP[,'UMAP_1'][pos:(pos_1)]
        UMAP.comp.y_1 = UMAP[,'UMAP_2'][pos:(pos_1)]
        ID.comp_1 <- as.character(rep(unique[i+1], times=length(pt.comp_1)))
        age.comp_1 <- as.character(rep(unique_age[i+1], times=length(pt.comp_1)))
        sex.comp_1 <- as.character(rep(unique_sex[i+1], times=length(pt.comp_1)))
        region_broad.comp_1 <- as.character(rep(unique_region_broad[i+1], times=length(pt.comp_1)))
      }
      pt.comp[[i+1]] <- pt.comp_1
      UMAP.comp.x[[i+1]] = UMAP.comp.x_1
      UMAP.comp.y[[i+1]] = UMAP.comp.y_1
      ID.comp[[i+1]] <- ID.comp_1
      age.comp[[i+1]] <- age.comp_1
      sex.comp[[i+1]] <- sex.comp_1
      region_broad.comp[[i+1]] <- region_broad.comp_1
    }
  }
  eta <- N-sum(lengths(pt.comp))
  if(eta > 0){
    max <- which.max(len)
    pos <- match(unique[max], ID)
    pos_1 <- match(unique[max+1], ID)
    window_new <- (pos_1-pos)/(n[max]+eta+1)
    #step_new <- (pos_1-pos-window_new)/(n[max]+eta)
    step_new <- window_new
    pt.comp[[max]] <- SlidingWindow("mean", pt[pos:(pos_1-1)], window_new, step_new)
    UMAP.comp.x[[max]] <- SlidingWindow("mean", UMAP[,'UMAP_1'][pos:(pos_1-1)], window_new, step_new)
    UMAP.comp.y[[max]] <- SlidingWindow("mean", UMAP[,'UMAP_2'][pos:(pos_1-1)], window_new, step_new)
    ID.comp[[max]] <- as.character(rep(unique[max], times=length(pt.comp[[max]])))
    age.comp[[max]] <- as.character(rep(unique_age[max], times=length(pt.comp[[max]])))
    sex.comp[[max]] <- as.character(rep(unique_sex[max], times=length(pt.comp[[max]])))
    region_broad.comp[[max]] <- as.character(rep(unique_region_broad[max], times=length(pt.comp[[max]])))
  }
  pt.comp <- unlist(pt.comp)
  UMAP.comp.x <- unlist(UMAP.comp.x)
  UMAP.comp.y <- unlist(UMAP.comp.y)
  ID.comp <- unlist(ID.comp)
  age.comp <- unlist(age.comp)
  sex.comp <- unlist(sex.comp)
  region_broad.comp <- unlist(region_broad.comp)
  max.pt = max(pt.comp)
  
  #use sliding window to compress expression with ID information
  len <- c()
  exp.comp <- vector(mode = "list", length = length)
  print(paste0("Compressing lineage ", lineage))
  mat <- as.data.frame(exp[,9:ncol(exp)])
  exp.comp <- vector(mode = "list", length = length)
  exp.comp = pbsapply(mat, compress_3, length = length, unique = unique, ID = ID, n = n, N = N, l = 108, cl = cores)
  exp_data <- cbind(pt.comp, ID.comp, age.comp, sex.comp, region_broad.comp, UMAP.comp.x, UMAP.comp.y, exp.comp)
  exp_data <- as.data.frame(exp_data)
  exp_data$pt.comp <- as.numeric(exp_data$pt.comp)
  exp_data_ordered <- exp_data[order(exp_data$pt.comp), ]
  mat <- exp_data_ordered[,8:(ncol(exp_data_ordered))]
  d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
  print("Fitting curves")
  fit = pbsapply(mat, fit.m3_3, pt = d, max.pt = max(d), N = N, cl = cores)
  fit = apply(fit, 2, as.numeric)
  return(list("expression" = exp_data_ordered, "expectation" = fit, "pseudotime" = d))
  exp$expression[exp$expression < 0] <- 0
  exp$expectation[exp$expectation < 0] <- 0
  exp
}

compress_3 <- function(df, length, unique, ID, n, N, l){
  len <- c()
  exp.comp <- vector(mode = "list", length = l)
  for (i in 1:(length-1)){
    pos <- match(unique[i], ID)
    pos_1 <- match(unique[i+1], ID)
    n_1 <- n[i]
    window <- (pos_1-pos)/n_1
    step <- (pos_1-pos-window)/n_1
    len <- c(len, (pos_1-pos))
    if (length(df[pos:(pos_1-1)]) <= n_1){
      exp.comp_1 = df[pos:(pos_1-1)]
    }
    else{
      exp.comp_1 = SlidingWindow("mean", df[pos:(pos_1-1)], window, step)
    }
    exp.comp[[i]] <- exp.comp_1
    if (i == (length-1)){
      pos <- match(unique[i+1], ID)
      pos_1 <- length(ID)
      n_1 <- n[(length(n))]
      window <- (pos_1-pos+1)/n_1
      step = ((pos_1-pos+1-window)/n_1)
      len <- c(len, (pos_1-pos))
      if (length(df[pos:pos_1]) <= n_1){
        exp.comp_1 = df[pos:pos_1]
      }
      else{
        exp.comp_1 = SlidingWindow("mean", df[pos:pos_1], window, step)
      }
      exp.comp[[i+1]] <- exp.comp_1
    }
  }
  eta <- N-sum(lengths(exp.comp))
  if(eta > 0){
    max <- which.max(len)
    pos <- match(unique[max], ID)
    pos_1 <- match(unique[max+1], ID)
    window_new <- (pos_1-pos)/(n[max]+eta)
    step_new <- (pos_1-pos-window_new)/(n[max]+eta)
    exp.comp[[max]] <- SlidingWindow("mean", df[pos:(pos_1-1)], window_new, step_new)
  }
  exp.comp <- unlist(exp.comp)
}
