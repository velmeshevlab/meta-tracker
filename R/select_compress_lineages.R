#extend monocle3 class to add additional slots
setClass("cell_data_set_ext", contains = "cell_data_set", slots=c(graphs = "list", lineages="list", expression="list", expectation="list", pseudotime="list")) -> cell_data_set_ext

#' @export
#import monocle object
import_monocle <-function(cds){
cds <- as(cds,"cell_data_set_ext")
return(cds)
}

#' @export
#generate node plot
#filter = T will display only the nodes at branch points and at the ends of trajectories
#N controls the density of nodes to display if filter = T (larger values = less dense, N = 1 displays all nodes)
node_plot <- function(cds, filter = F, N = 50){
Y <- cds@principal_graph_aux[["UMAP"]]$dp_mst
d = as.data.frame(t(Y))
if(filter == T){
g = principal_graph(cds)[["UMAP"]]
dd = degree(g)
names1 = names(dd[dd > 2 | dd == 1])
names2 = names(dd[dd == 2])
names2 = sample(names2, length(names2)/N, replace = F)
d.f = d[c(names1, names2),]
ggplot(data=d, aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.01) + geom_text_repel(data=d.f, aes(x=UMAP_1, y=UMAP_2), label=rownames(d.f), size=0.3, hjust = 2, color = "red", max.overlaps = Inf, segment.size = 0.1) + monocle_theme_opts()
}
else{
ggplot(data=d, aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.01) + geom_text_repel(data=d, aes(x=UMAP_1, y=UMAP_2), label=rownames(d), size=0.3, hjust = 2, color = "red", max.overlaps = Inf, segment.size = 0.1) + monocle_theme_opts()
}
}

compress_lineage_v3 <- function(cds, lineage, start, gene = FALSE, N = 500, cores = F){
  cds_name = deparse(substitute(cds))
  if(gene == FALSE){
    input = paste0("compress_expression_v3(",cds_name,", lineage = '", lineage, "', start = ", start, ", gene = ", gene, ", N = ", N, ", cores = ", cores, ")")
  }
  else{
    input = paste0("compress_expression_v3(",cds_name,", lineage = '", lineage, "', start = ", start, ", gene = '", gene, "', N = ", N, ", cores = ", cores, ")")
  }
  exp = eval(parse(text=input))
  input = paste0(cds_name, "@expression$", lineage, " <- exp$expression")
  eval(parse(text=input))
  input = paste0(cds_name, "@expectation$", lineage, " <- exp$expectation")
  eval(parse(text=input))
  input = paste0(cds_name, "@pseudotime$", lineage, " <- exp$pseudotime")
  eval(parse(text=input))
  eval(parse(text=paste0("return(",cds_name, ")")))
}

compress_expression_v3 <- function(cds, lineage, start, gene = FALSE, N = 500, cores = F){
  cds_name = deparse(substitute(cds))
  if(cores != F){
    cl <- makeCluster(cores)
    clusterEvalQ(cl, c(library(evobiR)))
  }
  input = paste0("get_lineage_object(",cds_name,", lineage = '", lineage, "', start = ", start, ")")
  cds_subset = eval(parse(text=input))
  family = stats::quasipoisson()
  model = "expression ~ splines::ns(pseudotime, df=3)"
  names(cds_subset) <- rowData(cds_subset)$gene_short_name
  exp = as.data.frame(as_matrix(exprs(cds_subset)))
  exp = t(t(exp) /  pData(cds_subset)[, 'Size_Factor'])
  pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt <- as.data.frame(pt)
  colnames(pt) <- c("pseudotime")
  pt$cell <- rownames(pt)
  rownames <- pt$cell
  pt$ID <- cds_subset@colData[rownames,]$sample
  exp = cbind(pt, t(exp))
  #use sliding window to compress pseudotime with ID information
  #add the mean values of pt if length of the ID is smaller than n
  exp = exp[order(exp$ID, exp$pseudotime),]
  pt = exp[,"pseudotime"]
  length <- length(unique(exp$ID))
  UMAP <- as.data.frame(subset2@reductions$umap@cell.embeddings)
  UMAP <- UMAP[rownames(exp),]
  pt.comp <- vector(mode = "list", length = length)
  ID.comp <- vector(mode = "list", length = length)
  UMAP.comp.x <- vector(mode = "list", length = length)
  UMAP.comp.y <- vector(mode = "list", length = length)
  len <- c()
  unique <- unique(exp$ID)
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
    window <- (pos_1-pos)/(n_1+1)
    #step <- (pos_1-pos-window)/n_1
    step <- window
    len <- c(len, (pos_1-pos))
    pt.comp_1 = SlidingWindow("mean", pt[pos:(pos_1-1)], window, step)
    UMAP.comp.x_1 = SlidingWindow("mean", UMAP$UMAP_1[pos:(pos_1-1)], window, step)
    UMAP.comp.y_1 = SlidingWindow("mean", UMAP$UMAP_2[pos:(pos_1-1)], window, step)
    ID.comp_1 <- as.character(rep(unique[i], times=length(pt.comp_1)))
    if (length(pt[pos:(pos_1-1)]) <= n_1){
      pt.comp_1 = pt[pos:(pos_1-1)]
      UMAP.comp.x_1 = UMAP$UMAP_1[pos:(pos_1-1)]
      UMAP.comp.y_1 = UMAP$UMAP_2[pos:(pos_1-1)]
      ID.comp_1 <- as.character(rep(unique[i], times=length(pt.comp_1)))
    }
    pt.comp[[i]] <- pt.comp_1
    UMAP.comp.x[[i]] = UMAP.comp.x_1
    UMAP.comp.y[[i]] = UMAP.comp.y_1
    ID.comp[[i]] <- ID.comp_1
    if (i == (length-1)){
      pos <- match(unique[i+1], ID)
      pos_1 <- length(ID)
      n_1 <- n[(length(n))]
      window <- (pos_1-pos+1)/(n_1+1)
      #step = ((pos_1-pos+1-window)/n_1)
      step <- window
      pt.comp_1 = SlidingWindow("mean", pt[pos:pos_1], window, step)
      UMAP.comp.x_1 = SlidingWindow("mean", UMAP$UMAP_1[pos:(pos_1)], window, step)
      UMAP.comp.y_1 = SlidingWindow("mean", UMAP$UMAP_2[pos:(pos_1)], window, step)
      ID.comp_1 <- as.character(rep(unique[i+1], times=length(pt.comp_1)))
      if (length(pt[pos:pos_1]) <= n_1){
        pt.comp_1 = pt[pos:pos_1]
        UMAP.comp.x_1 = UMAP$UMAP_1[pos:(pos_1)]
        UMAP.comp.y_1 = UMAP$UMAP_2[pos:(pos_1)]
        ID.comp_1 <- as.character(rep(unique[i+1], times=length(pt.comp_1)))
      }
      pt.comp[[i+1]] <- pt.comp_1
      UMAP.comp.x[[i+1]] = UMAP.comp.x_1
      UMAP.comp.y[[i+1]] = UMAP.comp.y_1
      ID.comp[[i+1]] <- ID.comp_1
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
    UMAP.comp.x[[max]] <- SlidingWindow("mean", UMAP$UMAP_1[pos:(pos_1-1)], window_new, step_new)
    UMAP.comp.y[[max]] <- SlidingWindow("mean", UMAP$UMAP_2[pos:(pos_1-1)], window_new, step_new)
    ID.comp[[max]] <- as.character(rep(unique[max], times=length(pt.comp[[max]])))
  }
  pt.comp <- unlist(pt.comp)
  UMAP.comp.x <- unlist(UMAP.comp.x)
  UMAP.comp.y <- unlist(UMAP.comp.y)
  ID.comp <- unlist(ID.comp)
  max.pt = max(pt.comp)
  
  #use sliding window to compress expression with ID information
  len <- c()
  exp.comp <- vector(mode = "list", length = length)
  if(gene != F){
    exp_test = exp[,c("pseudotime", "ID", gene)]
    for (i in 1:(length-1)){
      pos <- match(unique[i], ID)
      pos_1 <- match(unique[i+1], ID)
      n_1 <- n[i]
      window <- (pos_1-pos)/(n_1+1)
      #step <- (pos_1-pos-window)/n_1
      step <- window
      len <- c(len, (pos_1-pos))
      exp.comp_1 = SlidingWindow("mean", exp_test$gene[pos:(pos_1-1)], window, step)
      if (length(exp_test$gene[pos:(pos_1-1)]) <= n_1){
        exp.comp_1 = exp_test$gene[pos:(pos_1-1)]
      }
      exp.comp[[i]] <- exp.comp_1
      if (i == (length-1)){
        pos <- match(unique[i+1], ID)
        pos_1 <- length(ID)
        n_1 <- n[(length(n))]
        window <- (pos_1-pos+1)/(n_1+1)
        #step = ((pos_1-pos+1-window)/n_1)
        step <- window
        exp.comp_1 = SlidingWindow("mean", exp_test$gene[pos:pos_1], window, step)
        if (length(exp_test$gene[pos:pos_1]) <= n_1){
          exp.comp_1 = exp_test$gene[pos:pos_1]
        }
        exp.comp[[i+1]] <- exp.comp_1
      }
    }
    eta <- N-sum(lengths(exp.comp))
    if(eta > 0){
      max <- which.max(len)
      pos <- match(unique[max], ID)
      pos_1 <- match(unique[max+1], ID)
      window_new <- (pos_1-pos)/(n[max]+eta+1)
      #step_new <- (pos_1-pos-window_new)/(n[max]+eta)
      step_new <- window_new
      exp.comp[[max]] <- SlidingWindow("mean", exp_test$gene[pos:(pos_1-1)], window_new, step_new)
    }
    exp.comp <- unlist(exp.comp)
    
  }
  else{
    print(paste0("Compressing lineage ", lineage, " and fitting curves"))
    mat <- as.data.frame(exp[,4:ncol(exp)])
    exp.comp = pbapply(mat, 2, compress_3, length = length, unique = unique, ID = ID, n = n, N = 500, l = 108, cl = cl)
  }
  
  if(gene != F){
    exp_data.sel = cbind(pt.comp, exp.comp)
    exp_data.sel = as.data.frame(exp_data.sel)
    colnames(exp_data.sel) <- c("pseudotime", "expression")
    exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
    exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
    d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    colnames(d) <- c("pseudotime")
    tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
    fit = predict(fit, newdata = d, type='response')
    }, error=function(cond) {fit = as.data.frame(rep(0, N))})
    exp = as.data.frame(cbind(exp.comp, fit, d))
    colnames(exp) <- c("expression", "expectation", "pseudotime")
    exp$expression[exp$expression < 0] <- 0
    exp$expectation[exp$expectation < 0] <- 0
  }
  else{
    exp_data <- cbind(pt.comp, ID.comp, UMAP.comp.x, UMAP.comp.y, exp.comp)
    exp_data <- as.data.frame(exp_data)
    exp_data$pt.comp <- as.numeric(exp_data$pt.comp)
    exp_data_ordered <- exp_data[order(exp_data$pt.comp), ]
    mat <- exp_data_ordered[,5:(ncol(exp_data_ordered))]
    d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    if(cores != F){
      fit = pbsapply(mat, fit.m3, pt = d, max.pt = max(d), N = N, cl = cl)
    }
    else{
      fit = pbsapply(mat, fit.m3, pt = d, max.pt = max(d), N = N)
    }
    fit = apply(fit, 2, as.numeric)
    return(list("expression" = exp_data_ordered, "expectation" = fit, "pseudotime" = d))
  }
  exp$expression[exp$expression < 0] <- 0
  exp$expectation[exp$expectation < 0] <- 0
  if(cores != F){
    stopCluster(cl)
  }
  return(exp)
}

compress_3 <- function(df, length, unique, ID, n, N=500, l){
  len <- c()
  exp.comp <- vector(mode = "list", length = l)
  for (i in 1:(length-1)){
    pos <- match(unique[i], ID)
    pos_1 <- match(unique[i+1], ID)
    n_1 <- n[i]
    window <- (pos_1-pos)/n_1
    step <- (pos_1-pos-window)/n_1
    len <- c(len, (pos_1-pos))
    exp.comp_1 = SlidingWindow("mean", df[pos:(pos_1-1)], window, step)
    if (length(df[pos:(pos_1-1)]) <= n_1){
      exp.comp_1 = df[pos:(pos_1-1)]
    }
    exp.comp[[i]] <- exp.comp_1
    if (i == (length-1)){
      pos <- match(unique[i+1], ID)
      pos_1 <- length(ID)
      n_1 <- n[(length(n))]
      window <- (pos_1-pos+1)/n_1
      step = ((pos_1-pos+1-window)/n_1)
      len <- c(len, (pos_1-pos))
      exp.comp_1 = SlidingWindow("mean", df[pos:pos_1], window, step)
      if (length(df[pos:pos_1]) <= n_1){
        exp.comp_1 = df[pos:pos_1]
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

fit.m3 <- function(exp.sel, pt, max.pt, model = "expression ~ splines::ns(pseudotime, df=3)", N = 500){
  require(speedglm)
  family = stats::quasipoisson()
  exp_data.sel = cbind(pt, exp.sel)
  colnames(exp_data.sel) <- c("pseudotime","expression")
  exp_data.sel = as.data.frame(exp_data.sel)
  exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
  exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
  tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
  d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
  colnames(d) <- c("pseudotime")
  fit = stats::predict(fit, newdata=d, type="response")
  return(fit)
  }, error=function(cond) {return(rep("NA", N))})
}
cds_name = deparse(substitute(cds_subset))
input = paste0(cds_name, "@expression$", lineage, " <- exp_c$expression")
eval(parse(text=input))
input = paste0(cds_name, "@expectation$", lineage, " <- exp_c$expectation")
eval(parse(text=input))
input = paste0(cds_name, "@pseudotime$", lineage, " <- exp_c$pseudotime")
eval(parse(text=input))

