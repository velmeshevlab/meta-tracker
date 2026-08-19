# Lineage expression compression into meta-cells + smoothed expectation curves.

# Cross-platform parallel backend: MulticoreParam (fork) on Unix, SnowParam
# (PSOCK) on Windows, SerialParam when cores <= 1.
.compress_bpparam <- function(cores) {
  if (is.null(cores) || cores <= 1) return(BiocParallel::SerialParam())
  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = cores, progressbar = TRUE)
  } else {
    BiocParallel::MulticoreParam(workers = cores, progressbar = TRUE)
  }
}

# Quasipoisson spline fit for one gene's meta-cell counts, predicted on a
# regular pseudotime grid. Returns a length-N numeric vector (NA on failure).
.fit_m3 <- function(exp.sel, pt, size_factor, predict_pt, lineage, model, N) {
  if (!requireNamespace("speedglm", quietly = TRUE)) stop("speedglm not installed")
  exp_data.sel <- data.frame(
    pseudotime  = as.numeric(pt),
    size_factor = as.numeric(size_factor),
    expression  = as.numeric(exp.sel))
  tryCatch({
    fit_model <- speedglm::speedglm(model, data = exp_data.sel,
                                    family = quasipoisson(), acc = 1e-3,
                                    model = FALSE, y = FALSE)
    newdata <- data.frame(pseudotime = as.numeric(predict_pt), size_factor = 1)
    predict(fit_model, newdata = newdata, type = "response")
  }, error = function(e) rep(NA_real_, N))
}

# Compress a lineage's cells into N pseudotime-ordered meta-cells and fit a
# smoothed expectation curve per gene.
.compress_expression <- function(cds, lineage, N, cores = 1, method = "sum", ID = FALSE){
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
  colnames(UMAP) <- c("umap_1", "umap_2")
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
    # Cross-platform parallel fit. Capture the fit function and inputs as objects
    # so PSOCK (Windows) workers don't need the package namespace loaded.
    fitfun <- .fit_m3
    FUN <- function(i) fitfun(exp.sel = mat_m[, i], pt = d, size_factor = size_factor,
                              predict_pt = predict_pt, lineage = lineage, model = model, N = N)
    res <- BiocParallel::bplapply(seq_along(genes), FUN, BPPARAM = .compress_bpparam(cores))
    fit_list_2 <- do.call(cbind, res)
    colnames(fit_list_2) <- genes
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

#' Compress a lineage into meta-cells with smoothed expectation curves
#'
#' Bins a lineage's cells into \code{N} pseudotime-ordered meta-cells, sums (or
#' also means) their expression, and fits a quasipoisson spline per gene to give
#' a smoothed expectation. Gene fits run in parallel across cores on Windows,
#' macOS, and Linux.
#'
#' @param cds A \code{metatracker_data_set} with the lineage isolated.
#' @param lineage Lineage name.
#' @param N Number of meta-cells (and prediction grid points).
#' @param method "sum" (default) or any other value to also compute the mean matrix.
#' @param cores Worker processes for the per-gene fits (default 1 = serial).
#' @param ID Passed through to the compression routine (default FALSE).
#' @return The \code{cds} with \code{@lineages}, \code{@expression},
#'   \code{@expectation}, and \code{@pseudotime} populated for \code{lineage}.
#' @export
compress_lineage <- function(cds, lineage, N, method = "sum", cores = 1, ID = FALSE){
  exp = .compress_expression(cds, lineage = lineage, method = method, N = N, cores = cores, ID = ID)
  cds@lineages[[lineage]]    <- exp$lineage
  cds@expression[[lineage]]  <- exp$expression
  cds@expectation[[lineage]] <- exp$expectation
  cds@pseudotime[[lineage]]  <- exp$pseudotime
  cds
}
