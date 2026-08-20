# Lineage expression compression into meta-cells + smoothed expectation curves.

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
.compress_expression <- function(cds, lineage, N, method = "sum", ID = FALSE, progress = TRUE){
  #Extract the lineage object
  cds_sub <- get_lineage_object(cds, lineage)
  #Get pseudotime from the principal graph aux slot
  updated_pt <- cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  #Store back in cds@lineages as a list (normalise first so re-compression of an
  #already-compressed lineage doesn't nest a list inside `name`).
  cell_barcodes <- .lineage_cells(cds@lineages[[lineage]])
  cds@lineages[[lineage]] <- list(
    name = cell_barcodes,
    updated_pt = updated_pt)
  sel.cells = cell_barcodes
  sel.cells = sel.cells[sel.cells %in% colnames(cds)]
  if (length(sel.cells) == 0)
    stop("Lineage '", lineage, "' has no cells present in cds.", call. = FALSE)
  cds_subset = cds[, sel.cells]
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
    predict_pt <- seq(0, 1, length.out = N)
    mat_m <- as.matrix(mat)
    genes <- colnames(mat_m)
    # Per-gene fit. Progress bar shown for single-lineage use; suppressed when
    # called from compress_lineages() (the outer per-lineage bar is shown there,
    # and forked children would otherwise garble the console).
    .fitcol <- function(i) {
      .fit_m3(exp.sel = mat_m[, i], pt = d, size_factor = size_factor,
              predict_pt = predict_pt, lineage = lineage, model = model, N = N)
    }
    fit_list_2 <- if (isTRUE(progress)) {
      op <- pbapply::pboptions(type = "timer")   # force the bar on for this fit
      on.exit(pbapply::pboptions(op), add = TRUE)
      pbapply::pbsapply(seq_along(genes), .fitcol)
    } else {
      sapply(seq_along(genes), .fitcol)
    }
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
#' a smoothed expectation, with a progress bar over genes.
#'
#' @param cds A \code{metatracker_data_set} with the lineage isolated.
#' @param lineage Lineage name.
#' @param N Number of meta-cells (and prediction grid points).
#' @param method "sum" (default) or any other value to also compute the mean matrix.
#' @param ID Passed through to the compression routine (default FALSE).
#' @param progress Show the per-gene progress bar (default TRUE). Set FALSE when
#'   fitting many lineages in parallel (see \code{compress_lineages}).
#' @return The \code{cds} with \code{@lineages}, \code{@expression},
#'   \code{@expectation}, and \code{@pseudotime} populated for \code{lineage}.
#' @export
compress_lineage <- function(cds, lineage, N, method = "sum", ID = FALSE, progress = TRUE){
  exp = .compress_expression(cds, lineage = lineage, method = method, N = N,
                             ID = ID, progress = progress)
  cds@lineages[[lineage]]    <- exp$lineage
  cds@expression[[lineage]]  <- exp$expression
  cds@expectation[[lineage]] <- exp$expectation
  cds@pseudotime[[lineage]]  <- exp$pseudotime
  cds
}

#' Compress every lineage, optionally in parallel
#'
#' Runs \code{compress_lineage} across several lineages and merges the results
#' into one object. Parallelised over lineages with \code{pbapply::pblapply},
#' using the same \code{cl} convention as \code{isolate_lineage}: an integer
#' forks on Unix/macOS and runs serially on Windows, or pass a cluster object to
#' parallelise anywhere.
#'
#' Because lineages are usually dispatched all at once, a per-lineage bar would
#' jump straight from 0 to 100\%. Instead, the first lineage streams its own
#' per-gene progress bar as a representative indicator of how the run is going;
#' the other lineages fit silently. (Visible under forking; with a PSOCK cluster
#' worker output is not forwarded to the console.)
#'
#' @param cds A \code{metatracker_data_set}.
#' @param lineages Lineage names to compress (default: all in \code{cds@lineages}).
#' @param N Number of meta-cells per lineage.
#' @param method "sum" (default) or any other value to also compute the mean matrix.
#' @param ID Passed through to \code{compress_lineage} (default FALSE).
#' @param cl Passed to \code{pbapply::pblapply}: an integer worker count (fork on
#'   Unix, serial on Windows) or a \code{parallel::makeCluster()} object.
#'   Default 1 (serial). Same convention as \code{isolate_lineage}.
#' @return The \code{cds} with all requested lineages compressed.
#' @export
compress_lineages <- function(cds, lineages = names(cds@lineages), N,
                              method = "sum", ID = FALSE, cl = 1){
  if (length(lineages) == 0) stop("No lineages to compress.", call. = FALSE)
  first_lin <- lineages[1]

  worker <- function(lin) {
    # Only the first lineage streams its per-gene bar; the rest fit silently.
    tmp <- compress_lineage(cds, lineage = lin, N = N, method = method,
                            ID = ID, progress = identical(lin, first_lin))
    list(lineage     = tmp@lineages[[lin]],
         expression  = tmp@expression[[lin]],
         expectation = tmp@expectation[[lin]],
         pseudotime  = tmp@pseudotime[[lin]])
  }

  # Suppress pblapply's own per-lineage bar (uninformative when all lineages
  # start at once); the first lineage's inner bar is shown instead.
  op <- pbapply::pboptions(type = "none")
  on.exit(pbapply::pboptions(op), add = TRUE)
  res <- pbapply::pblapply(lineages, worker, cl = cl)
  names(res) <- lineages

  for (lin in lineages) {
    r <- res[[lin]]
    if (inherits(r, "try-error"))
      stop("compress_lineage failed for lineage '", lin, "': ",
           conditionMessage(attr(r, "condition")), call. = FALSE)
    cds@lineages[[lin]]    <- r$lineage
    cds@expression[[lin]]  <- r$expression
    cds@expectation[[lin]] <- r$expectation
    cds@pseudotime[[lin]]  <- r$pseudotime
  }
  cds
}
