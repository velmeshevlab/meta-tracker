# Lineage expression compression into meta-cells + smoothed expectation curves.

# Start a background process that tails `logfile` and shows new lines as they
# are written. On Windows the reader opens in its OWN console window (inherited
# stdout does not surface in the R console there); on Unix its stdout inherits
# the R console. Returns a processx handle to kill (Unix) or NULL, plus cleans
# up its own window on Windows via the temp script. Returns a handle/list or
# NULL if it can't be started (the run still works, just without the stream).
.start_progress_reader <- function(logfile) {
  tail_code <- paste(
    'a <- commandArgs(TRUE); f <- a[1]',
    'while (!file.exists(f)) Sys.sleep(0.2)',
    'con <- file(f, "r")',
    'repeat {',
    '  l <- readLines(con, warn = FALSE)',
    '  if (length(l)) { cat(l, sep = "\\n"); cat("\\n"); flush(stdout()) }',
    '  Sys.sleep(0.5)',
    '}', sep = "\n")
  rscript <- file.path(R.home("bin"),
                       if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

  if (.Platform$OS.type == "windows") {
    # Open a separate console window running the tailer. Write the code to a
    # temp .R file and launch it in a new window; the returned closure closes it.
    script <- tempfile(fileext = ".R")
    writeLines(tail_code, script)
    title <- "compress_lineages progress"
    ok <- tryCatch({
      system2("cmd", c("/c", "start", shQuote(title),
                       shQuote(rscript), shQuote(script), shQuote(logfile)),
              wait = FALSE)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) {
      message("Could not open a progress window; running without the live stream.")
      return(NULL)
    }
    message("Live progress in a separate window titled '", title,
            "'. (Manual fallback: Get-Content '", logfile, "' -Wait)")
    # Return a closer that kills that window by title and removes the temp script.
    list(kill = function() {
      try(system2("taskkill", c("/FI", shQuote(paste0("WINDOWTITLE eq ", title, "*")), "/T", "/F"),
                  stdout = FALSE, stderr = FALSE), silent = TRUE)
      try(unlink(script), silent = TRUE)
    })
  } else {
    if (!requireNamespace("processx", quietly = TRUE)) {
      message("Install 'processx' for live progress; running without the console stream.")
      return(NULL)
    }
    p <- tryCatch(
      processx::process$new(rscript, c("-e", tail_code, logfile),
                            stdout = "", stderr = ""),
      error = function(e) { message("Could not start progress reader: ",
                                    conditionMessage(e)); NULL })
    if (is.null(p)) return(NULL)
    list(kill = function() try(p$kill(), silent = TRUE))
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
# Append a progress line to `logfile` (read by the background reader), or print
# directly when there is no logfile (serial runs). Package-level so it can be
# shipped to PSOCK workers without dragging in the caller's environment.
.compress_note <- function(logfile, txt) {
  line <- sprintf("[%s pid %d] %s", format(Sys.time(), "%H:%M:%S"), Sys.getpid(), txt)
  if (!is.null(logfile)) cat(line, "\n", file = logfile, append = TRUE)
  else { cat(line, "\n"); utils::flush.console() }
}

# One lineage's fit, run in a worker. `task` carries the small per-lineage
# subset (cds_sub) prepared by the parent, so the full cds is never shipped.
# Package-level so parLapply serialises it against the namespace, not the
# compress_lineages() frame (which holds cds).
.compress_worker <- function(task, N, method, ID, logfile) {
  t0 <- Sys.time()
  .compress_note(logfile, sprintf("compressing lineage '%s' ...", task$lineage))
  r <- .compress_fit(task$cds_sub, task$lineage, task$updated_pt,
                     task$cell_barcodes, N = N, method = method,
                     ID = ID, progress = FALSE)
  .compress_note(logfile, sprintf("finished lineage '%s' (%s)", task$lineage,
                 format(round(difftime(Sys.time(), t0), 1))))
  r
}

# Compress a lineage: project + subset via get_lineage_object, then fit.
.compress_expression <- function(cds, lineage, N, method = "sum", ID = FALSE, progress = TRUE){
  cds_sub       <- get_lineage_object(cds, lineage)
  updated_pt    <- cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  cell_barcodes <- .lineage_cells(cds@lineages[[lineage]])
  .compress_fit(cds_sub, lineage, updated_pt, cell_barcodes,
                N = N, method = method, ID = ID, progress = progress)
}

.compress_fit <- function(cds_sub, lineage, updated_pt, cell_barcodes, N, method = "sum", ID = FALSE, progress = TRUE){
  # Meta-cell binning + per-gene fitting on a single lineage's subset (cds_sub).
  # cds_sub is produced by get_lineage_object() upstream; operating on it (not the
  # full cds) is what keeps memory small when workers run in parallel.
  sel.cells = cell_barcodes[cell_barcodes %in% colnames(cds_sub)]
  if (length(sel.cells) == 0)
    stop("Lineage '", lineage, "' has no cells present in the subset.", call. = FALSE)
  cds_subset = cds_sub[, sel.cells]
  lineage_entry <- list(name = cell_barcodes, updated_pt = updated_pt)
  #prepare raw count matrix
  exp = as.data.frame(as.matrix(exprs(cds_subset)))
  exp_sum <- t(exp)
  exp_sum = exp_sum[, rownames(cds_sub)]
  #prepare size factor
  size_factor <- (pData(cds_subset)[, 'Size_Factor'])
  size_factor <- size_factor[rownames(exp_sum)]
  #prepare pseudotime
  pt <- updated_pt
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
      exp_mean = exp_mean[,rownames(cds_sub)]
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
        "lineage" = lineage_entry, "expression" = list("sum" = meta_sum_ordered, "mean" = meta_mean_ordered), "expectation" = fit_list_2, "pseudotime"  = list("real" = meta_sum_ordered$pseudotime, "scaled" = d)))
    }
    return(list("lineage" = lineage_entry, "expression" = meta_sum_ordered, "expectation" = fit_list_2, "pseudotime"  = list("real" = meta_sum_ordered$pseudotime, "scaled" = d)))
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
#' into one object. Lineages are processed in chunks sized to the worker count,
#' printing cumulative progress from the parent (\code{"6/11 lineages are being
#' processed"}) before each chunk. Parallelism depends on \code{cl}: an integer
#' \code{> 1} forks on Unix/macOS; on Windows (no fork) an integer \code{> 1}
#' auto-creates a PSOCK cluster of that size; or pass your own
#' \code{makeCluster()} object. PSOCK workers each receive a copy of \code{cds}
#' and do not forward console output, so per-lineage detail on Windows comes
#' from the parent progress prints and the optional \code{log} file.
#'
#' @param cds A \code{metatracker_data_set}.
#' @param lineages Lineage names to compress (default: all in \code{cds@lineages}).
#' @param N Number of meta-cells per lineage.
#' @param method "sum" (default) or any other value to also compute the mean matrix.
#' @param ID Passed through to \code{compress_lineage} (default FALSE).
#' @param cl Integer worker count (fork on Unix, serial on Windows) or a
#'   \code{parallel::makeCluster()} object. Default 1 (serial). Same convention
#'   as \code{isolate_lineage}.
#' @param log Optional path to a progress log file. Each worker appends a line
#'   when it starts and finishes a lineage (with PID and elapsed time). Useful
#'   Live monitoring uses a temporary progress file that each worker appends to;
#'   a background reader echoes it to the console and the file is deleted when
#'   the run finishes.
#' @param monitor Show live per-lineage progress from a background reader
#'   (default TRUE). Requires the \code{processx} package; if unavailable, the
#'   run still works but without the live console stream.
#' @return The \code{cds} with all requested lineages compressed.
#' @export
compress_lineages <- function(cds, lineages = names(cds@lineages), N,
                              method = "sum", ID = FALSE, cl = 1, monitor = TRUE){
  if (length(lineages) == 0) stop("No lineages to compress.", call. = FALSE)
  cat(sprintf("Compressing %d lineage(s): %s\n",
              length(lineages), paste(lineages, collapse = ", ")))
  utils::flush.console()

  parallel_run <- inherits(cl, "cluster") || (is.numeric(cl) && cl > 1)

  # Progress transport: workers can't print to the console (esp. PSOCK), but they
  # can append to a file that a background reader echoes to the console.
  logfile <- NULL
  reader  <- NULL
  if (parallel_run && isTRUE(monitor)) {
    logfile <- tempfile("compress_progress_", fileext = ".log")
    file.create(logfile)
    reader  <- .start_progress_reader(logfile)
  }
  on.exit({
    if (!is.null(reader)) { Sys.sleep(1); try(reader$kill(), silent = TRUE) }
    if (!is.null(logfile)) unlink(logfile)   # remove the log at the end
  }, add = TRUE)

  # Parent-side projection + subset for one lineage (cheap: ~seconds). Produces a
  # small per-lineage object so workers never receive the full cds.
  prepare <- function(lin) {
    cds_sub       <- get_lineage_object(cds, lin)
    updated_pt    <- cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
    cell_barcodes <- .lineage_cells(cds@lineages[[lin]])
    list(lineage = lin, cds_sub = cds_sub,
         updated_pt = updated_pt, cell_barcodes = cell_barcodes)
  }

  # Quiet pbapply bars during the run; progress comes from the reader + prints.
  op <- pbapply::pboptions(type = "none")
  on.exit(pbapply::pboptions(op), add = TRUE)

  # Resolve a worker backend.
  #  - a cluster passed as `cl`         -> use it as-is
  #  - integer cl > 1 on Unix/macOS     -> fork (mclapply, fresh process per lineage)
  #  - integer cl > 1 on Windows        -> auto-create a PSOCK cluster (no fork on Windows)
  #  - otherwise                        -> serial
  clobj <- if (inherits(cl, "cluster")) cl else NULL
  if (is.null(clobj) && is.numeric(cl) && cl > 1 && .Platform$OS.type == "windows") {
    clobj <- parallel::makeCluster(as.integer(cl))
    parallel::clusterEvalQ(clobj, suppressMessages(library(metatracker)))
    on.exit(parallel::stopCluster(clobj), add = TRUE)   # only stop clusters we created
    message("Windows: created a ", as.integer(cl),
            "-worker PSOCK cluster (each worker holds a copy of cds).")
  }

  n <- length(lineages)
  n_workers <- if (!is.null(clobj)) length(clobj)
               else if (is.numeric(cl)) max(1L, as.integer(cl)) else 1L
  chunks <- split(lineages, ceiling(seq_along(lineages) / n_workers))

  run_chunk <- function(tasks) {
    if (!is.null(clobj)) {
      # parLapply serializes the FUN's environment; use the package-level
      # .compress_worker (env = namespace, no cds) so only the small tasks ship.
      parallel::parLapply(clobj, tasks, .compress_worker,
                          N = N, method = method, ID = ID, logfile = logfile)
    } else if (is.numeric(cl) && cl > 1 && .Platform$OS.type != "windows") {
      parallel::mclapply(tasks, .compress_worker,
                         N = N, method = method, ID = ID, logfile = logfile,
                         mc.cores = cl, mc.preschedule = FALSE)
    } else {
      lapply(tasks, .compress_worker,
             N = N, method = method, ID = ID, logfile = logfile)
    }
  }

  res <- list()
  done <- 0L
  for (ch in chunks) {
    done <- done + length(ch)
    cat(sprintf("%d/%d lineages are being processed\n", done, n)); utils::flush.console()
    tasks <- lapply(ch, prepare)          # parent: project + subset each lineage (cheap)
    res   <- c(res, run_chunk(tasks))     # ship only the small subsets to workers
  }
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
