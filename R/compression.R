# Lineage expression compression into meta-cells + smoothed expectation curves.

# Start a background process that tails `logfile` and shows new lines as they
# are written. On Windows the reader opens in its OWN console window (inherited
# stdout does not surface in the R console there); on Unix its stdout inherits
# the R console. Returns a processx handle to kill (Unix) or NULL, plus cleans
# up its own window on Windows via the temp script. Returns a handle/list or
# NULL if it can't be started (the run still works, just without the stream).
.start_progress_reader <- function(logfile, title = "compress_lineages progress") {
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
.compress_worker <- function(prep, N, method, ID, logfile) {
  t0 <- Sys.time()
  .compress_note(logfile, sprintf("compressing lineage '%s' ...", prep$lineage))
  r <- .compress_fit_prep(prep, N = N, method = method, ID = ID,
                          progress = FALSE, progress_log = logfile)
  .compress_note(logfile, sprintf("finished lineage '%s' (%s)", prep$lineage,
                 format(round(difftime(Sys.time(), t0), 1))))
  r
}

# Compress a lineage: aggregate to meta-cells (parent), then fit (worker).
# compress_lineage() (single lineage) runs both here; compress_lineages() calls
# .compress_prep() in the parent and .compress_fit_prep() in workers, so workers
# receive only the small aggregated matrices, never the per-cell data.
.compress_expression <- function(cds, lineage, N, method = "sum", ID = FALSE, progress = TRUE){
  prep <- .compress_prep(cds, lineage, N = N, method = method)
  .compress_fit_prep(prep, N = N, method = method, ID = ID, progress = progress)
}

# Parent-side: project + subset, then bin cells into N meta-cells and aggregate
# gene counts with a SPARSE matrix multiply (counts %*% indicator) -- no dense
# cells x genes matrix is ever formed. Returns the small N x genes meta matrices.
.compress_prep <- function(cds, lineage, N, method = "sum"){
  cds_sub       <- get_lineage_object(cds, lineage)
  updated_pt    <- cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  cell_barcodes <- .lineage_cells(cds@lineages[[lineage]])
  sel.cells     <- cell_barcodes[cell_barcodes %in% colnames(cds_sub)]
  if (length(sel.cells) == 0)
    stop("Lineage '", lineage, "' has no cells present in the subset.", call. = FALSE)
  cds_subset <- cds_sub[, sel.cells]

  counts <- exprs(cds_subset)                    # genes x cells (kept sparse)
  counts <- counts[rownames(cds_sub), ]          # gene order (rows)
  cells  <- colnames(counts)
  Ccell  <- length(cells)
  if (N >= Ccell)
    stop(sprintf("out of boundary: N (%d) must be smaller than number of cells (%d)", N, Ccell))

  # per-cell vectors aligned to counts columns
  sf   <- pData(cds_subset)[, "Size_Factor"]; names(sf) <- colnames(cds_subset); sf <- sf[cells]
  ptv  <- as.numeric(updated_pt[cells])
  umap <- reducedDims(cds_subset)[["UMAP"]][cells, , drop = FALSE]

  # meta-cell assignment: bin cells by pseudotime rank (order-independent, matches
  # the original cut(rank(pseudotime), N))
  meta_cell <- cut(rank(ptv), breaks = N, labels = FALSE)
  Ind <- Matrix::sparseMatrix(i = seq_len(Ccell), j = meta_cell, x = 1,
                              dims = c(Ccell, N))          # cells x N indicator
  n_cells  <- as.integer(Matrix::colSums(Ind))
  pt_mean  <- as.numeric(tapply(ptv,       meta_cell, mean)[as.character(seq_len(N))])
  sf_sum   <- as.numeric(tapply(as.numeric(sf), meta_cell, sum)[as.character(seq_len(N))])
  um1_mean <- as.numeric(tapply(umap[, 1], meta_cell, mean)[as.character(seq_len(N))])
  um2_mean <- as.numeric(tapply(umap[, 2], meta_cell, mean)[as.character(seq_len(N))])

  # summed gene counts per meta-cell: crossprod(Ind, counts_ct) = t(Ind) %*% (cells x genes)
  # = N x genes. Using Matrix::crossprod / Matrix::t keeps S4 dispatch explicit so it
  # can't fall through to base t.default on a sparse object.
  counts_ct <- Matrix::t(counts)                           # cells x genes (sparse)
  gene_sum  <- as.matrix(Matrix::crossprod(Ind, counts_ct))  # N x genes
  ord <- order(pt_mean)
  meta_sum_ordered <- data.frame(meta_cell = ord, n_cells = n_cells[ord],
                                 pseudotime = pt_mean[ord], umap_1 = um1_mean[ord],
                                 umap_2 = um2_mean[ord], size_factor = sf_sum[ord],
                                 check.names = FALSE)
  meta_sum_ordered <- cbind(meta_sum_ordered, gene_sum[ord, , drop = FALSE])

  meta_mean_ordered <- NULL
  if (method != "sum") {
    # per-cell normalisation then mean per meta-cell = (sum of counts/sf) / n_cells
    counts_norm_ct <- Matrix::Diagonal(x = 1 / as.numeric(sf)) %*% counts_ct  # scale each cell (row)
    gene_norm_sum  <- as.matrix(Matrix::crossprod(Ind, counts_norm_ct))       # N x genes
    gene_mean      <- gene_norm_sum / n_cells               # row m divided by its n_cells
    meta_mean_ordered <- data.frame(meta_cell = ord, n_cells = n_cells[ord],
                                    pseudotime = pt_mean[ord], umap_1 = um1_mean[ord],
                                    umap_2 = um2_mean[ord], check.names = FALSE)
    meta_mean_ordered <- cbind(meta_mean_ordered, gene_mean[ord, , drop = FALSE])
  }

  list(lineage = lineage, meta_sum = meta_sum_ordered, meta_mean = meta_mean_ordered,
       lineage_entry = list(name = cell_barcodes, updated_pt = updated_pt))
}

# Worker-side: per-gene quasipoisson fit on the (small) aggregated meta matrix.
.compress_fit_prep <- function(prep, N, method = "sum", ID = FALSE, progress = TRUE,
                               progress_log = NULL){
  lineage          <- prep$lineage
  lineage_entry    <- prep$lineage_entry
  meta_sum_ordered <- prep$meta_sum
  if (ID == FALSE) {
    mat   <- meta_sum_ordered[, 7:ncol(meta_sum_ordered)]
    model <- expression ~ splines::ns(pseudotime, df = 7) + offset(log(size_factor))
    sf_meta <- meta_sum_ordered$size_factor
    d <- (meta_sum_ordered$pseudotime - min(meta_sum_ordered$pseudotime)) /
         (max(meta_sum_ordered$pseudotime) - min(meta_sum_ordered$pseudotime))
    predict_pt <- seq(0, 1, length.out = N)
    mat_m <- as.matrix(mat)
    genes <- colnames(mat_m)
    .fitcol <- function(i) {
      .fit_m3(exp.sel = mat_m[, i], pt = d, size_factor = sf_meta,
              predict_pt = predict_pt, lineage = lineage, model = model, N = N)
    }
    fit_list_2 <- if (isTRUE(progress)) {
      op <- pbapply::pboptions(type = "timer")
      on.exit(pbapply::pboptions(op), add = TRUE)
      pbapply::pbsapply(seq_along(genes), .fitcol)
    } else if (!is.null(progress_log)) {
      ng <- length(genes); step <- max(1L, ng %/% 50L); last <- -1L
      out <- vector("list", ng)
      for (i in seq_len(ng)) {
        out[[i]] <- .fitcol(i)
        if (i %% step == 0L || i == ng) {
          pct <- as.integer(round(100 * i / ng))
          if (pct != last) { .compress_note(progress_log, sprintf("%s: %d%%", lineage, pct)); last <- pct }
        }
      }
      m <- do.call(cbind, out); colnames(m) <- genes; m
    } else {
      sapply(seq_along(genes), .fitcol)
    }
    colnames(fit_list_2) <- genes
    fit_list_2 <- apply(fit_list_2, 2, as.numeric)
    if (method != "sum") {
      return(list("lineage" = lineage_entry,
                  "expression" = list("sum" = meta_sum_ordered, "mean" = prep$meta_mean),
                  "expectation" = fit_list_2,
                  "pseudotime" = list("real" = meta_sum_ordered$pseudotime, "scaled" = d)))
    }
    return(list("lineage" = lineage_entry, "expression" = meta_sum_ordered,
                "expectation" = fit_list_2,
                "pseudotime" = list("real" = meta_sum_ordered$pseudotime, "scaled" = d)))
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

  # Quiet pbapply bars during the run; progress comes from the reader + prints.
  op <- pbapply::pboptions(type = "none")
  on.exit(pbapply::pboptions(op), add = TRUE)

  # Guard: parLapply serialises the worker with its environment. If this package
  # is source()'d into .GlobalEnv rather than installed, .compress_worker's
  # environment IS .GlobalEnv -- which holds `cds` -- and the full object would
  # be shipped to every worker (OOM). Detect and refuse rather than OOM silently.
  if ((inherits(cl, "cluster") || (is.numeric(cl) && cl > 1 &&
        .Platform$OS.type == "windows")) &&
      environmentName(environment(.compress_worker)) == "R_GlobalEnv") {
    stop("compress_lineages is running from .GlobalEnv (source()'d), which would ",
         "ship the payloads/cds to each PSOCK worker from the global env. Install ",
         "and library(metatracker) instead of source()-ing the files, then retry.",
         call. = FALSE)
  }

  # Resolve a worker backend.
  #  - a cluster passed as `cl`         -> use it as-is
  #  - integer cl > 1 on Unix/macOS     -> fork (mclapply)
  #  - integer cl > 1 on Windows        -> auto-create a PSOCK cluster (no fork on Windows)
  #  - otherwise                        -> serial
  clobj <- if (inherits(cl, "cluster")) cl else NULL
  if (is.null(clobj) && is.numeric(cl) && cl > 1 && .Platform$OS.type == "windows") {
    clobj <- parallel::makeCluster(as.integer(cl))
    parallel::clusterEvalQ(clobj, suppressMessages(library(metatracker)))
    on.exit(parallel::stopCluster(clobj), add = TRUE)   # only stop clusters we created
    message("Windows: created a ", as.integer(cl), "-worker PSOCK cluster ",
            "(each worker receives one lineage's exp_sum payload, not the full cds).")
  }

  n <- length(lineages)
  n_workers <- if (!is.null(clobj)) length(clobj)
               else if (is.numeric(cl)) max(1L, as.integer(cl)) else 1L
  # Chunk lineages to the worker count. Payloads are built per chunk and freed
  # before the next, so peak memory is ~2 x (n_workers x one exp_sum) + cds
  # rather than all payloads at once.
  chunks <- split(lineages, ceiling(seq_along(lineages) / n_workers))

  run_chunk <- function(preps) {
    if (!is.null(clobj)) {
      parallel::parLapply(clobj, preps, .compress_worker,
                          N = N, method = method, ID = ID, logfile = logfile)
    } else if (is.numeric(cl) && cl > 1 && .Platform$OS.type != "windows") {
      parallel::mclapply(preps, .compress_worker,
                         N = N, method = method, ID = ID, logfile = logfile,
                         mc.cores = cl, mc.preschedule = FALSE)
    } else {
      lapply(preps, .compress_worker,
             N = N, method = method, ID = ID, logfile = logfile)
    }
  }

  res  <- list()
  done <- 0L
  for (ch in chunks) {
    done <- done + length(ch)
    cat(sprintf("%d/%d lineages are being processed\n", done, n)); utils::flush.console()
    # Phase 1 (this chunk): build the small per-lineage payloads in the parent.
    preps <- lapply(ch, function(lin) .compress_prep(cds, lin, N = N, method = method))
    names(preps) <- ch
    # Phase 2 (this chunk): fit in parallel over the payloads.
    res <- c(res, run_chunk(preps))
    rm(preps); gc(verbose = FALSE)          # free this chunk before building the next
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
