# Between-lineage differential expression (tradeSeq).
#
# lineage_specific_genes_par() fits one tradeSeq GAM across all lineages'
# meta-cells, runs patternTest / diffEndTest, and derives per-lineage
# fold-changes. Relative to the first version: fitGAM can run in parallel
# (ncores > 1) by fitting gene-blocks in parallel and be checkpointed to disk
# (gam_file); AUC /
# fold-changes are predicted once per lineage (O(L)) rather than per pair
# (O(L^2)); an optional low-expression prefilter skips genes that cannot pass
# the downstream filters; and the fitted model is passed through to the
# per-lineage step (the original dropped it).


# --- Build the fitGAM inputs (counts / pseudotime / cellWeights / offset) -----
.prepare_gam_input <- function(cds,
                              lineages       = names(cds@lineages),
                              min_metacells  = 0,   # 0 = no filtering (original behaviour)
                              min_count      = 1) {

  counts_list     <- list()
  all_metacells   <- character(0)
  all_size_factor <- numeric(0)
  lineage_list    <- list()
  pt_list         <- list()

  for (lineage in lineages) {
    d_raw <- cds@expression[[lineage]]
    if (is.list(d_raw) && !is.null(d_raw[["sum"]])) d_raw <- d_raw[["sum"]]
    metacells <- paste0(lineage, "_", seq_len(nrow(d_raw)))

    d <- t(as.matrix(sapply(d_raw[, 7:ncol(d_raw)], as.numeric)))
    d <- d[rownames(cds), , drop = FALSE]
    colnames(d) <- metacells

    counts_list[[lineage]]  <- d
    all_metacells           <- c(all_metacells, metacells)
    all_size_factor         <- c(all_size_factor, d_raw$size_factor)
    lineage_list[[lineage]] <- metacells

    pt        <- d_raw$pseudotime
    pt        <- (pt - min(pt)) / (max(pt) - min(pt))
    names(pt) <- metacells
    pt_list[[lineage]] <- pt
  }

  counts <- do.call(cbind, counts_list)
  counts <- counts[rownames(cds), , drop = FALSE]

  cellWeights <- matrix(0, nrow = length(all_metacells), ncol = length(lineages),
                        dimnames = list(all_metacells, lineages))
  for (ln in lineages) {
    cellWeights[, ln] <- as.numeric(all_metacells %in% lineage_list[[ln]])
  }

  pseudotime <- matrix(0, nrow = length(all_metacells), ncol = length(lineages),
                       dimnames = list(all_metacells, lineages))
  for (ln in lineages) {
    pseudotime[names(pt_list[[ln]]), ln] <- pt_list[[ln]]
  }

  # --- optional prefilter -----------------------------------------------------
  # Keep a gene if it is detected in >= min_metacells metacells of at least one
  # lineage. Genes that are near-zero everywhere cost the same GAM fit time as
  # informative ones but can never pass the downstream filters.
  if (min_metacells > 0) {
    keep <- rep(FALSE, nrow(counts))
    for (ln in lineages) {
      n_det <- Matrix::rowSums(counts[, lineage_list[[ln]], drop = FALSE] >= min_count)
      keep  <- keep | (n_det >= min_metacells)
    }
    message(sprintf("Prefilter: keeping %d / %d genes (detected in >= %d metacells of >= 1 lineage)",
                    sum(keep), length(keep), min_metacells))
    counts <- counts[keep, , drop = FALSE]
  }

  list(counts      = counts,
       pseudotime  = pseudotime,
       cellWeights = cellWeights,
       offset      = log(all_size_factor))
}


# --- Fit the GAM, in parallel, with a checkpoint ------------------------------
# Append a line to the shared progress log (or message() if none). Short appends
# to a file opened O_APPEND are effectively atomic, so worker lines don't tangle.
.block_log <- function(logfile, txt) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", txt, "\n")
  if (is.null(logfile)) message(txt) else cat(line, file = logfile, append = TRUE)
}

# Fit one gene-block serially and return its fitGAM SCE. Top-level (namespace)
# function so PSOCK workers resolve it by reference instead of serialising the
# caller's frame (which holds the full counts matrix). `job` carries the block's
# id and its counts so each worker can log which block it is on.
.fit_gam_block <- function(job, pseudotime, cellWeights, offset, nknots, family,
                           n_total = NA, logfile = NULL) {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
  id <- job$id; cb <- job$counts
  t0 <- Sys.time()
  .block_log(logfile, sprintf("block %d/%s START  (%d genes)", id, n_total, nrow(cb)))
  res <- tradeSeq::fitGAM(counts = cb, pseudotime = pseudotime,
                          cellWeights = cellWeights, U = NULL, nknots = nknots,
                          offset = offset, family = family, parallel = FALSE, sce = TRUE)
  .block_log(logfile, sprintf("block %d/%s DONE   (%d genes, %s)",
                              id, n_total, nrow(cb), format(round(Sys.time() - t0, 1))))
  res
}

.fit_lineage_gam <- function(cds,
                            lineages      = names(cds@lineages),
                            nknots        = 6,
                            ncores        = 1,
                            nblocks       = NULL,   # gene-blocks; default = ncores
                            gam_file      = NULL,   # e.g. "gamlist.rds"
                            log_file      = NULL,   # per-block progress log
                            min_metacells = 0,
                            min_count     = 1,
                            family        = "nb",
                            ...) {

  # Reuse an existing checkpoint if present.
  if (!is.null(gam_file) && file.exists(gam_file)) {
    message("Loading existing fitGAM object from ", gam_file)
    return(readRDS(gam_file))
  }

  inp <- .prepare_gam_input(cds, lineages = lineages,
                           min_metacells = min_metacells, min_count = min_count)

  ng <- nrow(inp$counts)
  message(sprintf("Testing %d genes across %d metacells, %d lineages",
                  ng, ncol(inp$counts), length(lineages)))

  # Serial BLAS in the parent too (workers set it themselves).
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }

  # Split genes (rows) into contiguous blocks. fitGAM's per-call overhead is
  # small next to fitting thousands of genes, so coarse blocking (~one per core)
  # keeps scheduling overhead low; raise nblocks for better load balance AND
  # more frequent progress updates (the log/bar advance once per finished block).
  if (is.null(nblocks)) nblocks <- max(1L, as.integer(ncores))
  nblocks <- min(nblocks, ng)
  grp     <- cut(seq_len(ng), breaks = nblocks, labels = FALSE)
  blocks  <- split(seq_len(ng), grp)
  # jobs bundle each block's id with its counts (references, no copy) so the id
  # travels with the data to whichever worker fits it.
  jobs <- lapply(seq_along(blocks),
                 function(k) list(id = k, counts = inp$counts[blocks[[k]], , drop = FALSE]))
  n_total <- length(jobs)
  ps <- inp$pseudotime; cw <- inp$cellWeights; off <- inp$offset
  inp$counts <- NULL  # free the parent's full copy; jobs hold the data now

  # Per-block progress log. Default to a temp file when running in parallel
  # (worker console output is invisible); tail it to watch blocks complete.
  if (is.null(log_file) && ncores > 1 && n_total > 1) {
    log_file <- tempfile("metatracker_fitgam_", fileext = ".log")
  }
  if (!is.null(log_file)) {
    file.create(log_file)
    tail_hint <- if (.Platform$OS.type == "windows")
      sprintf("Get-Content '%s' -Wait", log_file) else sprintf("tail -f '%s'", log_file)
    message(sprintf("Per-block progress log: %s\n  (watch it live with:  %s )", log_file, tail_hint))
    # Open a separate window that streams the log (own console on Windows,
    # inherited stdout via processx on Unix). Closed on exit.
    reader <- .start_progress_reader(log_file, title = "lineage_specific_genes progress")
    if (!is.null(reader)) on.exit(reader$kill(), add = TRUE)
  }


  t0 <- Sys.time()
  if (ncores > 1 && length(jobs) > 1) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterEvalQ(cl, {
        library(tradeSeq)
        if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
          RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1)
        }
      })
      sce_list <- pbapply::pblapply(jobs, .fit_gam_block, pseudotime = ps,
                                    cellWeights = cw, offset = off, nknots = nknots,
                                    family = family, n_total = n_total,
                                    logfile = log_file, cl = cl)
    } else {
      # fork: workers share the blocks via copy-on-write, no serialisation cost
      sce_list <- pbapply::pblapply(jobs, .fit_gam_block, pseudotime = ps,
                                    cellWeights = cw, offset = off, nknots = nknots,
                                    family = family, n_total = n_total,
                                    logfile = log_file, cl = ncores)
    }
  } else {
    sce_list <- pbapply::pblapply(jobs, .fit_gam_block, pseudotime = ps,
                                  cellWeights = cw, offset = off, nknots = nknots,
                                  family = family, n_total = n_total,
                                  logfile = log_file)
  }
  message("fitGAM finished in ", format(Sys.time() - t0))

  # Reassemble the full-gene SCE. Per-cell design (colData/metadata: dm, X,
  # knots) is identical across blocks because every block saw the same
  # pseudotime/cellWeights; rbind stacks the per-gene rowData (beta, Sigma).
  gamlist <- if (length(sce_list) == 1) sce_list[[1]] else do.call(rbind, sce_list)

  if (!is.null(gam_file)) {
    saveRDS(gamlist, gam_file)
    message("Saved fitGAM object to ", gam_file)
  }
  gamlist
}


# --- Predict each lineage's AUC ONCE (instead of once per pair) ---------------
.compute_lineage_auc <- function(models, lineages, genes, N = 1000, knots = NULL) {

  dm <- colData(models)$tradeSeq$dm
  X  <- colData(models)$tradeSeq$X
  slingshotColData <- colData(models)$crv
  pseudotime <- slingshotColData[, grep(x = colnames(slingshotColData),
                                        pattern = "pseudotime")]
  betaMat   <- rowData(models)$tradeSeq$beta[[1]]
  beta      <- as.matrix(betaMat[genes, , drop = FALSE])
  knotPoints <- S4Vectors::metadata(models)$tradeSeq$knots

  time_points <- seq(from = 0, to = 1, length.out = N)
  idx <- seq_len(N)
  if (!is.null(knots)) {
    t1  <- unname(knotPoints[knots[[1]]])
    t2  <- unname(knotPoints[knots[[2]]])
    idx <- which(time_points >= t1 & time_points <= t2)
  }
  tp <- time_points[idx]

  auc <- matrix(NA_real_, nrow = length(genes), ncol = length(lineages),
                dimnames = list(genes, lineages))

  for (lin in lineages) {
    lineage_id <- which(lineages == lin)
    df   <- tradeSeq:::.getPredictRangeDf(dm, lineage_id, nPoints = N)
    Xdf  <- tradeSeq:::predictGAM(lpmatrix = X, df = df, pseudotime = pseudotime)
    yhat <- exp(Xdf %*% t(beta) + df$offset)
    yhat <- yhat[idx, , drop = FALSE]
    auc[, lin] <- apply(yhat, 2, function(y) pracma::trapz(tp, y))
  }
  auc
}


# Reproduces .calculate_FC_tradeseq()'s return value from a precomputed AUC matrix
.fc_from_auc <- function(auc, test_lineage, lineages, fc_names = NULL) {
  others <- lineages[lineages != test_lineage]   # same order as the original loop
  out <- log2(auc[, rep(test_lineage, length(others)), drop = FALSE] /
                auc[, others, drop = FALSE])
  if (length(lineages) == 2) {
    colnames(out) <- paste0("log2FC_", test_lineage, "vs", others)
  } else {
    colnames(out) <- paste0(fc_names, "_pattern")
  }
  out
}


# --- Patched v2: accepts an optional precomputed AUC matrix -------------------
# Identical to your .lineage_specific_genes_v2 except that the two
# .calculate_FC_tradeseq() calls become .fc_from_auc() when `auc` is supplied.
# If auc = NULL it falls back to the original behaviour using `model`.
.lineage_specific_genes_v2 <- function(test_lineage, cds, model, pattern, diffend,
                                      genes = NULL, lineages, auc = NULL, N = 1000) {
  print(paste0("Testing ", test_lineage))

  .get_FCs <- function(gene_names, fc_names = NULL) {
    if (!is.null(auc)) {
      .fc_from_auc(auc[gene_names, , drop = FALSE], test_lineage, lineages, fc_names)
    } else {
      .calculate_FC_tradeseq(models = model, test_lineage = test_lineage,
                            lineages = lineages, fc_names = fc_names,
                            genes = gene_names, N = N)
    }
  }

  if (length(lineages) > 2) {
    index <- which(test_lineage == lineages)
    p_list <- wd_list <- fc_list <- c()
    p_names <- wd_names <- fc_names <- c()
    flip_fc <- c()
    for (lineage in lineages) {
      index2 <- which(lineage == lineages)
      if (index != index2) {
        if (index < index2) {
          p_name  <- paste0("pvalue_",   index,  "vs", index2)
          wd_name <- paste0("waldStat_", index,  "vs", index2)
          fc_name <- paste0("logFC",     index,  "_",  index2)
          flip_fc <- c(flip_fc, FALSE)
        } else {
          p_name  <- paste0("pvalue_",   index2, "vs", index)
          wd_name <- paste0("waldStat_", index2, "vs", index)
          fc_name <- paste0("logFC",     index2, "_",  index)
          flip_fc <- c(flip_fc, TRUE)
        }
        p_list   <- c(p_list, p_name)
        wd_list  <- c(wd_list, wd_name)
        fc_list  <- c(fc_list, fc_name)
        p_names  <- c(p_names,  paste0("pvalue_",   test_lineage, "vs", lineage))
        wd_names <- c(wd_names, paste0("waldStat_", test_lineage, "vs", lineage))
        fc_names <- c(fc_names, paste0("log2FC_",   test_lineage, "vs", lineage))
      }
    }

    # ---- pattern test ----
    p_values_p <- pattern[, p_list]
    wd_values  <- pattern[, wd_list]
    ref        <- rownames(pattern)
    p_values_p <- p_values_p[ref, , drop = FALSE]
    wd_values  <- wd_values[ref, , drop = FALSE]
    stopifnot(identical(ref, rownames(p_values_p)),
              identical(ref, rownames(wd_values)))
    pattern_sel <- cbind(pattern[, c(1, 3), drop = FALSE], p_values_p, wd_values)
    colnames(pattern_sel) <- c("waldStat_combined", "pvalue_combined", p_names, wd_names)
    colnames(pattern_sel) <- paste0(colnames(pattern_sel), "_pattern")

    # ---- diffend test ----
    p_values_d <- diffend[, p_list]
    wd_values  <- diffend[, wd_list]
    fc_values  <- diffend[, fc_list]
    fc_values[, flip_fc] <- -fc_values[, flip_fc]
    fc_values  <- as.data.frame(fc_values / log(2))
    average_FC <- apply(fc_values, 1, .get_average_FC)

    ref        <- rownames(diffend)
    p_values_d <- p_values_d[ref, , drop = FALSE]
    wd_values  <- wd_values[ref, , drop = FALSE]
    fc_values  <- fc_values[ref, , drop = FALSE]
    average_FC <- average_FC[ref]
    stopifnot(identical(ref, rownames(p_values_d)),
              identical(ref, rownames(wd_values)),
              identical(ref, rownames(fc_values)),
              identical(ref, names(average_FC)))

    diffend_sel <- cbind(diffend[, c(1, 3), drop = FALSE], p_values_d,
                         wd_values, fc_values, average_FC)
    colnames(diffend_sel) <- c("waldStat_combined", "pvalue_combined",
                               p_names, wd_names, fc_names, "average_FC")
    colnames(diffend_sel) <- paste0(colnames(diffend_sel), "_diffend")
    pattern_sel <- pattern_sel[rownames(diffend_sel), ]

    gene_names <- rownames(pattern_sel)
    FCs_sel    <- .get_FCs(gene_names, fc_names)
    pattern_sel <- pattern_sel[rownames(FCs_sel), ]
    average_FC  <- apply(FCs_sel, 1, .get_average_FC)[rownames(FCs_sel)]
    if (!identical(rownames(FCs_sel), names(average_FC))) {
      stop("Error: rownames(FCs_sel) and names(average_FC) do not match.")
    }
    pattern_fin <- cbind(pattern_sel, FCs_sel, average_FC)
    colnames(pattern_fin) <- c(colnames(pattern_sel), colnames(FCs_sel), "average_FC_pattern")
    list("pattern_test" = as.data.frame(pattern_fin),
         "diffend_test" = as.data.frame(diffend_sel))

  } else {
    pattern_sel <- pattern[, c("waldStat", "pvalue")]
    diffend_sel <- diffend[, c("waldStat", "pvalue", "logFC1_2")]
    diffend_sel$logFC1_2 <- diffend_sel$logFC1_2 / log(2)
    if (test_lineage == lineages[2]) diffend_sel$logFC1_2 <- -diffend_sel$logFC1_2

    gene_names  <- rownames(pattern_sel)
    FCs_sel     <- .get_FCs(gene_names)
    pattern_sel <- pattern_sel[rownames(FCs_sel), ]
    diffend_sel <- diffend_sel[rownames(FCs_sel), ]
    pattern_fin <- as.data.frame(cbind(pattern_sel, FCs_sel))
    diffend_fin <- as.data.frame(diffend_sel)
    colnames(pattern_fin) <- c("waldStat_combined_pattern", "pvalue_combined_pattern", "average_FC_pattern")
    colnames(diffend_fin) <- c("waldStat_combined_diffend", "pvalue_combined_diffend", "average_FC_diffend")
    list("pattern_test" = pattern_fin, "diffend_test" = diffend_fin)
  }
}


# --- Main entry point ---------------------------------------------------------
#' @export
lineage_specific_genes_par <- function(cds,
                                       lineages       = names(cds@lineages),
                                       nknots         = 6,
                                       ncores         = 1,
                                       nblocks        = NULL,
                                       gam_file       = NULL,
                                       log_file       = NULL,
                                       gamlist        = NULL,
                                       min_metacells  = 0,
                                       min_count      = 1,
                                       precompute_auc = TRUE,
                                       N              = 1000,
                                       family         = "nb") {

  if (is.null(gamlist)) {
    gamlist <- .fit_lineage_gam(cds, lineages = lineages, nknots = nknots,
                               ncores = ncores, nblocks = nblocks, gam_file = gam_file,
                               log_file = log_file, min_metacells = min_metacells,
                               min_count = min_count, family = family)
  }

  message("Running patternTest ...")
  pattern <- tradeSeq::patternTest(models = gamlist, global = TRUE, pairwise = TRUE)
  message("Running diffEndTest ...")
  diffend <- tradeSeq::diffEndTest(models = gamlist, global = TRUE, pairwise = TRUE)

  auc <- NULL
  if (precompute_auc) {
    message("Precomputing per-lineage AUC ...")
    common <- intersect(rownames(pattern), rownames(diffend))
    auc <- .compute_lineage_auc(gamlist, lineages = lineages, genes = common, N = N)
  }

  out <- sapply(lineages, .lineage_specific_genes_v2,
                cds      = cds,
                model    = gamlist,   # <-- THE FIX
                pattern  = pattern,
                diffend  = diffend,
                lineages = lineages,
                auc      = auc,
                N        = N,
                simplify = FALSE)

  cds@lineage_genes <- out
  cds
}


# --- Fallback FC path (used only when precompute_auc = FALSE) -----------------
.calculate_FC_tradeseq <- function(models, test_lineage, lineages, fc_names = NULL,
                                   genes, knots = NULL, N){
  dm <- colData(models)$tradeSeq$dm
  X  <- colData(models)$tradeSeq$X
  slingshotColData <- colData(models)$crv
  pseudotime <- slingshotColData[, grep(x = colnames(slingshotColData),
                                        pattern = "pseudotime")]
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  beta    <- betaMat[genes, ]
  knotPoints <- S4Vectors::metadata(models)$tradeSeq$knots
  lookup <- setNames(seq_along(lineages), lineages)
  time_points <- seq(from = 0, to = 1, length.out = N)

  lineage_id <- unname(lookup[test_lineage])
  df  <- tradeSeq:::.getPredictRangeDf(dm, lineage_id, nPoints = N)
  Xdf <- tradeSeq:::predictGAM(lpmatrix = X, df = df, pseudotime = pseudotime)
  yhat_mat_test   <- exp(Xdf %*% t(beta) + df$offset)
  auc_vector_test <- apply(yhat_mat_test, 2, function(y) pracma::trapz(time_points, y))
  if (!is.null(knots)) {
    t1 <- unname(knotPoints[knots[[1]]]); t2 <- unname(knotPoints[knots[[2]]])
    indices <- which(time_points >= t1 & time_points <= t2)
    yhat_mat_test   <- yhat_mat_test[indices, ]
    auc_vector_test <- apply(yhat_mat_test, 2, function(y) pracma::trapz(time_points[indices], y))
  }
  all_auc_log_diff <- c()
  for (lin in lineages) {
    if (lin != test_lineage) {
      lineage_id <- unname(lookup[lin])
      df  <- tradeSeq:::.getPredictRangeDf(dm, lineage_id, nPoints = N)
      Xdf <- tradeSeq:::predictGAM(lpmatrix = X, df = df, pseudotime = pseudotime)
      yhat_mat   <- exp(Xdf %*% t(beta) + df$offset)
      auc_vector <- apply(yhat_mat, 2, function(y) pracma::trapz(time_points, y))
      if (!is.null(knots)) {
        t1 <- unname(knotPoints[knots[[1]]]); t2 <- unname(knotPoints[knots[[2]]])
        indices <- which(time_points >= t1 & time_points <= t2)
        yhat_mat   <- yhat_mat[indices, ]
        auc_vector <- apply(yhat_mat, 2, function(y) pracma::trapz(time_points[indices], y))
      }
      all_auc_log_diff <- cbind(all_auc_log_diff, log2(auc_vector_test / auc_vector))
    }
  }
  if (length(lineages) == 2) {
    colnames(all_auc_log_diff) <- paste0("log2FC_", test_lineage, "vs", lineages[lineages != test_lineage])
  } else {
    colnames(all_auc_log_diff) <- paste0(fc_names, "_pattern")
  }
  all_auc_log_diff
}

.get_average_FC <- function(FC){
  linear_values <- 2^FC
  log2(mean(linear_values, na.rm = TRUE))
}

.get_meta_p <- function(Ps) {
  Ps <- Ps[!is.na(Ps)]
  Ps[Ps == 0] <- .Machine$double.xmin
  if (length(Ps) == 0) return(NA_real_)
  stat <- -2 * sum(log(Ps)); df <- 2 * length(Ps)
  pchisq(stat, df = df, lower.tail = FALSE)
}
