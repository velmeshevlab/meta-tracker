#' Compile metatracker results into pattern, diffend, and branch tables
#'
#' Produces three tidy data.frames:
#'   * lineage_pattern - every lineage's \code{filtered$pattern_filtered}, kept
#'     columns + gene + lineage;
#'   * lineage_diffend - every lineage's \code{filtered$diffend_filtered}, same;
#'   * branch_genes    - every \code{metadata(cds)$branch_genes} table, rbound.
#'
#' Within each lineage (and each branch) the rows are sorted by p-value
#' (ascending) then fold change (descending) BEFORE the lineages/branches are
#' concatenated, so each lineage/branch block stays internally ranked.
#'
#' @param cds A metatracker cds with \code{@lineage_genes[[lin]]$filtered$*} and
#'   \code{metadata(cds)$branch_genes} populated.
#' @param pattern_cols,diffend_cols Columns to keep from each filtered table.
#' @param sort_fc Which fold-change column to sort by, "median" or "average".
#' @param branch_p_col,branch_fc_col p-value / fold-change columns in the branch
#'   tables used for sorting (defaults match format_branch_specific_genes()).
#' @return list(lineage_pattern, lineage_diffend, branch_genes).
compile_metatracker_results <- function(cds,
                                        pattern_cols = c("pvalue_combined_pattern",
                                                         "average_FC_pattern",
                                                         "median_FC_pattern"),
                                        diffend_cols = c("pvalue_combined_diffend",
                                                         "average_FC_diffend",
                                                         "median_FC_diffend"),
                                        sort_fc = c("median", "average"),
                                        branch_p_col  = "meta_p",
                                        branch_fc_col = "median_FC") {
  sort_fc  <- match.arg(sort_fc)
  patt_pcol <- "pvalue_combined_pattern"; patt_fcol <- paste0(sort_fc, "_FC_pattern")
  diff_pcol <- "pvalue_combined_diffend"; diff_fcol <- paste0(sort_fc, "_FC_diffend")

  .grab <- function(df, cols, lineage) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df   <- as.data.frame(df)
    keep <- intersect(cols, colnames(df))
    out  <- df[, keep, drop = FALSE]
    out$gene    <- rownames(df)
    out$lineage <- lineage
    out
  }
  # sort by p ascending, fold change descending; rows with NA keys go last
  .sort_df <- function(df, p_col, fc_col) {
    if (is.null(df) || nrow(df) == 0) return(df)
    p  <- if (p_col  %in% colnames(df)) df[[p_col]]  else rep(NA_real_, nrow(df))
    fc <- if (fc_col %in% colnames(df)) df[[fc_col]] else rep(NA_real_, nrow(df))
    df[order(p, -fc, na.last = TRUE), , drop = FALSE]
  }
  .finish <- function(df, ord) {
    if (nrow(df) == 0) return(df)
    df <- df[, intersect(ord, colnames(df)), drop = FALSE]
    rownames(df) <- NULL
    df
  }

  # ---- lineage pattern / diffend, sorted per lineage then stacked ----
  patt_list <- list(); diff_list <- list()
  for (lin in names(cds@lineage_genes)) {
    filt <- cds@lineage_genes[[lin]][["filtered"]]
    if (is.null(filt)) next
    p <- .grab(filt[["pattern_filtered"]], pattern_cols, lin)
    d <- .grab(filt[["diffend_filtered"]], diffend_cols, lin)
    if (!is.null(p)) patt_list[[lin]] <- .sort_df(p, patt_pcol, patt_fcol)
    if (!is.null(d)) diff_list[[lin]] <- .sort_df(d, diff_pcol, diff_fcol)
  }
  lineage_pattern <- .finish(.rbind_fill(patt_list), c("gene", pattern_cols, "lineage"))
  lineage_diffend <- .finish(.rbind_fill(diff_list), c("gene", diffend_cols, "lineage"))

  # ---- branch genes, sorted per branch then stacked ----
  bg <- metadata(cds)$branch_genes
  br_list <- list()
  for (nm in names(bg)) {
    df <- bg[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    df <- as.data.frame(df)
    df$gene <- rownames(df)
    br_list[[nm]] <- .sort_df(df, branch_p_col, branch_fc_col)
  }
  branch_genes <- .rbind_fill(br_list)
  if (nrow(branch_genes) > 0) {
    front <- intersect("gene", colnames(branch_genes))
    branch_genes <- branch_genes[, c(front, setdiff(colnames(branch_genes), front)),
                                 drop = FALSE]
    rownames(branch_genes) <- NULL
  }

  list(lineage_pattern = lineage_pattern,
       lineage_diffend = lineage_diffend,
       branch_genes    = branch_genes)
}

# rbind a list of data.frames that may not share all columns: union the columns,
# fill missing ones with NA, then rbind. Empty data.frame if the list is empty.
.rbind_fill <- function(dfs) {
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (length(dfs) == 0) return(data.frame())
  all_cols <- unique(unlist(lapply(dfs, colnames)))
  dfs <- lapply(dfs, function(df) {
    for (m in setdiff(all_cols, colnames(df))) df[[m]] <- NA
    df[, all_cols, drop = FALSE]
  })
  do.call(rbind, dfs)
}
