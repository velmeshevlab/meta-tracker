#' Compile metatracker lineage- and branch-specific gene results into two tables
#'
#' Flattens the per-lineage filtered gene tables and the per-branch gene tables
#' into two tidy data.frames.
#'
#' Lineage genes: for every lineage, the kept columns of \code{pattern_filtered}
#' and \code{diffend_filtered} are merged by gene (full outer join), then all
#' lineages are row-bound with the lineage name in the last column.
#'
#' Branch genes: the per-branch tables in \code{metadata(cds)$branch_genes} are
#' row-bound as-is (gene moved to its own column).
#'
#' @param cds A metatracker cds with \code{@lineage_genes[[lin]]$filtered$*} and
#'   \code{metadata(cds)$branch_genes} populated.
#' @param pattern_cols,diffend_cols Columns to keep from each filtered table.
#' @return list(lineage_genes = <data.frame>, branch_genes = <data.frame>).
compile_metatracker_results <- function(cds,
                                        pattern_cols = c("pvalue_combined_pattern",
                                                         "average_FC_pattern",
                                                         "median_FC_pattern"),
                                        diffend_cols = c("pvalue_combined_diffend",
                                                         "average_FC_diffend",
                                                         "median_FC_diffend")) {

  .grab <- function(df, cols) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df   <- as.data.frame(df)
    cols <- intersect(cols, colnames(df))
    out  <- df[, cols, drop = FALSE]
    out$gene <- rownames(df)
    out
  }

  # ---- lineage-specific genes ----
  lin_list <- list()
  for (lin in names(cds@lineage_genes)) {
    filt <- cds@lineage_genes[[lin]][["filtered"]]
    if (is.null(filt)) next
    p <- .grab(filt[["pattern_filtered"]], pattern_cols)
    d <- .grab(filt[["diffend_filtered"]], diffend_cols)
    if (is.null(p) && is.null(d)) next
    merged <- if (is.null(d)) p
              else if (is.null(p)) d
              else merge(p, d, by = "gene", all = TRUE)
    merged$lineage <- lin
    lin_list[[lin]] <- merged
  }
  lineage_genes <- .rbind_fill(lin_list)
  if (nrow(lineage_genes) > 0) {
    ord <- c("gene",
             intersect(c(pattern_cols, diffend_cols), colnames(lineage_genes)),
             "lineage")
    lineage_genes <- lineage_genes[, ord, drop = FALSE]
    rownames(lineage_genes) <- NULL
  }

  # ---- branch-specific genes ----
  bg <- metadata(cds)$branch_genes
  br_list <- list()
  for (nm in names(bg)) {
    df <- bg[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    df <- as.data.frame(df)
    df$gene <- rownames(df)
    br_list[[nm]] <- df
  }
  branch_genes <- .rbind_fill(br_list)
  if (nrow(branch_genes) > 0) {
    front <- intersect("gene", colnames(branch_genes))
    branch_genes <- branch_genes[, c(front, setdiff(colnames(branch_genes), front)),
                                 drop = FALSE]
    rownames(branch_genes) <- NULL
  }

  list(lineage_genes = lineage_genes, branch_genes = branch_genes)
}

# rbind a list of data.frames that may not share all columns: union the columns,
# fill missing ones with NA, then rbind. Returns an empty data.frame if the list
# is empty.
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
