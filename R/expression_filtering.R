# Expression-based gene filtering.

# Normalise a @lineages entry to a character vector of cell barcodes,
# whether it is stored as a plain vector or as a data.frame/list with a
# `name` column.
.lineage_cells <- function(x) {
  if (is.data.frame(x) || is.list(x)) {
    if (!is.null(x[["name"]])) return(as.character(x[["name"]]))
    return(as.character(unlist(x, use.names = FALSE)))
  }
  as.character(x)
}

#' Genes expressed above a threshold in any lineage
#'
#' @param cds A \code{metatracker_data_set} with populated \code{@lineages}.
#' @param mode "number" (absolute cell count) or "ratio" (fraction of lineage cells).
#' @param N Cell-count cutoff when \code{mode = "number"}.
#' @param ratio Fraction cutoff when \code{mode = "ratio"}.
#' @return Character vector of gene names.
#' @export
filter_by_expression <- function(cds, mode = "number", N = 100, ratio = 0.01){
  data = counts(cds)
  lineages = names(cds@lineages)
  all.expressed_genes = c()
  for(lineage in lineages){
    cells = .lineage_cells(cds@lineages[[lineage]])
    missing = setdiff(cells, colnames(data))
    if (length(missing) > 0)
      warning(sprintf("Lineage '%s': %d/%d cells not in cds; dropping them.",
                      lineage, length(missing), length(cells)), call. = FALSE)
    cells = cells[cells %in% colnames(data)]
    if (length(cells) == 0) next
    data.sub = data[, cells, drop = FALSE]
    if(mode == "ratio"){
      cutoff = ratio*ncol(data.sub)
    } else {
      cutoff = N
    }
    expressed_genes = rownames(data.sub)[rowSums(data.sub > 0) > cutoff]
    all.expressed_genes = append(all.expressed_genes, expressed_genes)
  }
  all.expressed_genes = unique(all.expressed_genes)
  all.expressed_genes
}

.filter_by_expression_lineage <- function(cds, lineage, mode = "number", N = 100, ratio = 0.01){
  data = counts(cds)
  cells = .lineage_cells(cds@lineages[[lineage]])
  cells = cells[cells %in% colnames(data)]
  if (length(cells) == 0) return(character(0))
  data.sub = data[, cells, drop = FALSE]
  if(mode == "ratio"){
    cutoff = ratio*ncol(data.sub)
  } else {
    cutoff = N
  }
  expressed_genes = rownames(data.sub)[rowSums(data.sub > 0) > cutoff]
  expressed_genes
}
