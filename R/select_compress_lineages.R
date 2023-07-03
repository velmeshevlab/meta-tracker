setClass("cell_data_set_ext", contains = "cell_data_set", slots=c(graphs = "list", lineages="list", expression="list", expectation="list", pseudotime="list")) -> cell_data_set_ext

#' @export
import_monocle <-function(cds){
cds <- as(cds,"cell_data_set_ext")
return(cds)
}
