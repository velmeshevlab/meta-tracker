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
