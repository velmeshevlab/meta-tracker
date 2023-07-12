monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(plot.title = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.line.y = element_line(size=1, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

#extend monocle3 class to add additional slots
#' @export
setClass("cell_data_set_ext", contains = "cell_data_set", slots=c(graphs = "list", lineages="list", expression="list", expectation="list", pseudotime="list")) -> cell_data_set_ext

#' @export
import_monocle <-function(cds){
cds <- as(cds,"cell_data_set_ext")
return(cds)
}

#' @export
#generate node plot
#filter = T will display only the nodes at branch points and at the ends of trajectories
#N controls the density of nodes to display if filter = T (larger values = less dense, N = 1 displays all nodes)
node_plot <- function(cds, filter = F, N = 50, size = 0.5){
Y <- cds@principal_graph_aux[["UMAP"]]$dp_mst
d = as.data.frame(t(Y))
if(filter == T){
g = principal_graph(cds)[["UMAP"]]
dd = degree(g)
names1 = names(dd[dd > 2 | dd == 1])
names2 = names(dd[dd == 2])
names2 = sample(names2, length(names2)/N, replace = F)
d.f = d[c(names1, names2),]
ggplot(data=d, aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.01) + geom_text_repel(data=d.f, aes(x=UMAP_1, y=UMAP_2), label=rownames(d.f), size=size, hjust = 2, color = "red", max.overlaps = Inf, segment.size = 0.1) + monocle_theme_opts()
}
else{
ggplot(data=d, aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.01) + geom_text(data=d, aes(x=UMAP_1, y=UMAP_2), label=rownames(d), size=size, hjust = 1, color = "red") + monocle_theme_opts()
}
}

#' @export
isolate_graph <- function(cds, start, end, lineage, include_nodes = NULL){
#get lineage graph
cds_name = deparse(substitute(cds))
sub.graph = isolate_graph_sub(cds, start, end, lineage, include_nodes = include_nodes)
input = paste0(cds_name, "@graphs$", lineage, " <- make_graph(sub.graph)")
eval(parse(text=input))
eval(parse(text=paste0("return(", cds_name, ")")))
}

#' @export
isolate_lineage <- function(cds, lineage, sel_clusters = NULL, start_regions = F, starting_clusters = F, subset = FALSE, N = 5, cl = 1){
cds_name = deparse(substitute(cds))
sel.cells = isolate_lineage_sub(cds, lineage, sel_clusters = sel_clusters, start_regions = start_regions, starting_clusters = starting_clusters, subset = subset, N = N, cl = cl)
input = paste0(cds_name, "@lineages$", lineage, " <- sel.cells")
eval(parse(text=input))
return(cds)
}

#' @export
combine_objects <- function(obj1, obj2, name1, name2){
  cds_new = new("cell_data_set_ext")
  #cds_new@'preprocess_aux'<-obj1@'preprocess_aux'
  cds_new@'reduce_dim_aux'<-obj1@'reduce_dim_aux'
  cds_new@'principal_graph_aux'<-obj1@'principal_graph_aux'
  cds_new@'principal_graph'<-obj1@'principal_graph'
  cds_new@'clusters'<-obj1@'clusters'
  cds_new@'int_elementMetadata'<-obj1@'int_elementMetadata'
  cds_new@'int_colData'<-obj1@'int_colData'
  cds_new@'int_metadata'<-obj1@'int_metadata'
  cds_new@'rowRanges'<-obj1@'rowRanges'
  cds_new@'colData'<-obj1@'colData'
  cds_new@'assays'<-obj1@'assays'
  cds_new@'NAMES'<-obj1@'NAMES'
  cds_new@'elementMetadata'<-obj1@'elementMetadata'
  cds_new@'metadata'<-obj1@'metadata'
  cds_new@'graphs'<-c(obj1@'graphs', obj2@'graphs')
  cds_new@'lineages'<-c(obj1@'lineages', obj2@'lineages')
  cds_new@'expression'<-c(obj1@'expression', obj2@'expression')
  cds_new@'expectation'<-c(obj1@'expectation', obj2@'expectation')
  cds_new@'pseudotime'<-c(obj1@'pseudotime', obj2@'pseudotime')
  names(cds_new@'graphs') <- c(paste(names(obj1@'graphs'), name1, sep = ""), paste(names(obj2@'graphs'), name2, sep = ""))
  names(cds_new@'lineages') <- c(paste(names(obj1@'lineages'), name1, sep = ""), paste(names(obj2@'lineages'), name2, sep = ""))
  names(cds_new@'expression') <- c(paste(names(obj1@'expression'), name1, sep = ""), paste(names(obj2@'expression'), name2, sep = ""))
  names(cds_new@'expectation') <- c(paste(names(obj1@'expectation'), name1, sep = ""), paste(names(obj2@'expectation'), name2, sep = ""))
  names(cds_new@'pseudotime') <- c(paste(names(obj1@'pseudotime'), name1, sep = ""), paste(names(obj2@'pseudotime'), name2, sep = ""))
  cds_new
  }

#' @export
combine_lineages <- function(cds, start){
  cds_name = deparse(substitute(cds))
  lineage = names(cds@lineages)[1]
  input = paste0(cds_name, "@graphs$", lineage)
  if(length(names(cds@lineages)) > 1){
    for(lineage in names(cds@lineages)[2:length(names(cds@lineages))]){
      input = paste0(input, ",", cds_name,"@graphs$", lineage)
    }
    input = paste0("igraph::union(", input, ")")
  }
  g = eval(parse(text=input))
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  principal_graph(cds)[["UMAP"]] <- g
  cds@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(g))]
  cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
  closest_vertex = apply(cells_UMAP[,c("UMAP_1", "UMAP_2")], 1, calculate_closest_vertex, nodes = as.matrix(nodes_UMAP[,names(V(g))]))
  closest_vertex = as.data.frame(closest_vertex)
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex
  source_url("https://raw.githubusercontent.com/cole-trapnell-lab/monocle3/master/R/learn_graph.R")
  cds <- project2MST(cds, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP[,names(V(g))])
  cds <- order_cells(cds, root_pr_nodes = as.character(paste0("Y_",start)))
  return(cds)
}

#' @export
path.distance <- function(path){
dists=c()
for(i in 2:nrow(path)){
x1 = path[i-1,1]
y1 = path[i-1,2]
x2 = path[i,1]
y2 = path[i,2]
d.x = x2 - x1
d.y = y2 - y1
dist = sqrt(d.x*d.x + d.y*d.y)
dists = append(dist,dists)
}
return(mean(dists))
}

cell.selector_sub2 <- function(cell, coords, r){
x2 = cell[1]
y2 = cell[2]
d.x = x2 - coords[1]
d.y = y2 - coords[2]
dist = sqrt(d.x*d.x + d.y*d.y)
if(dist <= r){
return(TRUE)
}
else{
return(FALSE)
}
}

#' @export
selector_sub <- function(node, cells, r){
x1 = node[1]
y1 = node[2]
res = apply(cells, 1, cell.selector_sub2, coords = c(x1, y1), r = r, simplify = T)
res = names(res[res == TRUE])
return(res)
}

#' @export
cell.selector <- function(path, cells, r, cl){
sel.cells = c()
sel.cells = pbapply(path, 1, selector_sub, cells = cells, r = r, cl = cl, simplify = T)
return(unique(unlist(sel.cells)))
}

#' @export
make_graph <- function(sub.graph){
edges = names(sub.graph)
start.edges = c()
end.edges = c()
for(i in 1:(length(edges)-1)){
start.edges = append(start.edges, edges[i])
end.edges = append(end.edges, edges[i+1])
}
d = cbind(start.edges, end.edges)
g = graph_from_data_frame(d, directed = F)
return(g)
}

#' @export
included <- function(graph, include_nodes){
all(include_nodes %in% names(graph))
}

isolate_graph_sub <- function(cds, start, end, lineage, include_nodes = NULL){
#get lineage graph
reduction_method = "UMAP"
graph = cds@principal_graph[[reduction_method]]
#select cells that are 1) progenitor cells from the region of interest (MGE, CGE) or 2) lineage-committed cells
sub.graph = all_simple_paths(graph, paste0("Y_", start), paste0("Y_", end))
if(length(include_nodes) > 0){
sub.graph = sub.graph[sapply(sub.graph, included, include_nodes = include_nodes)]
}
lengths = lengths(sub.graph)
#get the shortest path
n = which(lengths==min(lengths))[1]
sub.graph = sub.graph[[n]]
}

get_lineage_object <- function(cds, lineage = FALSE, start, N = FALSE)
{
cds_name = deparse(substitute(cds))
if(lineage != FALSE){
input = paste0("sub.graph = ",cds_name,"@graphs$", lineage)
eval(parse(text=input))
input = paste0("sel.cells = ",cds_name,"@lineages$", lineage)
eval(parse(text=input))
}
else{
sel.cells = colnames(cds)
}
sel.cells = sel.cells[sel.cells %in% colnames(cds)]
nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
if(N != FALSE){
if(N < length(sel.cells)){
sel.cells = sample(sel.cells, N)
}
}
#subset the moncole object
cds_subset = cds[,sel.cells]
#set the graph, node and cell UMAP coordinates
if(lineage == FALSE){
sub.graph = principal_graph(cds_subset)[["UMAP"]]
}
principal_graph(cds_subset)[["UMAP"]] <- sub.graph
cds_subset@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(sub.graph))]
cds_subset@clusters[["UMAP"]]$partitions <- cds_subset@clusters[["UMAP"]]$partitions[colnames(cds_subset)]
#recalculate closest vertex for the selected cells
cells_UMAP = as.data.frame(reducedDims(cds_subset)["UMAP"])
closest_vertex = apply(cells_UMAP[,c("UMAP_1", "UMAP_2")], 1, calculate_closest_vertex, nodes = as.matrix(nodes_UMAP[,names(V(sub.graph))]))
closest_vertex = as.data.frame(closest_vertex)
cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex
source_url("https://raw.githubusercontent.com/cole-trapnell-lab/monocle3/master/R/learn_graph.R")
cds_subset <- project2MST(cds_subset, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP[,names(V(sub.graph))])
cds_subset <- order_cells(cds_subset, root_pr_nodes = c(paste0("Y_", as.character(start))))
return(cds_subset)
}

isolate_lineage_sub <- function(cds, lineage, sel_clusters = NULL, start_regions = NULL, starting_clusters = NULL, subset = FALSE, N = 5, cl = 1){
  cds_name = deparse(substitute(cds))
  input = paste0("sub.graph = ",cds_name,"@graphs$", lineage)
  eval(parse(text=input))
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  if(subset == F){
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names(V(sub.graph))]))
  }
  else{
    g = principal_graph(cds)[["UMAP"]]
    dd = degree(g)
    names1 = names(dd[dd > 2 | dd == 1])
    names2 = names(dd[dd == 2])
    names2 = sample(names2, length(names2)/subset, replace = F)
    names = c(names1, names2)
    names = intersect(names(V(sub.graph)), names)
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names]))
  }
  #select cells along the graph
  mean.dist = path.distance(nodes_UMAP.sub)
  r = mean.dist*N
  cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
  cells_UMAP = cells_UMAP[,c("UMAP_1", "UMAP_2")]
  sel.cells = cell.selector(nodes_UMAP.sub, cells_UMAP, r, cl = cl)
  #only keep cells in the progenitor and lineage-specific clusters
  sel.cells1 = c()
  sel.cells2 = sel.cells
  if(length(starting_clusters) > 0){
    sel.cells1 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% starting_clusters])
  }
  if(length(start_regions) > 0){
    sel.cells1 = sel.cells1[sel.cells1 %in% rownames(cds@colData[cds@colData$region %in% start_regions,])]
  }
  if(length(sel_clusters) > 0){
    sel.cells2 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% sel_clusters])
  }
  cells = unique(c(sel.cells1, sel.cells2))
  sel.cells = sel.cells[sel.cells %in% cells]
  return(sel.cells)
}

#' @export
calculate_closest_vertex <- function(cells, nodes){
new.pos = as.numeric(cells)
nearest.idx <- which.min(colSums((nodes - new.pos)^2))
out = as.integer(gsub("Y_", "", names(nearest.idx)))
}

#' @export
connect_nodes <- function(cds, node1, node2, add_node = F){
graph.old = cds@principal_graph[["UMAP"]]
if(add_node == F){
graph.new <- add_edges(graph.old, c(node1, node2))
}
else{
node_coords = cds@principal_graph_aux[["UMAP"]]$dp_mst
node_X = (node_coords[1,node1] + node_coords[1,node2])/2
node_Y = (node_coords[2,node1] + node_coords[2,node2])/2
new_name = paste0("Y_", as.character(length(names(V(graph.old)))+1))
node_coords = as.data.frame(c(node_X, node_Y))
colnames(node_coords) = new_name
rownames(node_coords) = c("UMAP_1", "UMAP_2")
cds@principal_graph_aux[["UMAP"]]$dp_mst <- cbind(cds@principal_graph_aux[["UMAP"]]$dp_mst, node_coords)
graph.new <- add_vertices(graph.old, 1,attr = list(name = new_name))
graph.new <- add_edges(graph.new, c(node1, new_name))
graph.new <- add_edges(graph.new, c(new_name, node2))
}
cds@principal_graph[["UMAP"]] <- graph.new
return(cds)
}

compress_lineage_v3 <- function(cds, lineage, start, gene = FALSE, N = 500, cores = F){
  cds_name = deparse(substitute(cds))
  if(gene == FALSE){
    input = paste0("compress_expression_v3(",cds_name,", lineage = '", lineage, "', start = ", start, ", gene = ", gene, ", N = ", N, ", cores = ", cores, ")")
  }
  else{
    input = paste0("compress_expression_v3(",cds_name,", lineage = '", lineage, "', start = ", start, ", gene = '", gene, "', N = ", N, ", cores = ", cores, ")")
  }
  exp = eval(parse(text=input))
  input = paste0(cds_name, "@expression$", lineage, " <- exp$expression")
  eval(parse(text=input))
  input = paste0(cds_name, "@expectation$", lineage, " <- exp$expectation")
  eval(parse(text=input))
  input = paste0(cds_name, "@pseudotime$", lineage, " <- exp$pseudotime")
  eval(parse(text=input))
  eval(parse(text=paste0("return(",cds_name, ")")))
}

compress_expression_v3 <- function(cds, lineage, start, gene = FALSE, N = 500, cores = F){
  cds_name = deparse(substitute(cds))
  if(cores != F){
    cl <- makeCluster(cores)
    clusterEvalQ(cl, c(library(evobiR)))
  }
  input = paste0("get_lineage_object(",cds_name,", lineage = '", lineage, "', start = ", start, ")")
  cds_subset = eval(parse(text=input))
  family = stats::quasipoisson()
  model = "expression ~ splines::ns(pseudotime, df=3)"
  names(cds_subset) <- rowData(cds_subset)$gene_short_name
  exp = as.data.frame(as_matrix(exprs(cds_subset)))
  exp = t(t(exp) /  pData(cds_subset)[, 'Size_Factor'])
  pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt <- as.data.frame(pt)
  colnames(pt) <- c("pseudotime")
  pt$cell <- rownames(pt)
  rownames <- pt$cell
  pt$ID <- cds_subset@colData[rownames,]$sample
  exp = cbind(pt, t(exp))
  #use sliding window to compress pseudotime with ID information
  #add the mean values of pt if length of the ID is smaller than n
  exp = exp[order(exp$ID, exp$pseudotime),]
  pt = exp[,"pseudotime"]
  length <- length(unique(exp$ID))
  UMAP <- as.data.frame(subset2@reductions$umap@cell.embeddings)
  UMAP <- UMAP[rownames(exp),]
  pt.comp <- vector(mode = "list", length = length)
  ID.comp <- vector(mode = "list", length = length)
  UMAP.comp.x <- vector(mode = "list", length = length)
  UMAP.comp.y <- vector(mode = "list", length = length)
  len <- c()
  unique <- unique(exp$ID)
  ID <- exp$ID
  if (length > N)
  {
    return(FALSE)
  }
  n <- round(N/length)
  abs <- (N-n*length)
  if (abs < 0){
    n <- c(rep(n,length-abs(abs)), rep(n-1,abs(abs)))
    n <- sample(n, size = length, replace = FALSE)
  }
  if (abs > 0){
    n <- c(rep(n,length-abs(abs)), rep(n+1,abs(abs)))
    n <- sample(n, size = length, replace = FALSE)
  }
  for (i in 1:(length-1)){
    pos <- match(unique[i], ID)
    pos_1 <- match(unique[i+1], ID)
    n_1 <- n[i]
    window <- (pos_1-pos)/(n_1+1)
    #step <- (pos_1-pos-window)/n_1
    step <- window
    len <- c(len, (pos_1-pos))
    pt.comp_1 = SlidingWindow("mean", pt[pos:(pos_1-1)], window, step)
    UMAP.comp.x_1 = SlidingWindow("mean", UMAP$UMAP_1[pos:(pos_1-1)], window, step)
    UMAP.comp.y_1 = SlidingWindow("mean", UMAP$UMAP_2[pos:(pos_1-1)], window, step)
    ID.comp_1 <- as.character(rep(unique[i], times=length(pt.comp_1)))
    if (length(pt[pos:(pos_1-1)]) <= n_1){
      pt.comp_1 = pt[pos:(pos_1-1)]
      UMAP.comp.x_1 = UMAP$UMAP_1[pos:(pos_1-1)]
      UMAP.comp.y_1 = UMAP$UMAP_2[pos:(pos_1-1)]
      ID.comp_1 <- as.character(rep(unique[i], times=length(pt.comp_1)))
    }
    pt.comp[[i]] <- pt.comp_1
    UMAP.comp.x[[i]] = UMAP.comp.x_1
    UMAP.comp.y[[i]] = UMAP.comp.y_1
    ID.comp[[i]] <- ID.comp_1
    if (i == (length-1)){
      pos <- match(unique[i+1], ID)
      pos_1 <- length(ID)
      n_1 <- n[(length(n))]
      window <- (pos_1-pos+1)/(n_1+1)
      #step = ((pos_1-pos+1-window)/n_1)
      step <- window
      pt.comp_1 = SlidingWindow("mean", pt[pos:pos_1], window, step)
      UMAP.comp.x_1 = SlidingWindow("mean", UMAP$UMAP_1[pos:(pos_1)], window, step)
      UMAP.comp.y_1 = SlidingWindow("mean", UMAP$UMAP_2[pos:(pos_1)], window, step)
      ID.comp_1 <- as.character(rep(unique[i+1], times=length(pt.comp_1)))
      if (length(pt[pos:pos_1]) <= n_1){
        pt.comp_1 = pt[pos:pos_1]
        UMAP.comp.x_1 = UMAP$UMAP_1[pos:(pos_1)]
        UMAP.comp.y_1 = UMAP$UMAP_2[pos:(pos_1)]
        ID.comp_1 <- as.character(rep(unique[i+1], times=length(pt.comp_1)))
      }
      pt.comp[[i+1]] <- pt.comp_1
      UMAP.comp.x[[i+1]] = UMAP.comp.x_1
      UMAP.comp.y[[i+1]] = UMAP.comp.y_1
      ID.comp[[i+1]] <- ID.comp_1
    }
  }
  eta <- N-sum(lengths(pt.comp))
  if(eta > 0){
    max <- which.max(len)
    pos <- match(unique[max], ID)
    pos_1 <- match(unique[max+1], ID)
    window_new <- (pos_1-pos)/(n[max]+eta+1)
    #step_new <- (pos_1-pos-window_new)/(n[max]+eta)
    step_new <- window_new
    pt.comp[[max]] <- SlidingWindow("mean", pt[pos:(pos_1-1)], window_new, step_new)
    UMAP.comp.x[[max]] <- SlidingWindow("mean", UMAP$UMAP_1[pos:(pos_1-1)], window_new, step_new)
    UMAP.comp.y[[max]] <- SlidingWindow("mean", UMAP$UMAP_2[pos:(pos_1-1)], window_new, step_new)
    ID.comp[[max]] <- as.character(rep(unique[max], times=length(pt.comp[[max]])))
  }
  pt.comp <- unlist(pt.comp)
  UMAP.comp.x <- unlist(UMAP.comp.x)
  UMAP.comp.y <- unlist(UMAP.comp.y)
  ID.comp <- unlist(ID.comp)
  max.pt = max(pt.comp)
  
  #use sliding window to compress expression with ID information
  len <- c()
  exp.comp <- vector(mode = "list", length = length)
  if(gene != F){
    exp_test = exp[,c("pseudotime", "ID", gene)]
    for (i in 1:(length-1)){
      pos <- match(unique[i], ID)
      pos_1 <- match(unique[i+1], ID)
      n_1 <- n[i]
      window <- (pos_1-pos)/(n_1+1)
      #step <- (pos_1-pos-window)/n_1
      step <- window
      len <- c(len, (pos_1-pos))
      exp.comp_1 = SlidingWindow("mean", exp_test$gene[pos:(pos_1-1)], window, step)
      if (length(exp_test$gene[pos:(pos_1-1)]) <= n_1){
        exp.comp_1 = exp_test$gene[pos:(pos_1-1)]
      }
      exp.comp[[i]] <- exp.comp_1
      if (i == (length-1)){
        pos <- match(unique[i+1], ID)
        pos_1 <- length(ID)
        n_1 <- n[(length(n))]
        window <- (pos_1-pos+1)/(n_1+1)
        #step = ((pos_1-pos+1-window)/n_1)
        step <- window
        exp.comp_1 = SlidingWindow("mean", exp_test$gene[pos:pos_1], window, step)
        if (length(exp_test$gene[pos:pos_1]) <= n_1){
          exp.comp_1 = exp_test$gene[pos:pos_1]
        }
        exp.comp[[i+1]] <- exp.comp_1
      }
    }
    eta <- N-sum(lengths(exp.comp))
    if(eta > 0){
      max <- which.max(len)
      pos <- match(unique[max], ID)
      pos_1 <- match(unique[max+1], ID)
      window_new <- (pos_1-pos)/(n[max]+eta+1)
      #step_new <- (pos_1-pos-window_new)/(n[max]+eta)
      step_new <- window_new
      exp.comp[[max]] <- SlidingWindow("mean", exp_test$gene[pos:(pos_1-1)], window_new, step_new)
    }
    exp.comp <- unlist(exp.comp)
    
  }
  else{
    print(paste0("Compressing lineage ", lineage, " and fitting curves"))
    mat <- as.data.frame(exp[,4:ncol(exp)])
    exp.comp = pbapply(mat, 2, compress_3, length = length, unique = unique, ID = ID, n = n, N = 500, l = 108, cl = cl)
  }
  
  if(gene != F){
    exp_data.sel = cbind(pt.comp, exp.comp)
    exp_data.sel = as.data.frame(exp_data.sel)
    colnames(exp_data.sel) <- c("pseudotime", "expression")
    exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
    exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
    d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    colnames(d) <- c("pseudotime")
    tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
    fit = predict(fit, newdata = d, type='response')
    }, error=function(cond) {fit = as.data.frame(rep(0, N))})
    exp = as.data.frame(cbind(exp.comp, fit, d))
    colnames(exp) <- c("expression", "expectation", "pseudotime")
    exp$expression[exp$expression < 0] <- 0
    exp$expectation[exp$expectation < 0] <- 0
  }
  else{
    exp_data <- cbind(pt.comp, ID.comp, UMAP.comp.x, UMAP.comp.y, exp.comp)
    exp_data <- as.data.frame(exp_data)
    exp_data$pt.comp <- as.numeric(exp_data$pt.comp)
    exp_data_ordered <- exp_data[order(exp_data$pt.comp), ]
    mat <- exp_data_ordered[,5:(ncol(exp_data_ordered))]
    d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    if(cores != F){
      fit = pbsapply(mat, fit.m3, pt = d, max.pt = max(d), N = N, cl = cl)
    }
    else{
      fit = pbsapply(mat, fit.m3, pt = d, max.pt = max(d), N = N)
    }
    fit = apply(fit, 2, as.numeric)
    return(list("expression" = exp_data_ordered, "expectation" = fit, "pseudotime" = d))
  }
  exp$expression[exp$expression < 0] <- 0
  exp$expectation[exp$expectation < 0] <- 0
  if(cores != F){
    stopCluster(cl)
  }
  return(exp)
}

compress_3 <- function(df, length, unique, ID, n, N=500, l){
  len <- c()
  exp.comp <- vector(mode = "list", length = l)
  for (i in 1:(length-1)){
    pos <- match(unique[i], ID)
    pos_1 <- match(unique[i+1], ID)
    n_1 <- n[i]
    window <- (pos_1-pos)/n_1
    step <- (pos_1-pos-window)/n_1
    len <- c(len, (pos_1-pos))
    exp.comp_1 = SlidingWindow("mean", df[pos:(pos_1-1)], window, step)
    if (length(df[pos:(pos_1-1)]) <= n_1){
      exp.comp_1 = df[pos:(pos_1-1)]
    }
    exp.comp[[i]] <- exp.comp_1
    if (i == (length-1)){
      pos <- match(unique[i+1], ID)
      pos_1 <- length(ID)
      n_1 <- n[(length(n))]
      window <- (pos_1-pos+1)/n_1
      step = ((pos_1-pos+1-window)/n_1)
      len <- c(len, (pos_1-pos))
      exp.comp_1 = SlidingWindow("mean", df[pos:pos_1], window, step)
      if (length(df[pos:pos_1]) <= n_1){
        exp.comp_1 = df[pos:pos_1]
      }
      exp.comp[[i+1]] <- exp.comp_1
    }
  }
  eta <- N-sum(lengths(exp.comp))
  if(eta > 0){
    max <- which.max(len)
    pos <- match(unique[max], ID)
    pos_1 <- match(unique[max+1], ID)
    window_new <- (pos_1-pos)/(n[max]+eta)
    step_new <- (pos_1-pos-window_new)/(n[max]+eta)
    exp.comp[[max]] <- SlidingWindow("mean", df[pos:(pos_1-1)], window_new, step_new)
  }
  exp.comp <- unlist(exp.comp)
}

fit.m3 <- function(exp.sel, pt, max.pt, model = "expression ~ splines::ns(pseudotime, df=3)", N = 500){
  require(speedglm)
  family = stats::quasipoisson()
  exp_data.sel = cbind(pt, exp.sel)
  colnames(exp_data.sel) <- c("pseudotime","expression")
  exp_data.sel = as.data.frame(exp_data.sel)
  exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
  exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
  tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
  d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
  colnames(d) <- c("pseudotime")
  fit = stats::predict(fit, newdata=d, type="response")
  return(fit)
  }, error=function(cond) {return(rep("NA", N))})
}
