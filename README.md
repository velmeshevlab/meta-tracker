This is the analysis guide starting from cleaned up trajectories to identify lineage-specific genes and compare conditions within the same lineage.


#isolate lineages
#L2-3
lineage = "L2_3"
start = 1334
end = 2344
inc.node = c("Y_1853", "Y_572", "Y_2177")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("1", "28", "0", "9", "13", "8", "2", "27", "29", "39", "5")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 5)
