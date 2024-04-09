# meta-tracker basic tutorial
This is the analysis guide starting from a monocle 3 object to identify lineage-specific genes and compare conditions within the same lineage.
## Step-by-step tutorial
1. Import the monocle 3 object.

```cds = import_monocle(cds)```

3. Generate a node plot to visualize node names for downstream analysis.

node_plot(cds)

3. Isolate cells along specific lineages. You need to specify the metatracker object, start and end of the trajectory as integer node numbers, name you want to assign to the trajectory, as well as an optional parameter of nodes to include in the trajectory (in case you want your trajectory to pass through specific points).

#isolate lineages
#L2-3
lineage = "L2_3"
start = 1334
end = 2344
inc.node = c("Y_1853", "Y_572", "Y_2177")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("1", "28", "0", "9", "13", "8", "2", "27", "29", "39", "5")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 5)
