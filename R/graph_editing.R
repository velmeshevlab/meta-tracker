# Import, interactive principal-graph editing, and lineage isolation.

#' @export
import_monocle <-function(cds){
cds <- as(cds,"metatracker_data_set")
return(cds)
}

# Interactive graph modification with gene expression coloring
#' @export
graph_mod_interactive <- function(cds,
                                  reduction_method = "UMAP",
                                  point_size       = 1,
                                  node_size        = 3,
                                  segment_size     = 1,
                                  min_expr=0.1,
                                  N                = 10) {
  g        <- cds@principal_graph[[reduction_method]]
  Y        <- cds@principal_graph_aux[[reduction_method]]$dp_mst
  nodes_df <- as.data.frame(t(Y), stringsAsFactors = FALSE)
  colnames(nodes_df) <- c("x", "y")
  nodes_df$node <- rownames(nodes_df)
  nodes_df$col  <- "black"
  
  cell_coords <- as.data.frame(reducedDims(cds)[[reduction_method]], stringsAsFactors = FALSE)
  colnames(cell_coords) <- c("x", "y")
  cell_coords$cell_id <- rownames(cell_coords)
  metadata <- as.data.frame(colData(cds), stringsAsFactors = FALSE)
  
  set.seed(42)
  total_cells <- nrow(cell_coords)
  sample_size <- ceiling(total_cells / N)
  sampled_cells <- cell_coords %>% slice(sample(seq_len(total_cells), sample_size))
  
  cds_exprs_all <- as.matrix(SingleCellExperiment::counts(cds)[ , sampled_cells$cell_id, drop = FALSE])
  cds_exprs_all <- t(t(cds_exprs_all) / size_factors(cds)[sampled_cells$cell_id])
  
  el <- as.data.frame(get.edgelist(g), stringsAsFactors = FALSE)
  colnames(el) <- c("from", "to")
  edges_df <- data.frame(
    col = "black",
    from = el$from,
    to   = el$to,
    x    = nodes_df$x[match(el$from, nodes_df$node)],
    y    = nodes_df$y[match(el$from, nodes_df$node)],
    xend = nodes_df$x[match(el$to,   nodes_df$node)],
    yend = nodes_df$y[match(el$to,   nodes_df$node)],
    stringsAsFactors = FALSE
  )
  
  ui <- fluidPage(
    titlePanel("Interactive Graph Modifier"),
    fluidRow(
      column(
        width = 12, align = "center",
        selectInput("metadata", "Metadata:", choices = colnames(metadata)),
        selectInput("gene", "Gene:", choices = c("", rownames(cds_exprs_all))),
        actionButton("zoom_in",    label = NULL, icon = icon("search-plus")),
        actionButton("zoom_out",   label = NULL, icon = icon("search-minus")),
        actionButton("reset_zoom", "Reset Zoom"),
        actionButton("cancel",     label = NULL, icon = icon("undo")),
        tags$span("Pan:"),
        actionButton("pan_left",   label = NULL, icon = icon("arrow-left")),
        actionButton("pan_up",     label = NULL, icon = icon("arrow-up")),
        actionButton("pan_down",   label = NULL, icon = icon("arrow-down")),
        actionButton("pan_right",  label = NULL, icon = icon("arrow-right")),
        actionButton("done",       "Done")
      )
    ),
    fluidRow(column(width = 12, plotOutput("plot", click = "clickposition", hover = "hoverpos", width = "100%", height = "600px"))),
    fluidRow(column(width = 12, verbatimTextOutput("status")))
  )
  
  server <- function(input, output, session) {
    rv <- reactiveValues(
      nodes   = nodes_df,
      edges   = edges_df,
      stage   = 0,
      new_pt  = NULL,
      hovered = NULL,
      selected_node = NULL,
      xlim    = range(nodes_df$x),
      ylim    = range(nodes_df$y)
    )
    
    observeEvent(input$cancel, {
      if (rv$stage == 1 && !is.null(rv$new_pt)) {
        rv$new_pt <- NULL
        rv$stage <- 0
      } else if (rv$stage == 2 && !is.null(rv$selected_node)) {
        rv$nodes$col[rv$nodes$node == rv$selected_node] <- "black"
        rv$selected_node <- NULL
        rv$stage <- 0
      } else if (nrow(rv$edges) > nrow(edges_df)) {
        last_edge <- tail(rv$edges, 1)
        rv$edges <- rv$edges[-nrow(rv$edges), ]
        new_node <- last_edge$from
        if (!new_node %in% c(edges_df$from, edges_df$to)) {
          rv$nodes <- rv$nodes[rv$nodes$node != new_node, ]
        }
        rv$stage <- 0
        rv$new_pt <- NULL
        rv$selected_node <- NULL
      }
    })
    
    observeEvent(input$done, {
      stopApp(list(nodes = rv$nodes, edges = rv$edges))
    })
    
    output$status <- renderText({
      if (rv$stage == 0) {
        "Step 1: Click to add node or select existing node to connect."
      } else if (rv$stage == 1) {
        "Step 2: Click an existing node to connect to new node."
      } else if (rv$stage == 2) {
        paste("Step 2: Click second node to connect with", rv$selected_node)
      }
    })
    
    observeEvent(input$clickposition, {
      pt <- input$clickposition
      # If a node is already highlighted (green), connect it to clicked node
      if (!is.null(rv$hovered) && is.null(rv$new_pt)) {
        if (rv$stage == 0) {
          # Select the first node to connect from
          rv$selected_node <- rv$hovered
          rv$nodes$col[rv$nodes$node == rv$selected_node] <- "red"
          rv$hovered <- NULL
          rv$stage <- 2
        } else if (rv$stage == 2) {
          # Connect selected node to second hovered node
          from <- rv$selected_node
          to <- rv$hovered
          if (from != to) {
            rv$edges <- rbind(rv$edges, data.frame(
              from = from, to = to,
              x = rv$nodes$x[rv$nodes$node == from],
              y = rv$nodes$y[rv$nodes$node == from],
              xend = rv$nodes$x[rv$nodes$node == to],
              yend = rv$nodes$y[rv$nodes$node == to],
              col = "cyan",
              stringsAsFactors = FALSE
            ))
          }
          rv$nodes$col[rv$nodes$node == rv$selected_node] <- "black"
          rv$stage <- 0
          rv$selected_node <- NULL
          rv$new_pt <- NULL
        }
      } else {
        # If no node is hovered, add a new node and prepare to connect
        if (rv$stage == 0) {
          rv$new_pt  <- c(pt$x, pt$y)
          rv$stage   <- 1
        } else if (rv$stage == 1) {
          d2 <- (rv$nodes$x - pt$x)^2 + (rv$nodes$y - pt$y)^2
          targ <- rv$nodes$node[which.min(d2)]
          max_node = max(as.numeric(str_split_i(rv$nodes$node, "_", 2)))
          new_name <- paste0("Y_", max_node + 1)
          rv$nodes <- rbind(rv$nodes, data.frame(x = rv$new_pt[1], y = rv$new_pt[2], node = new_name, col = "cyan", stringsAsFactors = FALSE))
          rv$edges <- rbind(rv$edges, data.frame(from = new_name, to = targ,
                                                 x = rv$new_pt[1], y = rv$new_pt[2],
                                                 xend = rv$nodes$x[rv$nodes$node == targ],
                                                 yend = rv$nodes$y[rv$nodes$node == targ],
                                                 col = "cyan",
                                                 stringsAsFactors = FALSE))
          rv$stage   <- 0
          rv$new_pt  <- NULL
        }
      }
    })
    
    observeEvent(input$hoverpos, {
      hv <- input$hoverpos
      d2 <- (rv$nodes$x - hv$x)^2 + (rv$nodes$y - hv$y)^2
      i  <- which.min(d2)
      if (!is.null(rv$selected_node) && rv$nodes$node[i] == rv$selected_node) {
        rv$hovered <- NULL  # prevent self-hover during connection
      } else if (sqrt(d2[i]) < diff(rv$xlim)*0.01) {
        rv$hovered <- rv$nodes$node[i]
      } else if (sqrt(d2[i]) < diff(rv$xlim)*0.01) {
        rv$hovered <- rv$nodes$node[i]
      } else {
        rv$hovered <- NULL
      }
    })
    
    observeEvent(input$reset_zoom, { rv$xlim <- range(rv$nodes$x); rv$ylim <- range(rv$nodes$y) })
    observeEvent(input$zoom_in,   { dx <- diff(rv$xlim); dy <- diff(rv$ylim); rv$xlim <- rv$xlim + c(0.1 * dx, -0.1 * dx); rv$ylim <- rv$ylim + c(0.1 * dy, -0.1 * dy) })
    observeEvent(input$zoom_out,  { dx <- diff(rv$xlim); dy <- diff(rv$ylim); rv$xlim <- rv$xlim + c(-0.1 * dx, 0.1 * dx); rv$ylim <- rv$ylim + c(-0.1 * dy, 0.1 * dy) })
    observeEvent(input$pan_left,  { dx <- diff(rv$xlim); rv$xlim <- rv$xlim + c(-0.1 * dx, -0.1 * dx) })
    observeEvent(input$pan_right, { dx <- diff(rv$xlim); rv$xlim <- rv$xlim + c( 0.1 * dx,  0.1 * dx) })
    observeEvent(input$pan_up,    { dy <- diff(rv$ylim); rv$ylim <- rv$ylim + c( 0.1 * dy,  0.1 * dy) })
    observeEvent(input$pan_down,  { dy <- diff(rv$ylim); rv$ylim <- rv$ylim + c(-0.1 * dy, -0.1 * dy) })
    
    observe({
      gene_selected <- !is.null(input$gene) && input$gene != "" && input$gene %in% rownames(cds_exprs_all)
      if (gene_selected) {
        output$plot <- renderPlot({
          expr_values <- cds_exprs_all[input$gene, sampled_cells$cell_id]
          sampled <- sampled_cells %>% mutate(expr = expr_values)
          p <- ggplot() +
            geom_point(data = sampled, aes(x = x, y = y, color = log10(expr + min_expr)), size = point_size) +
            scale_color_gradient(low = adjustcolor("gray75", alpha.f = 0.2), high = "red")
          
          p <- p +
            geom_segment(data = rv$edges, aes(x = x, y = y, xend = xend, yend = yend), size = segment_size, color = rv$edges$col) +
            geom_point(data = rv$nodes, aes(x = x, y = y), size = node_size, color = rv$nodes$col)
          
          if (!is.null(rv$hovered)) {
            hrow <- rv$nodes[rv$nodes$node == rv$hovered, ]
            color_highlight <- if (!is.null(rv$selected_node) && rv$hovered == rv$selected_node) "red" else "green"
            p <- p + geom_point(data = hrow, aes(x = x, y = y), size = node_size * 1.5, color = color_highlight)
          }
          
          if (!is.null(rv$new_pt)) {
            temp_df <- data.frame(x = rv$new_pt[1], y = rv$new_pt[2])
            p <- p + geom_point(data = temp_df, aes(x = x, y = y), size = node_size * 1.5, color = "red")
          }
          
          p + coord_cartesian(xlim = rv$xlim, ylim = rv$ylim) + theme_minimal()
        })
      }
      else if (!gene_selected && !is.null(input$metadata) && input$metadata %in% colnames(metadata)) {
        output$plot <- renderPlot({
          meta_values <- metadata[sampled_cells$cell_id, input$metadata]
          sampled <- sampled_cells %>% mutate(meta_value = meta_values)
          
          is_discrete <- is.factor(meta_values) || is.character(meta_values)
          p <- ggplot() +
            geom_point(data = sampled, aes(x = x, y = y, color = meta_value), size = point_size) +
            if (is_discrete) {
              scale_color_discrete(name = input$metadata)
            } else {
              scale_color_gradient(low = "lightblue", high = "firebrick4", name = input$metadata, trans = "log10")
            }
          
          p <- p +
            geom_segment(data = rv$edges, aes(x = x, y = y, xend = xend, yend = yend), size = segment_size, color = rv$edges$col) +
            geom_point(data = rv$nodes, aes(x = x, y = y), size = node_size, color = rv$nodes$col)
          
          if (!is.null(rv$hovered)) {
            hrow <- rv$nodes[rv$nodes$node == rv$hovered, ]
            p <- p + geom_point(data = hrow, aes(x = x, y = y), size = node_size * 1.5, color = "green")
          }
          
          if (!is.null(rv$new_pt)) {
            temp_df <- data.frame(x = rv$new_pt[1], y = rv$new_pt[2])
            p <- p + geom_point(data = temp_df, aes(x = x, y = y), size = node_size * 1.5, color = "red")
          }
          
          p + coord_cartesian(xlim = rv$xlim, ylim = rv$ylim) + theme_minimal()
        })
      }
    })
    
    observeEvent(input$metadata, {
      updateSelectInput(session, "gene", selected = "")
    })
    
    observeEvent(input$done, { stopApp(list(nodes = rv$nodes, edges = rv$edges)) })
  }
  
  res <- runApp(shinyApp(ui, server))
  new_nodes <- res$nodes
  new_edges <- res$edges
  new_nodes = data.frame(name=new_nodes$node, x = new_nodes$x, y = new_nodes$y)
  rownames(new_nodes) <- new_nodes$name
  
  cds@principal_graph_aux[[reduction_method]]$dp_mst <- t(as.matrix(new_nodes[, c("x", "y")]))
  cds@principal_graph[[reduction_method]] <- graph_from_data_frame(new_edges[, c("from", "to")], vertices = new_nodes, directed = FALSE)
  
  return(cds)
  
}

# Depends on: shiny, igraph, dplyr, ggplot2, plotly, colorspace, SingleCellExperiment, Matrix
#' @export
graph_selection_interactive <- function(cds,
                                        lineage,
                                        point_size       = 0.5,
                                        node_size        = 1,
                                        reduction_method = "UMAP",
                                        segment_size     = 1,
                                        N                 = 10) {
  # Extract principal graph and node coordinates
  g <- cds@principal_graph[[reduction_method]]
  Y <- cds@principal_graph_aux[[reduction_method]]$dp_mst
  nodes <- as.data.frame(t(Y), stringsAsFactors = FALSE)
  colnames(nodes) <- c("x", "y")
  nodes$node <- rownames(nodes)
  
  # Cell embedding coordinates
  cell_coords <- as.data.frame(reducedDims(cds)[[reduction_method]], stringsAsFactors = FALSE)
  colnames(cell_coords) <- c("x", "y")
  cell_coords$cell_id <- rownames(cell_coords)
  metadata <- as.data.frame(colData(cds), stringsAsFactors = FALSE)
  
  # Sample cells
  set.seed(42)
  total_cells <- nrow(cell_coords)
  sample_size <- ceiling(total_cells / N)
  sampled_cells <- cell_coords %>%
    slice(sample(seq_len(total_cells), sample_size))
  
  # Precompute normalized expression for sampled cells
  cds_exprs_all <- as.matrix(counts(cds)[ , sampled_cells$cell_id, drop = FALSE])
  cds_exprs_all <- t(t(cds_exprs_all) / size_factors(cds)[sampled_cells$cell_id])
  
  # Build edges data frame
  el <- as.data.frame(get.edgelist(g), stringsAsFactors = FALSE)
  colnames(el) <- c("from", "to")
  coords <- nodes[, c("node","x","y")]
  edges <- data.frame(
    from = el$from,
    to   = el$to,
    x    = coords$x[match(el$from, coords$node)],
    y    = coords$y[match(el$from, coords$node)],
    xend = coords$x[match(el$to,   coords$node)],
    yend = coords$y[match(el$to,   coords$node)],
    stringsAsFactors = FALSE
  )
  
  # Helper: recompute selected path edges
  recalc <- function(clicks) {
    selected <- unique(clicks)
    highlight <- edges[0, , drop = FALSE]
    if (length(clicks) >= 2) {
      for (i in seq_len(length(clicks)-1)) {
        pth <- shortest_paths(
          g,
          from    = clicks[i],
          to      = clicks[i+1],
          weights = NA,
          output  = "vpath"
        )$vpath[[1]]
        seq_nodes <- names(pth)
        selected <- unique(c(selected, seq_nodes))
        pairs <- tibble(
          from = head(seq_nodes, -1),
          to   = tail(seq_nodes, -1)
        )
        segs <- edges %>% semi_join(pairs, by = c("from","to"))
        revs <- edges %>% semi_join(pairs, by = c("from"="to","to"="from"))
        highlight <- distinct(bind_rows(highlight, segs, revs))
      }
    }
    list(all_selected = selected, highlight_edges = highlight)
  }
  
  # Define Shiny app UI
  app <- shinyApp(
    ui = fluidPage(
      titlePanel("Interactive Graph Path Selector"),
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            inputId = "gene",
            label   = "Search gene:",
            choices = c("", rownames(counts(cds))),
            selected = character(0),
            multiple = FALSE,
            options = list(
              placeholder = 'Type a gene...',
              server = TRUE
            )
          ),
          selectInput(
            inputId = "color_by",
            label   = "Color cells by:",
            choices = colnames(metadata),
            selected = colnames(metadata)[1]
          ),
          actionButton("undo", "Undo Last Selection"),
          actionButton("done", "Finish & Return Selection"),
          br(), br(),
          verbatimTextOutput("sel_nodes")
        ),
        mainPanel(
          plotlyOutput("plot", height = "600px"),
          uiOutput(
            "hover_info",
            style = paste(
              "position:absolute; pointer-events:none;",
              "background: rgba(255,255,255,0.8); padding:4px;",
              "border:1px solid #ccc; border-radius:4px;"
            )
          )
        )
      )
    ),
    server = function(input, output, session) {
      rv <- reactiveValues(
        clicked         = character(),
        all_selected    = character(),
        highlight_edges = edges[0, , drop = FALSE]
      )
      # Clear gene selection when metadata changes
      observeEvent(input$color_by, {
        updateSelectizeInput(session, "gene", selected = character(0))
      })
      
      # Hover info
      output$hover_info <- renderUI({
        hov <- event_data("plotly_hover", source = "graph")
        if (is.null(hov)) return(NULL)
        key <- hov$key
        label <- if (key %in% nodes$node) key else metadata[key, input$color_by]
        div(
          style = sprintf("position:absolute; left:%dpx; top:%dpx;", hov$clientX+10, hov$clientY+10),
          strong(label)
        )
      })
      
      # Node click handling
      observeEvent(event_data("plotly_click", source = "graph"), {
        clk <- event_data("plotly_click", source = "graph")
        key <- clk$key
        if (!is.null(key) && key %in% nodes$node) {
          rv$clicked <- c(rv$clicked, key)
          rec <- recalc(rv$clicked)
          rv$all_selected    <- rec$all_selected
          rv$highlight_edges <- rec$highlight_edges
        }
      })
      
      # Undo last selection
      observeEvent(input$undo, {
        if (length(rv$clicked) > 0) {
          rv$clicked <- head(rv$clicked, -1)
          rec <- recalc(rv$clicked)
          rv$highlight_edges <- rec$highlight_edges
        }
      })
      
      # Plot rendering with conditional coloring
      output$plot <- renderPlotly({
        sampled <- sampled_cells
        if (nzchar(input$gene)) {
          expr_vals <- cds_exprs_all[input$gene, sampled$cell_id]
          sampled$expr_value <- expr_vals
          p <- ggplot() +
            geom_point(
              data = sampled,
              aes(x = x, y = y, key = cell_id,
                  text = sprintf("%s: %.3f", input$gene, expr_value),
                  color = expr_value),
              size  = point_size, alpha = 0.6
            ) +
            scale_colour_gradient(low = "grey75", high = "red",
                                  name = input$gene, limits = c(0, as.numeric(quantile(expr_vals, 0.99))))
        } else {
          sampled <- sampled_cells %>%
            mutate(meta_value = metadata[cell_id, input$color_by])
          vals <- unique(sampled$meta_value)
          pal <- qualitative_hcl(length(vals), palette = "Dark 3")
          names(pal) <- vals
          p <- ggplot() +
            geom_point(
              data = sampled,
              aes(x = x, y = y, key = cell_id,
                  text = meta_value, color = meta_value),
              size = point_size, alpha = 0.6
            ) +
            scale_color_manual(values = pal, na.value = "grey50",
                               name = input$color_by)
        }
        p <- p +
          geom_segment(
            data = edges,
            aes(x = x, y = y, xend = xend, yend = yend),
            size = segment_size, alpha = 0.5, color = "cyan"
          ) +
          {if (nrow(rv$highlight_edges) > 0) geom_segment(
            data = rv$highlight_edges,
            aes(x = x, y = y, xend = xend, yend = yend),
            size = segment_size, color = "red"
          )} +
          geom_point(
            data = nodes,
            aes(x = x, y = y, key = node, text = node),
            size = node_size, color = "black"
          ) +
          geom_point(
            data = subset(nodes, node %in% rv$clicked),
            aes(x = x, y = y, key = node, text = node),
            size = node_size*1.5, color = "orange"
          ) +
          theme_minimal()
        
        ggplotly(p, tooltip = "text", source = "graph") %>% layout(dragmode = "select")
      })
      
      # Selected nodes print
      output$sel_nodes <- renderPrint({
        if (length(rv$clicked) == 0) "No selections." else rv$all_selected
      })
      
      # Finish selection
      observeEvent(input$done, stopApp(list(selected_nodes = rv$all_selected)))
    }
  )
  
  # Run and return updated CDS with subgraph
  result <- runApp(app)
  sel <- result$selected_nodes
  subg <- induced_subgraph(g, vids = V(g)[name %in% sel])
  cds@graphs[[lineage]] <- subg
  return(cds)
}

#' @export
isolate_lineage <- function(cds, lineage, sel_clusters = NULL, start_regions = F, starting_clusters = F, subset = FALSE, N = 5, cl = 1){
sel.cells = .isolate_lineage_sub(cds, lineage, sel_clusters = sel_clusters, start_regions = start_regions, starting_clusters = starting_clusters, subset = subset, N = N, cl = cl)
cds@lineages[[lineage]] <- sel.cells
return(cds)
}

.isolate_lineage_sub <- function(cds, lineage, sel_clusters = NULL, start_regions = NULL, starting_clusters = NULL, subset = FALSE, N = 5, cl = 1){
  sub.graph = cds@graphs[[lineage]]
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
  mean.dist = .path_distance(nodes_UMAP.sub)
  r = mean.dist*N
  cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
  colnames(cells_UMAP) <- toupper(colnames(cells_UMAP))
  cells_UMAP = cells_UMAP[,c("UMAP_1", "UMAP_2")]
  sel.cells = .cell_selector(nodes_UMAP.sub, cells_UMAP, r, cl = cl)
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

.cell_selector <- function(path, cells, r, cl){
sel.cells = c()
sel.cells = pbapply(path, 1, .selector_sub, cells = cells, r = r, cl = cl, simplify = T)
return(unique(unlist(sel.cells)))
}

.selector_sub <- function(node, cells, r){
x1 = node[1]
y1 = node[2]
res = apply(cells, 1, .cell_selector_sub2, coords = c(x1, y1), r = r, simplify = T)
res = names(res[res == TRUE])
return(res)
}

.cell_selector_sub2 <- function(cell, coords, r){
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

.path_distance <- function(path){
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
#' Subset a metatracker object to one lineage and (re)compute pseudotime
#'
#' Extracts the cells of a lineage, restricts the principal graph and its
#' auxiliary coordinates to that lineage, and — unless \code{recalculate_pt =
#' FALSE} — reprojects cells onto the sub-graph and re-orders pseudotime from
#' the shared start node.
#'
#' @param cds A \code{metatracker_data_set}.
#' @param lineage Lineage name (e.g. "VIP"); \code{FALSE} keeps all cells.
#' @param N Optional cap: randomly downsample to \code{N} cells.
#' @param recalculate_pt Reproject and re-order pseudotime on the subset (default TRUE).
#' @return The subset \code{cds}.
#' @export
get_lineage_object <- function(cds, lineage = FALSE, N = FALSE, recalculate_pt = TRUE) {
  start = .find_start_node(cds)
  if (lineage != FALSE) {
    sub.graph = cds@graphs[[lineage]]
    sel.cells = .lineage_cells(cds@lineages[[lineage]])   # handles vector OR list($name)
  } else {
    sel.cells = colnames(cds)
  }
  sel.cells = sel.cells[sel.cells %in% colnames(cds)]
  if (length(sel.cells) == 0)
    stop("Lineage '", lineage, "' has no cells present in cds.", call. = FALSE)
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  if (N != FALSE) {
    if (N < length(sel.cells)) {
      sel.cells = sample(sel.cells, N)
    }
  }
  # subset the monocle object
  cds_subset = cds[, sel.cells]
  # set the graph, node and cell UMAP coordinates
  if (lineage == FALSE) {
    sub.graph = principal_graph(cds_subset)[["UMAP"]]
  }
  cds_subset@principal_graph[["UMAP"]] <- sub.graph
  cds_subset@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[, names(V(sub.graph))]
  cds_subset@clusters[["UMAP"]]$partitions <- cds_subset@clusters[["UMAP"]]$partitions[colnames(cds_subset)]
  # recalculate closest vertex for the selected cells (vectorised)
  cells_UMAP = as.data.frame(reducedDims(cds_subset)[["UMAP"]])
  colnames(cells_UMAP) <- toupper(colnames(cells_UMAP))
  closest_vertex = .assign_closest_vertex(
    as.matrix(cells_UMAP[, c("UMAP_1", "UMAP_2")]),
    as.matrix(nodes_UMAP[, names(V(sub.graph))]))
  closest_vertex = as.data.frame(closest_vertex)
  cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex
  if (isTRUE(recalculate_pt)) {
    # monocle3-internal projection helpers, pulled from the installed package
    # (the original source_url() download would fail on restricted networks).
    project2MST <- utils::getFromNamespace("project2MST", "monocle3")
    ppls        <- utils::getFromNamespace("project_point_to_line_segment", "monocle3")
    cds_subset  <- project2MST(cds_subset, ppls, FALSE, TRUE, "UMAP",
                               nodes_UMAP[, names(V(sub.graph))])
    cds_subset  <- order_cells(cds_subset, root_pr_nodes = start)
  }
  return(cds_subset)
}

# Most frequent degree-1 (leaf) node across lineage graphs = shared start.
.find_start_node <- function(cds) {
  nodes = c()
  for (name in names(cds@graphs)) {
    sub.graph = cds@graphs[[name]]
    start_end = V(sub.graph)[degree(sub.graph) == 1]$name
    nodes = append(nodes, start_end)
  }
  nodes = as.character(nodes)
  start = names(sort(table(nodes), decreasing = TRUE)[1])
  start
}

# Nearest graph vertex for every cell, vectorised.
# cells: M x 2 (cell UMAP coords); nodes: 2 x K (dims x graph nodes, Y_* names).
# Returns an integer vector (node numbers) named by cell, matching the old
# per-cell apply() output but computed with BLAS matrix ops in chunks.
.assign_closest_vertex <- function(cells, nodes, chunk = 20000) {
  nodesK     <- t(nodes)                 # K x 2
  node_sq    <- rowSums(nodesK^2)        # K
  node_names <- colnames(nodes)          # "Y_1", ...
  M   <- nrow(cells)
  out <- integer(M)
  for (s in seq(1, M, by = chunk)) {
    idx   <- s:min(s + chunk - 1L, M)
    # ||c - n||^2 argmin over n <=> argmax of (2 c.n - ||n||^2)
    score <- 2 * (cells[idx, , drop = FALSE] %*% t(nodesK))
    score <- sweep(score, 2, node_sq, "-")
    nn    <- max.col(score, ties.method = "first")
    out[idx] <- as.integer(gsub("Y_", "", node_names[nn]))
  }
  names(out) <- rownames(cells)
  out
}
