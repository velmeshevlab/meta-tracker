#extend monocle3 class to add additional slots
#' @export
setClass("metatracker_data_set", contains = "cell_data_set", slots=c(dynamic_genes="list", lineage_genes="list", graphs = "list", lineages="list", expression="list", expectation="list", pseudotime="list")) -> metatracker_data_set

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

#' @export
import_monocle <-function(cds){
cds <- as(cds,"metatracker_data_set")
return(cds)
}

#' @export
# Interactive graph modification with gene expression coloring
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
  
  cds_exprs_all <- SingleCellExperiment::counts(cds)[ , sampled_cells$cell_id, drop = FALSE]
  cds_exprs_all <- t(t(cds_exprs_all) / size_factors(cds)[sampled_cells$cell_id])
  
  el <- as.data.frame(get.edgelist(g), stringsAsFactors = FALSE)
  colnames(el) <- c("from", "to")
  edges_df <- data.frame(
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
        tags$span("Pan:"),
        actionButton("pan_left",   label = NULL, icon = icon("arrow-left")),
        actionButton("pan_up",     label = NULL, icon = icon("arrow-up")),
        actionButton("pan_down",   label = NULL, icon = icon("arrow-down")),
        actionButton("pan_right",  label = NULL, icon = icon("arrow-right")),
        actionButton("done",       "Done")
      )
    ),
    fluidRow(
      column(
        width = 12,
        plotOutput("plot", click = "clickposition", hover = "hoverpos", width = "100%", height = "600px")
      )
    ),
    fluidRow(column(width = 12, verbatimTextOutput("status")))
  )
  
  server <- function(input, output, session) {
    rv <- reactiveValues(
      nodes   = nodes_df,
      edges   = edges_df,
      stage   = 0,
      new_pt  = NULL,
      hovered = NULL,
      xlim    = range(nodes_df$x),
      ylim    = range(nodes_df$y)
    )
    
    output$status <- renderText({
      if (rv$stage == 0) {
        "Step 1: Click to place new node."
      } else {
        "Step 2: Hover and click an existing node to connect."
      }
    })
    
    observeEvent(input$clickposition, {
      pt <- input$clickposition
      if (rv$stage == 0) {
        rv$new_pt  <- c(pt$x, pt$y)
        rv$stage   <- 1
        rv$hovered <- NULL
      } else {
        targ <- if (!is.null(rv$hovered)) rv$hovered else {
          d2 <- (rv$nodes$x - pt$x)^2 + (rv$nodes$y - pt$y)^2
          rv$nodes$node[which.min(d2)]
        }
        max_node = max(as.numeric(str_split_i(rv$nodes$node, "_", 2)))
        new_name <- paste0("Y_", max_node + 1)
        rv$nodes <- rbind(rv$nodes, data.frame(x = rv$new_pt[1], y = rv$new_pt[2], node = new_name, col = "cyan", stringsAsFactors = FALSE))
        rv$edges <- rbind(rv$edges, data.frame(from = new_name, to = targ,
                                               x = rv$new_pt[1], y = rv$new_pt[2],
                                               xend = rv$nodes$x[rv$nodes$node == targ],
                                               yend = rv$nodes$y[rv$nodes$node == targ],
                                               stringsAsFactors = FALSE))
        rv$stage <- 0; rv$new_pt <- NULL; rv$hovered <- NULL
      }
    })
    
    observeEvent(input$hoverpos, {
      if (rv$stage == 1) {
        hv <- input$hoverpos
        d2 <- (rv$nodes$x - hv$x)^2 + (rv$nodes$y - hv$y)^2
        i  <- which.min(d2)
        if (sqrt(d2[i]) < diff(rv$xlim)*0.01 && rv$nodes$col[i] != "cyan") {
          rv$hovered <- rv$nodes$node[i]
        } else {
          rv$hovered <- NULL
        }
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
            geom_segment(data = rv$edges, aes(x = x, y = y, xend = xend, yend = yend), size = segment_size, color = "black") +
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
            geom_segment(data = rv$edges, aes(x = x, y = y, xend = xend, yend = yend), size = segment_size, color = "black") +
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

#' @export
#generate node plot
node_plot <- function(cds, point_size = 1, reduction_method = "UMAP", segment_size = 1){
# 1) Nodes data frame
g = cds@principal_graph[[reduction_method]]
Y <- cds@principal_graph_aux[[reduction_method]]$dp_mst
nodes = as.data.frame(t(Y))
colnames(nodes) <- c("x", "y")
nodes$node <- rownames(nodes)
# 2) Edges data frame
#   get.edgelist(g) returns a two‐column matrix of from/to vertex names (or indices)
el <- as.data.frame(get.edgelist(g), stringsAsFactors = FALSE)
colnames(el) <- c("from", "to")
# 3) Join coordinates
edges <- el %>%
  # join on ‘from’, grabs x,y
  left_join(nodes, by = c("from" = "node")) %>%
  # join on ‘to’, any overlapping names (x,y) get the “.to” suffix
  left_join(nodes, by = c("to"   = "node"), suffix = c("", ".to")) %>%
  # rename the .to columns into xend/yend
  rename(
    xend = x.to,
    yend = y.to
  ) %>%
  # now we have ‘from’, ‘to’, x, y, xend, yend
  select(from, to, x, y, xend, yend)
# 4) Plot
p <- ggplot() +
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    size = segment_size, alpha = 0.5,
    colour = "cyan"
  ) +
  geom_point(
    data = nodes,
    aes(x = x, y = y),
    size = point_size
  ) +
  monocle_theme_opts()
ggplotly(p)
}

#' @export
# Depends on: shiny, igraph, dplyr, ggplot2, plotly, colorspace, SingleCellExperiment, Matrix
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
  cds_exprs_all <- counts(cds)[ , sampled_cells$cell_id, drop = FALSE]
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
            size = segment_size, alpha = 0.3, color = "grey70"
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
sel.cells = isolate_lineage_sub(cds, lineage, sel_clusters = sel_clusters, start_regions = start_regions, starting_clusters = starting_clusters, subset = subset, N = N, cl = cl)
cds@lineages[[lineage]] <- sel.cells
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

find_start_node <- function(cds){
  nodes = c()
  for(name in names(cds@graphs)){
    sub.graph = cds_new@graphs[[name]]
    start_end = V(sub.graph)[degree(sub.graph) == 1]$name
    nodes = append(nodes, start_end)
  }
  nodes = as.character(nodes)
  start = names(sort(table(nodes),decreasing=TRUE)[1])
  start
}

get_lineage_object <- function(cds, lineage = FALSE, N = FALSE){
start = find_start_node(cds)
{
if(lineage != FALSE){
sub.graph = cds@graphs[[lineage]]
sel.cells = cds@lineages[[lineage]]
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
cds_subset@principal_graph[["UMAP"]] <- sub.graph
cds_subset@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(sub.graph))]
cds_subset@clusters[["UMAP"]]$partitions <- cds_subset@clusters[["UMAP"]]$partitions[colnames(cds_subset)]
#recalculate closest vertex for the selected cells
cells_UMAP = as.data.frame(cds_subset@reducedDims[["UMAP"]])
closest_vertex = apply(cells_UMAP[,c("UMAP_1", "UMAP_2")], 1, calculate_closest_vertex, nodes = as.matrix(nodes_UMAP[,names(V(sub.graph))]))
closest_vertex = as.data.frame(closest_vertex)
cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex
source_url("https://raw.githubusercontent.com/cole-trapnell-lab/monocle3/master/R/learn_graph.R")
cds_subset <- project2MST(cds_subset, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP[,names(V(sub.graph))])
cds_subset <- order_cells(cds_subset, root_pr_nodes = start)
return(lineage_cds)
}
  }

isolate_lineage_sub <- function(cds, lineage, sel_clusters = NULL, start_regions = NULL, starting_clusters = NULL, subset = FALSE, N = 5, cl = 1){
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
  mean.dist = path.distance(nodes_UMAP.sub)
  r = mean.dist*N
  cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
  colnames(cells_UMAP) <- toupper(colnames(cells_UMAP))
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
rownames(node_coords) = c("umap_1", "umap_2")
cds@principal_graph_aux[["UMAP"]]$dp_mst <- cbind(cds@principal_graph_aux[["UMAP"]]$dp_mst, node_coords)
graph.new <- add_vertices(graph.old, 1,attr = list(name = new_name))
graph.new <- add_edges(graph.new, c(node1, new_name))
graph.new <- add_edges(graph.new, c(new_name, node2))
}
cds@principal_graph[["UMAP"]] <- graph.new
return(cds)
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

#' @export
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

#' @export
#add a new node
add_node <- function(cds, N = 10, subset = TRUE, node1) {
  nodes_df <-t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  nodes_df = as.data.frame(nodes_df)
  nodes_df <- nodes_df %>% mutate(name = rownames(.))
  g = cds@principal_graph[["UMAP"]]
  edges_df <- get.data.frame(g, what = "edges") %>%
    left_join(nodes_df, by = c("from" = "name")) %>%
    rename(x_start = umap_1, y_start = umap_2) %>%
    left_join(nodes_df, by = c("to" = "name")) %>%
    rename(x_end = umap_1, y_end = umap_2)
  X <- reducedDims(cds)[["UMAP"]]
  if(subset == T){
    X = X[sample(rownames(X), round(nrow(X)/N)),]
  }
  ui <- fluidPage(
    plotlyOutput("scatter"),
    verbatimTextOutput("coords")
  )
  server <- function(input, output, session) {
    
    output$scatter <- renderPlotly({
      ggplotly(
        ggplot(data=X, aes(x=umap_1, y=umap_2), aes_string(umap_1, umap_2, key = seq_len(nrow(X)))) + geom_point(size=0.5) + geom_segment(data = edges_df,aes(x = x_start, y = y_start, xend = x_end, yend = y_end), color = "cyan") + monocle_theme_opts(),
        source = "scatter"
      )
    })
    observeEvent(event_data("plotly_click", source = "scatter"), {
      click <- event_data("plotly_click", source = "scatter")
      # return whichever columns you like:
      stopApp(list(x = click$x,
                   y = click$y
                   ))
    }, ignoreInit = TRUE)
  }
  coords = runApp(shinyApp(ui, server))
  graph.old = cds@principal_graph[["UMAP"]]
  new_name = paste0("Y_", as.character(length(names(V(graph.old)))+1))
  node_coords = as.data.frame(c(coords$x, coords$y))
  colnames(node_coords) = new_name
  rownames(node_coords) = c("umap_1", "umap_2")
  cds@principal_graph_aux[["UMAP"]]$dp_mst <- cbind(cds@principal_graph_aux[["UMAP"]]$dp_mst, node_coords)
  graph.new <- add_vertices(graph.old, 1,attr = list(name = new_name))
  graph.new <- add_edges(graph.new, c(node1, new_name))
  cds@principal_graph[["UMAP"]] <- graph.new
  return(cds)
}

