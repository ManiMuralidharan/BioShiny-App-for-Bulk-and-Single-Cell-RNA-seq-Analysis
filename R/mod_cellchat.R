# mod_cellchat.R ---------------------------------------------------------
# Cell-cell communication with CellChat on the processed Seurat object.
# Reads: store$sc. Writes: store$cellchat.
# CellChat is installed from GitHub (see install.R).

mod_cellchat_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 3, title = "Cell-cell communication (CellChat)", status = "primary", solidHeader = TRUE,
      uiOutput(ns("setup")),
      selectInput(ns("species"), "Ligand-receptor database", c("Human" = "human", "Mouse" = "mouse")),
      selectInput(ns("db_subset"), "Interaction types",
                  c("All" = "all", "Secreted signalling" = "Secreted Signaling",
                    "Cell-cell contact" = "Cell-Cell Contact", "ECM-receptor" = "ECM-Receptor")),
      numericInput(ns("min_cells"), "Min cells per group", 10, min = 1),
      actionButton(ns("run"), "Run CellChat", class = "btn-primary", icon = icon("comments")),
      hr(),
      uiOutput(ns("pair_ui")),
      downloadButton(ns("dl"), "Interactions CSV")
    ),
    box(width = 9, title = "Results", status = "info", solidHeader = TRUE,
      htmlOutput(ns("status")),
      tabsetPanel(
        tabPanel("Network (count)", plotOutput(ns("circle_n"), height = "560px")),
        tabPanel("Network (strength)", plotOutput(ns("circle_w"), height = "560px")),
        tabPanel("Heatmap", plotOutput(ns("heat"), height = "560px")),
        tabPanel("L-R bubble plot", plotOutput(ns("bubble"), height = "620px")),
        tabPanel("Interaction table", DTOutput(ns("table")))
      )
    )
  )
}

mod_cellchat_server <- function(id, store) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$setup <- renderUI({
      need_step(store$sc, "Process single-cell data in the Single-cell tab first.")
      md <- store$sc[[]]
      cols <- names(md)[sapply(md, function(x) is.factor(x) || is.character(x))]
      selectInput(ns("group_col"), "Group cells by", cols, selected = "seurat_clusters")
    })

    observeEvent(input$run, {
      if (is.null(store$sc)) { showNotification("Process single-cell data first.", type = "error"); return() }
      req(input$group_col)
      res <- run_step("Running CellChat (can take several minutes)", {
        need_pkg("CellChat")
        seu <- store$sc
        # CellChat dislikes labels such as "0"; prefix numeric cluster IDs.
        lab <- as.character(seu[[]][[input$group_col]])
        if (all(grepl("^[0-9]+$", lab))) lab <- paste0("C", lab)
        lv <- unique(lab)[order(suppressWarnings(as.numeric(sub("^C", "", unique(lab)))), unique(lab))]
        meta <- data.frame(labels = factor(lab, levels = lv), row.names = colnames(seu))
        data_norm <- SeuratObject::LayerData(seu, assay = "RNA", layer = "data")

        cc <- CellChat::createCellChat(object = data_norm, meta = meta, group.by = "labels")
        db <- if (input$species == "human") CellChat::CellChatDB.human else CellChat::CellChatDB.mouse
        if (input$db_subset != "all") db <- CellChat::subsetDB(db, search = input$db_subset, key = "annotation")
        cc@DB <- db
        cc <- CellChat::subsetData(cc)
        cc <- CellChat::identifyOverExpressedGenes(cc)
        cc <- CellChat::identifyOverExpressedInteractions(cc)
        cc <- CellChat::computeCommunProb(cc, type = "triMean")
        cc <- CellChat::filterCommunication(cc, min.cells = input$min_cells)
        cc <- CellChat::computeCommunProbPathway(cc)
        cc <- CellChat::aggregateNet(cc)
        tbl <- CellChat::subsetCommunication(cc)
        if (!nrow(tbl)) stop("No significant interactions found. Try 'All' interaction types or a lower min cells.")
        list(obj = cc, table = tbl, groups = levels(cc@idents), group_col = input$group_col)
      })
      req(res)
      store$cellchat <- res
    })

    cc <- reactive({ need_step(store$cellchat, "Click 'Run CellChat'."); store$cellchat })

    output$status <- renderUI({
      req(store$cellchat); t <- store$cellchat$table
      HTML(sprintf("<p><b>%d</b> significant ligand-receptor interactions across <b>%d</b> pathways between <b>%d</b> groups.</p>",
                   nrow(t), length(unique(t$pathway_name)), length(store$cellchat$groups)))
    })
    output$pair_ui <- renderUI({
      g <- cc()$groups
      tagList(
        selectInput(ns("sources"), "Bubble plot: sender groups", g, selected = g[1], multiple = TRUE),
        selectInput(ns("targets"), "Bubble plot: receiver groups", g, selected = g, multiple = TRUE)
      )
    })

    gsize <- reactive(as.numeric(table(cc()$obj@idents)))
    output$circle_n <- renderPlot({
      o <- cc()$obj
      CellChat::netVisual_circle(o@net$count, vertex.weight = gsize(), weight.scale = TRUE,
                                 label.edge = FALSE, title.name = "Number of interactions")
    })
    output$circle_w <- renderPlot({
      o <- cc()$obj
      CellChat::netVisual_circle(o@net$weight, vertex.weight = gsize(), weight.scale = TRUE,
                                 label.edge = FALSE, title.name = "Interaction strength")
    })
    output$heat <- renderPlot({
      ComplexHeatmap::draw(CellChat::netVisual_heatmap(cc()$obj, measure = "weight"))
    })
    output$bubble <- renderPlot({
      req(input$sources, input$targets)
      o <- cc()$obj; g <- cc()$groups
      CellChat::netVisual_bubble(o, sources.use = match(input$sources, g),
                                 targets.use = match(input$targets, g), remove.isolate = TRUE)
    })
    output$table <- renderDT({
      t <- cc()$table
      keep <- intersect(c("source", "target", "ligand", "receptor", "prob", "pval", "interaction_name", "pathway_name", "annotation"), names(t))
      datatable(t[, keep], rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 15)) %>%
        formatSignif("prob", 3)
    })

    output$dl <- downloadHandler(function() "cellchat_interactions.csv",
      function(file) { req(store$cellchat); write.csv(store$cellchat$table, file, row.names = FALSE) })
  })
}
