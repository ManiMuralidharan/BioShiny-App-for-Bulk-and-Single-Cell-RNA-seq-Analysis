# mod_trajectory.R -------------------------------------------------------
# Slingshot lineages + pseudotime on the processed Seurat object.
# Reads: store$sc (needs PCA, UMAP and a cluster column).
# Writes: store$traj (list: pseudotime matrix, lineages, curves, settings).

mod_trajectory_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 3, title = "Trajectory (Slingshot)", status = "primary", solidHeader = TRUE,
      uiOutput(ns("setup")),
      selectInput(ns("space"), "Fit in",
                  c("PCA (recommended)" = "pca", "UMAP (2-D; curves drawn exactly)" = "umap")),
      numericInput(ns("ndims"), "PCs (if PCA)", 10, min = 2, max = 50),
      actionButton(ns("run"), "Run Slingshot", class = "btn-primary", icon = icon("route")),
      hr(),
      helpText("Pseudotime is only meaningful if the selected cells form a continuous process (e.g. differentiation). Choose the root from biology, not convenience."),
      downloadButton(ns("dl"), "Pseudotime CSV")
    ),
    box(width = 9, title = "Results", status = "info", solidHeader = TRUE,
      htmlOutput(ns("status")),
      tabsetPanel(
        tabPanel("Pseudotime UMAP", uiOutput(ns("lin_ui")), plotOutput(ns("umap"), height = "560px")),
        tabPanel("Lineages", verbatimTextOutput(ns("lineages"))),
        tabPanel("Gene trends", selectizeInput(ns("genes"), "Genes (up to 4)", NULL, multiple = TRUE,
                                               options = list(maxItems = 4)),
                 plotOutput(ns("trend"), height = "480px")),
        tabPanel("Genes changing with pseudotime", helpText("Spearman correlation of each variable gene with pseudotime along the selected lineage. Exploratory ranking only; for formal testing use tradeSeq."),
                 DTOutput(ns("assoc")))
      )
    )
  )
}

mod_trajectory_server <- function(id, store) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Cluster choices come from the real Seurat object, not hardcoded labels.
    output$setup <- renderUI({
      need_step(store$sc, "Process single-cell data in the Single-cell tab first.")
      md <- store$sc[[]]
      cols <- names(md)[sapply(md, function(x) is.factor(x) || is.character(x))]
      tagList(
        selectInput(ns("cluster_col"), "Cluster labels", cols, selected = "seurat_clusters"),
        uiOutput(ns("root_ui"))
      )
    })
    output$root_ui <- renderUI({
      req(store$sc, input$cluster_col)
      lv <- levels(factor(store$sc[[]][[input$cluster_col]]))
      selectInput(ns("root"), "Root (start) cluster", lv, selected = lv[1])
    })

    observeEvent(input$run, {
      if (is.null(store$sc)) { showNotification("Process single-cell data first.", type = "error"); return() }
      req(input$cluster_col, input$root)
      res <- run_step("Running Slingshot", {
        need_pkg("slingshot")
        emb <- Seurat::Embeddings(store$sc, input$space)
        if (input$space == "pca") emb <- emb[, seq_len(min(input$ndims, ncol(emb))), drop = FALSE]
        cl <- as.character(store$sc[[]][[input$cluster_col]])
        sds <- slingshot::slingshot(emb, clusterLabels = cl, start.clus = input$root)
        pt <- slingshot::slingPseudotime(sds)
        list(pseudotime = pt, lineages = slingshot::slingLineages(sds),
             curves = if (input$space == "umap") slingshot::slingCurves(sds, as.df = TRUE) else NULL,
             space = input$space, cluster_col = input$cluster_col, root = input$root)
      })
      req(res)
      store$traj <- res
      updateSelectizeInput(session, "genes", choices = rownames(store$sc), server = TRUE,
                           selected = head(Seurat::VariableFeatures(store$sc), 2))
    })

    tr <- reactive({ need_step(store$traj, "Click 'Run Slingshot'."); store$traj })

    output$status <- renderUI({
      req(store$traj)
      HTML(sprintf("<p><b>%d lineage(s)</b> from root cluster <b>%s</b>, fitted in %s.</p>",
                   length(store$traj$lineages), store$traj$root, toupper(store$traj$space)))
    })
    output$lin_ui <- renderUI({
      t <- tr(); selectInput(ns("lineage"), "Lineage", colnames(t$pseudotime))
    })
    lin <- reactive({ t <- tr(); input$lineage %||% colnames(t$pseudotime)[1] })

    output$umap <- renderPlot({
      t <- tr()
      df <- as.data.frame(Seurat::Embeddings(store$sc, "umap")); names(df) <- c("UMAP1", "UMAP2")
      df$pseudotime <- t$pseudotime[, lin()]
      p <- ggplot(df, aes(UMAP1, UMAP2)) +
        geom_point(data = df[is.na(df$pseudotime), ], colour = "grey85", size = 0.7) +
        geom_point(data = df[!is.na(df$pseudotime), ], aes(colour = pseudotime), size = 0.9) +
        scale_colour_viridis_c() + theme_bw() + labs(title = paste(lin(), "(grey = cells not on this lineage)"))
      if (!is.null(t$curves)) {
        cv <- t$curves; names(cv)[1:2] <- c("UMAP1", "UMAP2")
        p <- p + geom_path(data = cv[order(cv$Lineage, cv$Order), ], aes(group = Lineage), linewidth = 1)
      }
      p
    })
    output$lineages <- renderPrint({
      t <- tr(); for (n in names(t$lineages)) cat(n, ":", paste(t$lineages[[n]], collapse = " -> "), "\n")
    })

    output$trend <- renderPlot({
      t <- tr(); validate(need(length(input$genes) > 0, "Choose one or more genes."))
      expr <- as.matrix(SeuratObject::LayerData(store$sc, assay = "RNA", layer = "data")[input$genes, , drop = FALSE])
      pt <- t$pseudotime[, lin()]
      df <- do.call(rbind, lapply(input$genes, function(g) data.frame(gene = g, pseudotime = pt, expr = expr[g, ])))
      df <- df[!is.na(df$pseudotime), ]
      ggplot(df, aes(pseudotime, expr)) + geom_point(alpha = 0.3, size = 0.7) +
        geom_smooth(method = "loess", se = TRUE, colour = "firebrick", formula = y ~ x) +
        facet_wrap(~gene, scales = "free_y") + theme_bw() + labs(y = "Log-normalised expression", title = lin())
    })

    output$assoc <- renderDT({
      t <- tr()
      pt <- t$pseudotime[, lin()]; on <- !is.na(pt)
      genes <- Seurat::VariableFeatures(store$sc)
      expr <- as.matrix(SeuratObject::LayerData(store$sc, assay = "RNA", layer = "data")[genes, on, drop = FALSE])
      rho <- apply(expr, 1, function(x) suppressWarnings(cor(x, pt[on], method = "spearman")))
      df <- data.frame(gene = genes, spearman_rho = round(rho, 3))
      df <- df[order(-abs(df$spearman_rho)), ]
      datatable(df, rownames = FALSE, options = list(pageLength = 15))
    })

    output$dl <- downloadHandler(function() "slingshot_pseudotime.csv",
      function(file) { req(store$traj); write.csv(store$traj$pseudotime, file) })
  })
}
