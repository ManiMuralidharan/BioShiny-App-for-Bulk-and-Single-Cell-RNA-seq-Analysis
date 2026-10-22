# mod_singlecell.R -------------------------------------------------------
# Seurat workflow: load -> QC -> filter/normalise/PCA/cluster/UMAP -> markers.
# Writes: store$sc (processed Seurat object), store$sc_markers.
# Downstream modules (trajectory, cellchat) read store$sc and use the
# cluster column `seurat_clusters` plus any other metadata columns.

mod_singlecell_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 3, title = "Single-cell (Seurat)", status = "primary", solidHeader = TRUE,
      h5(strong("1. Load")),
      radioButtons(ns("source"), NULL,
                   c("Upload" = "upload", "Demo: pbmc_small (80 cells)" = "demo")),
      conditionalPanel(sprintf("input['%s'] == 'upload'", ns("source")),
        fileInput(ns("file"), "10X .zip (matrix/barcodes/features), 10X .h5, or Seurat .rds",
                  accept = c(".zip", ".h5", ".rds"))),
      actionButton(ns("load"), "Load", class = "btn-primary", icon = icon("upload")),
      hr(),
      h5(strong("2. Filter & cluster")),
      numericInput(ns("min_feat"), "Min genes per cell", 200, min = 0),
      numericInput(ns("max_feat"), "Max genes per cell", 6000, min = 100),
      numericInput(ns("max_mt"), "Max % mitochondrial", 20, min = 0, max = 100),
      numericInput(ns("npcs"), "PCs to use", 20, min = 2, max = 100),
      sliderInput(ns("res"), "Clustering resolution", 0.1, 2, 0.5, step = 0.1),
      actionButton(ns("process"), "Process", class = "btn-danger", icon = icon("cogs")),
      hr(),
      h5(strong("3. Marker genes")),
      actionButton(ns("markers"), "Find markers", class = "btn-success", icon = icon("search")),
      br(), br(),
      downloadButton(ns("dl_markers"), "Markers CSV"),
      downloadButton(ns("dl_rds"), "Seurat .rds")
    ),
    box(width = 9, title = "Results", status = "info", solidHeader = TRUE,
      htmlOutput(ns("status")),
      tabsetPanel(id = ns("tabs"),
        tabPanel("QC", plotOutput(ns("qc_violin"), height = "380px"), plotOutput(ns("qc_scatter"), height = "320px")),
        tabPanel("Elbow", plotOutput(ns("elbow"), height = "420px")),
        tabPanel("UMAP", uiOutput(ns("umap_ctrl")), plotlyOutput(ns("umap"), height = "560px")),
        tabPanel("Markers", DTOutput(ns("marker_tbl")), plotOutput(ns("marker_dot"), height = "520px")),
        tabPanel("Gene expression", selectizeInput(ns("genes"), "Genes (up to 4)", NULL, multiple = TRUE,
                                                   options = list(maxItems = 4)),
                 plotOutput(ns("feature"), height = "560px"))
      )
    )
  )
}

# ---- loading helpers ----------------------------------------------------
read_sc_file <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  if (ext == "rds") {
    obj <- readRDS(path)
    if (inherits(obj, "SingleCellExperiment")) obj <- Seurat::as.Seurat(obj, counts = "counts", data = NULL)
    if (!inherits(obj, "Seurat")) stop("The .rds file does not contain a Seurat or SingleCellExperiment object.")
    return(obj)
  }
  if (ext == "h5") {
    need_pkg("hdf5r")
    m <- Seurat::Read10X_h5(path)
  } else if (ext == "zip") {
    d <- tempfile("tenx_"); utils::unzip(path, exdir = d)
    hit <- list.files(d, pattern = "^matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)
    if (!length(hit)) stop("No matrix.mtx(.gz) found inside the zip.")
    m <- Seurat::Read10X(dirname(hit[1]))
  } else stop("Unsupported file type: ", ext)
  if (is.list(m)) m <- m[["Gene Expression"]] %||% m[[1]]   # multi-modal 10X output
  Seurat::CreateSeuratObject(m, min.cells = 3, min.features = 1)
}

add_qc <- function(seu) {
  genes <- rownames(seu)
  pat <- if (any(grepl("^MT-", genes))) "^MT-" else if (any(grepl("^mt-", genes))) "^mt-" else NULL
  seu[["percent.mt"]] <- if (is.null(pat)) 0 else Seurat::PercentageFeatureSet(seu, pattern = pat)
  seu
}

process_seurat <- function(seu, min_feat, max_feat, max_mt, npcs, res) {
  keep <- seu$nFeature_RNA >= min_feat & seu$nFeature_RNA <= max_feat & seu$percent.mt <= max_mt
  if (sum(keep) < 30) stop(sprintf("Only %d cells pass the QC filters; relax the thresholds.", sum(keep)))
  seu <- subset(seu, cells = colnames(seu)[keep])
  npcs <- min(npcs, ncol(seu) - 1, 50)
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = min(2000, nrow(seu)), verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = npcs, verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, dims = 1:npcs, k.param = min(20, ncol(seu) - 1), verbose = FALSE)
  seu <- Seurat::FindClusters(seu, resolution = res, verbose = FALSE)
  seu <- Seurat::RunUMAP(seu, dims = 1:npcs, n.neighbors = min(30, ncol(seu) - 1), verbose = FALSE)
  seu
}

mod_singlecell_server <- function(id, store) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    raw <- reactiveVal(NULL)   # loaded, unfiltered object (local to this module)

    observeEvent(input$load, {
      obj <- run_step("Loading single-cell data", {
        need_pkg("Seurat")
        if (input$source == "demo") {
          need_pkg("SeuratObject")
          e <- new.env(); utils::data("pbmc_small", package = "SeuratObject", envir = e)
          o <- e$pbmc_small
          o <- Seurat::CreateSeuratObject(SeuratObject::LayerData(o, assay = "RNA", layer = "counts"),
                                          meta.data = o[[c("groups")]])
          o
        } else {
          req(input$file); read_sc_file(input$file$datapath, input$file$name)
        }
      })
      req(obj)
      obj <- add_qc(obj)
      raw(obj)
      store$sc <- NULL; store$sc_markers <- NULL; store$traj <- NULL; store$cellchat <- NULL
      if (input$source == "demo") {   # tiny demo needs looser defaults
        updateNumericInput(session, "min_feat", value = 10)
        updateNumericInput(session, "npcs", value = 10)
        updateSliderInput(session, "res", value = 1)
      }
      showNotification(sprintf("Loaded %d genes x %d cells", nrow(obj), ncol(obj)), type = "message")
    })

    observeEvent(input$process, {
      if (is.null(raw())) { showNotification("Load data first.", type = "error"); return() }
      seu <- run_step("Filtering, normalising, clustering, UMAP",
                      process_seurat(raw(), input$min_feat, input$max_feat, input$max_mt, input$npcs, input$res))
      req(seu)
      store$sc <- seu
      store$sc_markers <- NULL; store$traj <- NULL; store$cellchat <- NULL
      updateTabsetPanel(session, "tabs", "UMAP")
      updateSelectizeInput(session, "genes", choices = rownames(seu), server = TRUE,
                           selected = head(Seurat::VariableFeatures(seu), 2))
    })

    observeEvent(input$markers, {
      if (is.null(store$sc)) { showNotification("Process the data first.", type = "error"); return() }
      mk <- run_step("Finding cluster markers", {
        Seurat::Idents(store$sc) <- "seurat_clusters"
        Seurat::FindAllMarkers(store$sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                               max.cells.per.ident = 500, verbose = FALSE)
      })
      req(mk)
      store$sc_markers <- mk
      updateTabsetPanel(session, "tabs", "Markers")
    })

    output$status <- renderUI({
      txt <- if (!is.null(store$sc)) sprintf("<p><b>%d cells</b> after QC, <b>%d clusters</b>.</p>",
                                             ncol(store$sc), nlevels(store$sc$seurat_clusters))
             else if (!is.null(raw())) sprintf("<p>Loaded <b>%d cells</b>. Check QC, then click Process.</p>", ncol(raw()))
             else "<p>Load data to begin.</p>"
      HTML(txt)
    })

    output$qc_violin <- renderPlot({
      validate(need(!is.null(raw()), "Load data to see QC."))
      Seurat::VlnPlot(raw(), c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1) &
        ggplot2::theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    })
    output$qc_scatter <- renderPlot({
      req(raw())
      df <- raw()[[]]
      df$pass <- df$nFeature_RNA >= input$min_feat & df$nFeature_RNA <= input$max_feat & df$percent.mt <= input$max_mt
      ggplot(df, aes(nCount_RNA, nFeature_RNA, colour = pass)) + geom_point(size = 0.8, alpha = 0.6) +
        scale_colour_manual(values = c(`TRUE` = "grey30", `FALSE` = "red")) + scale_x_log10() + theme_bw() +
        labs(colour = "Passes filters", title = sprintf("%d of %d cells pass current filters", sum(df$pass), nrow(df)))
    })
    output$elbow <- renderPlot({
      need_step(store$sc, "Process the data first.")
      Seurat::ElbowPlot(store$sc, ndims = ncol(Seurat::Embeddings(store$sc, "pca")))
    })

    output$umap_ctrl <- renderUI({
      req(store$sc)
      md <- store$sc[[]]
      cols <- names(md)[sapply(md, function(x) is.factor(x) || is.character(x))]
      selectInput(ns("color_by"), "Colour by", cols, selected = "seurat_clusters")
    })
    output$umap <- renderPlotly({
      need_step(store$sc, "Process the data first.")
      req(input$color_by)
      df <- as.data.frame(Seurat::Embeddings(store$sc, "umap")); names(df) <- c("UMAP1", "UMAP2")
      df$group <- factor(store$sc[[]][[input$color_by]]); df$cell <- colnames(store$sc)
      ggplotly(ggplot(df, aes(UMAP1, UMAP2, colour = group, text = cell)) + geom_point(size = 0.8) +
                 theme_bw() + labs(colour = input$color_by), tooltip = c("text", "colour"))
    })

    output$marker_tbl <- renderDT({
      need_step(store$sc_markers, "Click 'Find markers'.")
      datatable(store$sc_markers, rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 10)) %>%
        formatSignif(c("p_val", "avg_log2FC", "p_val_adj"), 3)
    })
    output$marker_dot <- renderPlot({
      req(store$sc_markers, store$sc)
      top <- store$sc_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 3, with_ties = FALSE)
      Seurat::DotPlot(store$sc, features = unique(top$gene), group.by = "seurat_clusters") +
        Seurat::RotatedAxis() + labs(title = "Top 3 markers per cluster")
    })
    output$feature <- renderPlot({
      need_step(store$sc, "Process the data first.")
      validate(need(length(input$genes) > 0, "Choose one or more genes."))
      Seurat::FeaturePlot(store$sc, features = input$genes, order = TRUE)
    })

    output$dl_markers <- downloadHandler(function() "sc_markers.csv",
      function(file) { req(store$sc_markers); write.csv(store$sc_markers, file, row.names = FALSE) })
    output$dl_rds <- downloadHandler(function() "seurat_processed.rds",
      function(file) { req(store$sc); saveRDS(store$sc, file) })
  })
}
