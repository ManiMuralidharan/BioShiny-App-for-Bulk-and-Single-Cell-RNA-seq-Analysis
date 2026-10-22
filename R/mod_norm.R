# mod_norm.R -------------------------------------------------------------
# Reads: store$counts (or store$batch_corrected if present, for display).
# Writes: store$norm (log-scale matrix), store$norm_method.

mod_norm_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 3, title = "Normalisation", status = "primary", solidHeader = TRUE,
      selectInput(ns("method"), "Method",
                  c("logCPM, TMM (edgeR)" = "TMM", "log2 normalised counts, RLE (DESeq2)" = "RLE",
                    "VST (DESeq2)" = "VST")),
      helpText("This matrix is used for plots, heatmaps and PCA. DE always uses raw counts with its own internal normalisation."),
      actionButton(ns("run"), "Normalise", class = "btn-primary", icon = icon("balance-scale")),
      hr(),
      downloadButton(ns("download"), "Download matrix")
    ),
    box(width = 9, title = "Diagnostics", status = "info", solidHeader = TRUE,
      tabsetPanel(
        tabPanel("Boxplot", plotOutput(ns("box"), height = "450px")),
        tabPanel("Density", plotOutput(ns("dens"), height = "450px")),
        tabPanel("Mean-SD", plotOutput(ns("meansd"), height = "450px")),
        tabPanel("PCA", plotlyOutput(ns("pca"), height = "500px"))
      )
    )
  )
}

normalise_counts <- function(counts, method) {
  switch(method,
    TMM = {
      need_pkg("edgeR")
      d <- edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM")
      edgeR::cpm(d, log = TRUE, prior.count = 2)
    },
    RLE = {
      need_pkg("DESeq2")
      dds <- DESeq2::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(counts, data.frame(row.names = colnames(counts), x = rep(1, ncol(counts))), ~1))
      log2(DESeq2::counts(dds, normalized = TRUE) + 1)
    },
    VST = {
      need_pkg("DESeq2"); need_pkg("SummarizedExperiment")
      dds <- DESeq2::DESeqDataSetFromMatrix(counts, data.frame(row.names = colnames(counts), x = rep(1, ncol(counts))), ~1)
      fun <- if (nrow(counts) < 1000) DESeq2::varianceStabilizingTransformation else DESeq2::vst
      SummarizedExperiment::assay(fun(dds, blind = TRUE))
    })
}

mod_norm_server <- function(id, store) {
  moduleServer(id, function(input, output, session) {

    observeEvent(input$run, {
      need_step(store$counts, "Load data first.")
      src <- store$batch_corrected %||% store$counts
      m <- run_step(paste("Normalising:", input$method), normalise_counts(src, input$method))
      req(m)
      store$norm <- m
      store$norm_method <- paste0(input$method, if (!is.null(store$batch_corrected)) " (on ComBat-Seq counts)" else "")
      showNotification(paste("Normalised with", store$norm_method), type = "message")
    })

    ready <- reactive({ need_step(store$norm, "Click 'Normalise' to compute the normalised matrix."); TRUE })
    grp <- reactive({ if ("condition" %in% colnames(store$meta)) store$meta$condition else factor(store$meta[[2]]) })

    output$box <- renderPlot({ ready()
      boxplot(store$norm, las = 2, outline = FALSE, col = as.integer(grp()) + 1, main = store$norm_method) })
    output$dens <- renderPlot({ ready()
      m <- store$norm; d <- apply(m, 2, density)
      plot(NA, xlim = range(sapply(d, `[[`, "x")), ylim = c(0, max(sapply(d, function(z) max(z$y)))),
           xlab = "Expression", ylab = "Density", main = store$norm_method)
      for (i in seq_along(d)) lines(d[[i]], col = as.integer(grp())[i] + 1)
      legend("topright", levels(grp()), col = seq_along(levels(grp())) + 1, lty = 1, bty = "n") })
    output$meansd <- renderPlot({ ready()
      m <- store$norm; mu <- rowMeans(m); s <- apply(m, 1, sd)
      smoothScatter(rank(mu), s, xlab = "Rank of mean", ylab = "SD", main = "Mean-SD (flat = variance stabilised)")
      lines(lowess(rank(mu), s), col = "red", lwd = 2) })
    output$pca <- renderPlotly({ ready()
      df <- pca_df(store$norm, store$meta); df$group <- grp(); ve <- attr(df, "var")
      ggplotly(ggplot(df, aes(PC1, PC2, color = group, text = sample)) + geom_point(size = 4) + theme_bw() +
                 labs(x = sprintf("PC1 (%s%%)", ve[1]), y = sprintf("PC2 (%s%%)", ve[2])), tooltip = c("text", "colour")) })

    output$download <- downloadHandler(
      filename = function() paste0("normalised_", input$method, ".csv"),
      content = function(file) { req(store$norm); write.csv(store$norm, file) })
  })
}
