# mod_batch.R ------------------------------------------------------------
# Reads: store$counts, store$meta.
# Writes: store$batch_corrected (raw counts in store$counts are NEVER overwritten).
#
# Note: for differential expression, adding batch as a covariate in the DE tab
# is the recommended approach. ComBat-Seq corrected counts are mainly for
# visualisation (PCA, heatmaps) and for tools that cannot model covariates.

mod_batch_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 3, title = "Batch correction (ComBat-Seq)", status = "primary", solidHeader = TRUE,
      uiOutput(ns("vars_ui")),
      actionButton(ns("run"), "Apply ComBat-Seq", class = "btn-primary", icon = icon("broom")),
      actionButton(ns("clear"), "Discard correction", icon = icon("undo")),
      hr(),
      helpText("For DE, prefer adding batch as a covariate in the DE tab. After applying, re-run Normalisation so plots use the corrected counts."),
      downloadButton(ns("download"), "Download corrected counts")
    ),
    box(width = 9, title = "Before / after", status = "info", solidHeader = TRUE,
      htmlOutput(ns("status")),
      plotOutput(ns("pca"), height = "480px")
    )
  )
}

mod_batch_server <- function(id, store) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$vars_ui <- renderUI({
      need_step(store$meta, "Load data first.")
      cols <- setdiff(colnames(store$meta), "sample")
      cols <- cols[sapply(store$meta[cols], is.factor)]
      tagList(
        selectInput(ns("batch"), "Batch variable", cols, selected = if ("batch" %in% cols) "batch" else cols[1]),
        selectInput(ns("group"), "Biological group to preserve", cols, selected = if ("condition" %in% cols) "condition" else cols[1])
      )
    })

    observeEvent(input$run, {
      need_step(store$counts, "Load data first."); req(input$batch, input$group)
      if (input$batch == input$group) { showNotification("Batch and group must be different columns.", type = "error"); return() }
      b <- store$meta[[input$batch]]; g <- store$meta[[input$group]]
      if (nlevels(droplevels(b)) < 2) { showNotification("Batch variable has only one level.", type = "error"); return() }
      if (all(rowSums(table(b, g) > 0) == 1)) { showNotification("Batch and group are fully confounded (each batch has one group); correction would remove the biology.", type = "error"); return() }
      res <- run_step("Running ComBat-Seq", {
        need_pkg("sva")
        out <- suppressMessages(sva::ComBat_seq(as.matrix(store$counts), batch = b, group = g))
        dimnames(out) <- dimnames(store$counts); out
      })
      req(res)
      store$batch_corrected <- res
      store$batch_info <- c(batch = input$batch, group = input$group)
      showNotification("ComBat-Seq applied. Re-run Normalisation to update downstream plots.", type = "message")
    })
    observeEvent(input$clear, { store$batch_corrected <- NULL; store$batch_info <- NULL })

    output$status <- renderUI({
      HTML(if (is.null(store$batch_corrected)) "<p>No correction applied.</p>"
           else sprintf("<p>Corrected for <b>%s</b>, preserving <b>%s</b>.</p>", store$batch_info["batch"], store$batch_info["group"]))
    })

    output$pca <- renderPlot({
      need_step(store$counts, "Load data first."); req(input$batch, input$group)
      mk <- function(m, lab) { df <- pca_df(log2(m + 1), store$meta); df$panel <- lab; df }
      df <- mk(store$counts, "Before")
      if (!is.null(store$batch_corrected)) df <- rbind(df, mk(store$batch_corrected, "After ComBat-Seq"))
      df$panel <- factor(df$panel, c("Before", "After ComBat-Seq"))
      ggplot(df, aes(PC1, PC2, colour = .data[[input$batch]], shape = .data[[input$group]])) +
        geom_point(size = 4) + facet_wrap(~panel, scales = "free") + theme_bw() +
        labs(colour = input$batch, shape = input$group, title = "PCA of log2(count + 1), top 500 variable genes")
    })

    output$download <- downloadHandler(
      filename = function() "combatseq_corrected_counts.csv",
      content = function(file) { req(store$batch_corrected); write.csv(store$batch_corrected, file) })
  })
}
