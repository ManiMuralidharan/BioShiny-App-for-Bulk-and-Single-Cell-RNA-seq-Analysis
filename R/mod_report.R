# mod_report.R -----------------------------------------------------------
# Reads everything in the store and renders a real HTML report from
# inst/report.Rmd. The current results are saved to an .rds file and passed
# to the report as a parameter, so the report shows exactly what the app shows.

mod_report_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 4, title = "Analysis report", status = "primary", solidHeader = TRUE,
      textInput(ns("title"), "Title", "RNA-seq analysis report"),
      textInput(ns("author"), "Author", ""),
      textAreaInput(ns("notes"), "Notes (included at the top)", "", rows = 4),
      downloadButton(ns("html"), "Download HTML report"),
      br(), br(),
      downloadButton(ns("rds"), "Download results (.rds)")
    ),
    box(width = 8, title = "What will be included", status = "info", solidHeader = TRUE,
      verbatimTextOutput(ns("contents"))
    )
  )
}

snapshot_store <- function(store) {
  list(counts = store$counts, meta = store$meta, organism = store$organism,
       norm = store$norm, norm_method = store$norm_method,
       de = store$de, de_info = store$de_info, fdr = store$fdr, lfc = store$lfc,
       enrich = store$enrich, batch_info = store$batch_info,
       has_batch = !is.null(store$batch_corrected))
}

mod_report_server <- function(id, store) {
  moduleServer(id, function(input, output, session) {

    output$contents <- renderPrint({
      tick <- function(x) if (is.null(x)) "[ not run ]" else "[ included ]"
      cat("Data          ", tick(store$counts), "\n")
      cat("QC            ", tick(store$counts), "\n")
      cat("Normalisation ", tick(store$norm), if (!is.null(store$norm_method)) store$norm_method, "\n")
      cat("DE            ", tick(store$de), "\n")
      cat("Enrichment    ", tick(store$enrich), "\n")
      cat("Batch         ", tick(store$batch_info), "\n")
      cat("Session info   [ included ]\n")
    })

    output$html <- downloadHandler(
      filename = function() paste0(gsub("[^A-Za-z0-9]+", "_", input$title), ".html"),
      content = function(file) {
        need_pkg("rmarkdown")
        validate(need(!is.null(store$counts), "Load data first."))
        tmp <- tempfile("report_"); dir.create(tmp)
        rmd <- file.path(tmp, "report.Rmd")
        file.copy(file.path("inst", "report.Rmd"), rmd)
        rds <- file.path(tmp, "results.rds")
        saveRDS(snapshot_store(store), rds)
        id <- showNotification("Rendering report...", duration = NULL, type = "message")
        on.exit(removeNotification(id), add = TRUE)
        rmarkdown::render(rmd, output_file = file, envir = new.env(),
                          params = list(results = rds, title = input$title,
                                        author = input$author, notes = input$notes),
                          quiet = TRUE)
      })

    output$rds <- downloadHandler(
      filename = function() "bioshiny_results.rds",
      content = function(file) saveRDS(snapshot_store(store), file))
  })
}
