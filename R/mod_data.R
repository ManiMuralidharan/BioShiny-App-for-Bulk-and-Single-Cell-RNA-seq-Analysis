# mod_data.R -------------------------------------------------------------
# Upload + validation. Writes: store$counts, store$meta, store$organism.

mod_data_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 4, title = "Data import", status = "primary", solidHeader = TRUE,
      radioButtons(ns("source"), "Source",
                   c("Upload count matrix" = "upload", "Demo: airway (human, 8 samples)" = "demo")),
      conditionalPanel(sprintf("input['%s'] == 'upload'", ns("source")),
        fileInput(ns("counts_file"), "Counts (genes x samples; CSV/TSV, first column = gene ID)",
                  accept = c(".csv", ".tsv", ".txt")),
        fileInput(ns("meta_file"), "Metadata (CSV; one row per sample)", accept = ".csv"),
        uiOutput(ns("sample_col_ui"))
      ),
      selectInput(ns("organism"), "Organism", c("Human" = "Hs", "Mouse" = "Mm", "Other" = "other")),
      numericInput(ns("min_count"), "Filter: keep genes with at least this many counts...", 10, min = 0),
      numericInput(ns("min_samples"), "...in at least this many samples", 2, min = 1),
      actionButton(ns("load"), "Load & validate", class = "btn-primary", icon = icon("play"))
    ),
    box(width = 8, title = "Preview", status = "info", solidHeader = TRUE,
      tabsetPanel(
        tabPanel("Summary", verbatimTextOutput(ns("summary"))),
        tabPanel("Counts", DTOutput(ns("counts_tbl"))),
        tabPanel("Metadata", DTOutput(ns("meta_tbl")))
      )
    )
  )
}

mod_data_server <- function(id, store) {
  moduleServer(id, function(input, output, session) {

    meta_raw <- reactive({
      req(input$meta_file)
      read.csv(input$meta_file$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    })

    output$sample_col_ui <- renderUI({
      m <- meta_raw()
      selectInput(session$ns("sample_col"), "Metadata column holding sample IDs",
                  colnames(m), selected = if ("sample" %in% colnames(m)) "sample" else colnames(m)[1])
    })

    read_counts <- function(path, name) {
      sep <- if (tolower(tools::file_ext(name)) == "csv") "," else "\t"
      df <- read.delim(path, sep = sep, row.names = 1, check.names = FALSE)
      as.matrix(df)
    }

    observeEvent(input$load, {
      res <- run_step("Loading data", {
        if (input$source == "demo") {
          need_pkg("airway"); need_pkg("SummarizedExperiment")
          e <- new.env(); utils::data("airway", package = "airway", envir = e)
          se <- e$airway
          counts <- SummarizedExperiment::assay(se)
          meta <- as.data.frame(SummarizedExperiment::colData(se))
          meta$sample <- rownames(meta)
          meta$condition <- relevel(factor(meta$dex), ref = "untrt")
          meta$batch <- factor(meta$cell)   # cell line, useful for batch demo
          meta <- meta[, c("sample", "condition", "batch", setdiff(colnames(meta), c("sample", "condition", "batch")))]
        } else {
          req(input$counts_file, input$meta_file, input$sample_col)
          counts <- read_counts(input$counts_file$datapath, input$counts_file$name)
          meta <- meta_raw()
          meta$sample <- as.character(meta[[input$sample_col]])
        }

        # --- validation ---
        if (!is.numeric(counts)) stop("Count matrix contains non-numeric values.")
        if (any(counts < 0, na.rm = TRUE)) stop("Count matrix contains negative values.")
        if (anyNA(counts)) stop("Count matrix contains missing values.")
        if (any(abs(counts - round(counts)) > 1e-8)) {
          showNotification("Non-integer counts detected (e.g. Salmon/RSEM estimates); rounded to integers.",
                           type = "warning", duration = 8)
          counts <- round(counts)
        }
        storage.mode(counts) <- "integer"

        common <- intersect(colnames(counts), meta$sample)
        if (length(common) < 2) stop("Fewer than 2 sample names match between counts columns and metadata.")
        dropped <- setdiff(union(colnames(counts), meta$sample), common)
        if (length(dropped)) showNotification(paste("Unmatched samples dropped:", paste(dropped, collapse = ", ")),
                                              type = "warning", duration = 10)
        counts <- counts[, common, drop = FALSE]
        meta <- meta[match(common, meta$sample), , drop = FALSE]
        rownames(meta) <- meta$sample
        for (cn in colnames(meta)) if (is.character(meta[[cn]]) && cn != "sample") meta[[cn]] <- factor(meta[[cn]])

        keep <- rowSums(counts >= input$min_count) >= input$min_samples
        list(counts = counts[keep, , drop = FALSE], meta = meta, n_raw = nrow(counts))
      })
      req(res)

      # New data invalidates everything downstream.
      store$counts <- res$counts
      store$meta <- res$meta
      store$organism <- input$organism
      store$n_genes_raw <- res$n_raw
      store$norm <- NULL; store$norm_method <- NULL
      store$de <- NULL; store$de_info <- NULL
      store$enrich <- NULL; store$batch_corrected <- NULL
      showNotification(sprintf("Loaded %d genes x %d samples", nrow(res$counts), ncol(res$counts)), type = "message")
    })

    output$summary <- renderPrint({
      need_step(store$counts, "Load data to begin.")
      cat(sprintf("Genes before filter: %d\nGenes after filter:  %d\nSamples: %d\nGene ID type (guessed): %s\n\nMetadata columns:\n",
                  store$n_genes_raw, nrow(store$counts), ncol(store$counts), guess_keytype(rownames(store$counts))))
      for (cn in setdiff(colnames(store$meta), "sample")) {
        v <- store$meta[[cn]]
        if (is.factor(v)) cat(sprintf("  %-12s %s\n", cn, paste(sprintf("%s (n=%d)", levels(v), tabulate(v)), collapse = ", ")))
      }
    })
    output$counts_tbl <- renderDT({
      need_step(store$counts, "Load data to begin.")
      datatable(head(store$counts, 200), options = list(scrollX = TRUE, pageLength = 10))
    })
    output$meta_tbl <- renderDT({
      need_step(store$meta, "Load data to begin.")
      datatable(store$meta, options = list(scrollX = TRUE))
    })
  })
}
