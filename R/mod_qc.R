# mod_qc.R ---------------------------------------------------------------
# Reads: store$counts, store$meta. Writes nothing.

mod_qc_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 3, title = "QC options", status = "primary", solidHeader = TRUE,
      uiOutput(ns("color_ui")),
      helpText("Plots use log2(count + 1) of the filtered raw counts.")
    ),
    box(width = 9, title = "Sample QC", status = "info", solidHeader = TRUE,
      tabsetPanel(
        tabPanel("Library size", plotlyOutput(ns("libsize"), height = "450px")),
        tabPanel("Detected genes", plotlyOutput(ns("detected"), height = "450px")),
        tabPanel("Distribution", plotOutput(ns("boxplot"), height = "450px")),
        tabPanel("PCA", plotlyOutput(ns("pca"), height = "500px")),
        tabPanel("Correlation", plotOutput(ns("corr"), height = "550px"))
      )
    )
  )
}

mod_qc_server <- function(id, store) {
  moduleServer(id, function(input, output, session) {
    ready <- reactive({ need_step(store$counts, "Load data in the Data tab first."); TRUE })

    output$color_ui <- renderUI({
      req(store$meta)
      cols <- setdiff(colnames(store$meta), "sample")
      selectInput(session$ns("color_by"), "Colour samples by", cols,
                  selected = if ("condition" %in% cols) "condition" else cols[1])
    })
    grp <- reactive({ req(input$color_by); factor(store$meta[[input$color_by]]) })
    logc <- reactive({ ready(); log2(store$counts + 1) })

    output$libsize <- renderPlotly({
      ready()
      df <- data.frame(sample = colnames(store$counts), M = colSums(store$counts) / 1e6, group = grp())
      ggplotly(ggplot(df, aes(sample, M, fill = group)) + geom_col() + theme_bw() +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y = "Reads (millions)", x = NULL))
    })
    output$detected <- renderPlotly({
      ready()
      df <- data.frame(sample = colnames(store$counts), n = colSums(store$counts > 0), group = grp())
      ggplotly(ggplot(df, aes(sample, n, fill = group)) + geom_col() + theme_bw() +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y = "Genes with count > 0", x = NULL))
    })
    output$boxplot <- renderPlot({
      m <- logc()
      boxplot(m, las = 2, outline = FALSE, col = as.integer(grp()) + 1, ylab = "log2(count + 1)")
    })
    output$pca <- renderPlotly({
      df <- pca_df(logc(), store$meta); df$group <- grp(); ve <- attr(df, "var")
      ggplotly(ggplot(df, aes(PC1, PC2, color = group, text = sample)) + geom_point(size = 4) + theme_bw() +
                 labs(x = sprintf("PC1 (%s%%)", ve[1]), y = sprintf("PC2 (%s%%)", ve[2])), tooltip = c("text", "colour"))
    })
    output$corr <- renderPlot({
      need_pkg("pheatmap")
      ann <- data.frame(group = grp(), row.names = colnames(store$counts))
      pheatmap::pheatmap(cor(logc(), method = "pearson"), annotation_col = ann,
                         main = "Sample-sample Pearson correlation")
    })
  })
}
