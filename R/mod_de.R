# mod_de.R ---------------------------------------------------------------
# Reads: store$counts, store$meta (store$norm if available, else computes logCPM for plots).
# Writes: store$de (standardised table), store$de_info (settings used).
# Output columns, whatever the method: gene, baseMean, log2FoldChange, pvalue, padj.

mod_de_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 3, title = "DE setup", status = "primary", solidHeader = TRUE,
      uiOutput(ns("design_ui")),
      selectInput(ns("method"), "Method", c("DESeq2" = "deseq2", "edgeR quasi-likelihood" = "edger", "limma-voom" = "voom")),
      sliderInput(ns("fdr"), "FDR cutoff", 0.01, 0.25, 0.05, step = 0.01),
      sliderInput(ns("lfc"), "|log2FC| cutoff (for calling DEGs)", 0, 3, 1, step = 0.25),
      actionButton(ns("run"), "Run DE", class = "btn-danger", icon = icon("bolt")),
      hr(),
      downloadButton(ns("download"), "Download full table")
    ),
    box(width = 9, title = "Results", status = "info", solidHeader = TRUE,
      htmlOutput(ns("summary")),
      tabsetPanel(
        tabPanel("Table", DTOutput(ns("table"))),
        tabPanel("Volcano", plotlyOutput(ns("volcano"), height = "560px")),
        tabPanel("MA", plotlyOutput(ns("ma"), height = "560px")),
        tabPanel("Top-50 heatmap", plotOutput(ns("heatmap"), height = "750px")),
        tabPanel("Gene", uiOutput(ns("gene_ui")), plotOutput(ns("gene_plot"), height = "420px"))
      )
    )
  )
}

run_de <- function(counts, meta, factor_col, ref, test, covars, method) {
  meta <- droplevels(meta)
  meta$.grp <- relevel(factor(meta[[factor_col]]), ref = ref)
  keep <- meta$.grp %in% c(ref, test)          # only the two compared groups
  counts <- counts[, keep, drop = FALSE]; meta <- droplevels(meta[keep, , drop = FALSE])
  if (any(table(meta$.grp) < 2)) stop("Each group needs at least 2 samples.")
  covars <- setdiff(covars, factor_col)
  form <- as.formula(paste("~", paste(c(covars, ".grp"), collapse = " + ")))
  coef_name <- paste0(".grp", test)

  if (method == "deseq2") {
    need_pkg("DESeq2")
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, meta, form)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    r <- DESeq2::results(dds, contrast = c(".grp", test, ref))
    out <- data.frame(gene = rownames(r), baseMean = r$baseMean, log2FoldChange = r$log2FoldChange,
                      pvalue = r$pvalue, padj = r$padj)
  } else {
    need_pkg("edgeR"); need_pkg("limma")
    design <- model.matrix(form, meta)
    d <- edgeR::calcNormFactors(edgeR::DGEList(counts))
    if (method == "edger") {
      d <- edgeR::estimateDisp(d, design)
      fit <- edgeR::glmQLFit(d, design)
      tt <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = coef_name), n = Inf, sort.by = "none")$table
      out <- data.frame(gene = rownames(tt), baseMean = 2^tt$logCPM, log2FoldChange = tt$logFC,
                        pvalue = tt$PValue, padj = tt$FDR)
    } else {
      v <- limma::voom(d, design)
      fit <- limma::eBayes(limma::lmFit(v, design))
      tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
      out <- data.frame(gene = rownames(tt), baseMean = 2^tt$AveExpr, log2FoldChange = tt$logFC,
                        pvalue = tt$P.Value, padj = tt$adj.P.Val)
    }
  }
  # NA padj (DESeq2 independent filtering / outliers) kept as NA, treated as not significant.
  out[order(out$padj, na.last = TRUE), ]
}

mod_de_server <- function(id, store) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$design_ui <- renderUI({
      need_step(store$meta, "Load data first.")
      cols <- setdiff(colnames(store$meta), "sample")
      cols <- cols[sapply(store$meta[cols], function(x) is.factor(x) && nlevels(x) >= 2 && nlevels(x) < nrow(store$meta))]
      validate(need(length(cols) > 0, "No metadata column with 2+ groups found."))
      tagList(
        selectInput(ns("factor"), "Compare groups in", cols, selected = if ("condition" %in% cols) "condition" else cols[1]),
        uiOutput(ns("levels_ui")),
        selectInput(ns("covars"), "Adjust for (optional, e.g. batch)", setdiff(cols, input$factor), multiple = TRUE)
      )
    })
    output$levels_ui <- renderUI({
      req(input$factor); lv <- levels(store$meta[[input$factor]])
      tagList(
        selectInput(ns("ref"), "Reference (denominator)", lv, selected = lv[1]),
        selectInput(ns("test"), "Test (numerator)", lv, selected = lv[2])
      )
    })

    observeEvent(input$run, {
      need_step(store$counts, "Load data first.")
      req(input$factor, input$ref, input$test)
      if (input$ref == input$test) { showNotification("Reference and test must differ.", type = "error"); return() }
      res <- run_step(paste("Running", input$method),
                      run_de(store$counts, store$meta, input$factor, input$ref, input$test, input$covars, input$method))
      req(res)
      store$de <- res
      store$de_info <- list(method = input$method, factor = input$factor, ref = input$ref, test = input$test,
                            covars = input$covars)
      store$enrich <- NULL   # enrichment is now stale
    })

    # thresholds stay live so the user can explore without re-running DE;
    # Enrichment and Report read these.
    observe({ store$fdr <- input$fdr; store$lfc <- input$lfc })

    de_called <- reactive({
      need_step(store$de, "Set up the comparison and click 'Run DE'.")
      store$de %>% mutate(status = case_when(
        !is.na(padj) & padj < input$fdr & log2FoldChange >  input$lfc ~ "Up",
        !is.na(padj) & padj < input$fdr & log2FoldChange < -input$lfc ~ "Down",
        TRUE ~ "NS"))
    })

    output$summary <- renderUI({
      d <- de_called(); i <- store$de_info
      HTML(sprintf("<p><b>%s</b>: %s <b>%s</b> vs <b>%s</b>%s &nbsp;|&nbsp; <span style='color:#b2182b'>Up %d</span>, <span style='color:#2166ac'>Down %d</span> (FDR &lt; %s, |log2FC| &gt; %s)</p>",
        toupper(i$method), i$factor, i$test, i$ref,
        if (length(i$covars)) paste0(" adjusted for ", paste(i$covars, collapse = ", ")) else "",
        sum(d$status == "Up"), sum(d$status == "Down"), input$fdr, input$lfc))
    })

    output$table <- renderDT({
      d <- de_called()
      datatable(d, filter = "top", rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15)) %>%
        formatSignif(c("baseMean", "log2FoldChange", "pvalue", "padj"), 3)
    })

    cols3 <- c(Up = "#b2182b", Down = "#2166ac", NS = "grey75")
    output$volcano <- renderPlotly({
      d <- de_called() %>% filter(!is.na(pvalue))
      p <- ggplot(d, aes(log2FoldChange, -log10(pvalue), color = status, text = gene)) +
        geom_point(alpha = 0.6, size = 1.2) + scale_color_manual(values = cols3) +
        geom_vline(xintercept = c(-input$lfc, input$lfc), linetype = "dashed", colour = "grey40") +
        theme_bw() + labs(y = "-log10(p-value)", colour = NULL)
      ggplotly(p, tooltip = c("text", "x", "y"))
    })
    output$ma <- renderPlotly({
      d <- de_called()
      p <- ggplot(d, aes(log2(baseMean + 1), log2FoldChange, color = status, text = gene)) +
        geom_point(alpha = 0.6, size = 1.2) + scale_color_manual(values = cols3) +
        geom_hline(yintercept = 0) + theme_bw() + labs(x = "log2 mean expression", colour = NULL)
      ggplotly(p, tooltip = c("text", "x", "y"))
    })
    # Matrix for heatmap/gene plots: the Normalisation tab's result if it was run,
    # otherwise TMM logCPM computed here, so these plots never depend on tab order.
    expr_mat <- reactive({
      if (!is.null(store$norm)) return(list(m = store$norm, label = store$norm_method))
      list(m = normalise_counts(store$counts, "TMM"), label = "logCPM (TMM, computed automatically)")
    })

    output$heatmap <- renderPlot({
      d <- de_called(); need_pkg("pheatmap")
      top <- d %>% filter(status != "NS") %>% head(50)
      validate(need(nrow(top) >= 2, "Fewer than 2 significant genes at these thresholds."))
      em <- expr_mat()$m
      m <- em[intersect(top$gene, rownames(em)), , drop = FALSE]
      ann <- data.frame(group = store$meta[[store$de_info$factor]], row.names = colnames(m))
      pheatmap::pheatmap(m, scale = "row", annotation_col = ann, fontsize_row = 7,
                         color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
                         main = "Top significant genes (row z-score)")
    })
    output$gene_ui <- renderUI({
      d <- de_called()
      selectizeInput(ns("gene"), "Gene", choices = NULL)
    })
    observe({
      d <- de_called()
      updateSelectizeInput(session, "gene", choices = d$gene, selected = d$gene[1], server = TRUE)
    })
    output$gene_plot <- renderPlot({
      req(input$gene, store$de); e <- expr_mat()
      validate(need(input$gene %in% rownames(e$m), "Gene not found."))
      df <- data.frame(expr = e$m[input$gene, ], group = store$meta[[store$de_info$factor]], sample = colnames(e$m))
      ggplot(df, aes(group, expr, colour = group)) + geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.15, size = 3) + theme_bw() + labs(y = e$label, title = input$gene) +
        theme(legend.position = "none")
    })

    output$download <- downloadHandler(
      filename = function() with(store$de_info, sprintf("DE_%s_%s_vs_%s.csv", method, test, ref)),
      content = function(file) write.csv(de_called(), file, row.names = FALSE))
  })
}
