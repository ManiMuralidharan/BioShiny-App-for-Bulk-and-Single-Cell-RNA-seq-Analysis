# mod_enrich.R -----------------------------------------------------------
# Reads: store$de, store$de_info, store$counts (background), store$organism.
# Writes: store$enrich.

mod_enrich_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 3, title = "Functional enrichment", status = "primary", solidHeader = TRUE,
      selectInput(ns("db"), "Analysis",
                  c("GO Biological Process (ORA)" = "BP", "GO Molecular Function (ORA)" = "MF",
                    "GO Cellular Component (ORA)" = "CC", "KEGG (ORA)" = "KEGG", "GSEA on GO BP" = "GSEA_BP")),
      selectInput(ns("direction"), "Genes (ORA only)", c("Up + Down" = "both", "Up only" = "Up", "Down only" = "Down")),
      numericInput(ns("pcut"), "Adjusted p cutoff", 0.05, min = 0.001, max = 0.5, step = 0.01),
      helpText("Background for ORA = all genes that passed filtering (not the whole genome). KEGG needs internet access."),
      actionButton(ns("run"), "Run", class = "btn-primary", icon = icon("dna")),
      hr(),
      downloadButton(ns("download"), "Download table")
    ),
    box(width = 9, title = "Results", status = "info", solidHeader = TRUE,
      htmlOutput(ns("note")),
      tabsetPanel(
        tabPanel("Dot plot", plotOutput(ns("dot"), height = "620px")),
        tabPanel("Table", DTOutput(ns("table"))),
        tabPanel("Gene-concept network", plotOutput(ns("cnet"), height = "620px"))
      )
    )
  )
}

to_entrez <- function(ids, orgdb) {
  kt <- guess_keytype(ids)
  if (kt == "ENTREZID") return(setNames(ids, ids))
  if (kt == "ENSEMBL") ids2 <- strip_version(ids) else ids2 <- ids
  ez <- suppressMessages(AnnotationDbi::mapIds(orgdb, keys = unique(ids2), column = "ENTREZID",
                                               keytype = kt, multiVals = "first"))
  setNames(unname(ez[ids2]), ids)
}

mod_enrich_server <- function(id, store) {
  moduleServer(id, function(input, output, session) {

    observeEvent(input$run, {
      need_step(store$de, "Run DE first.")
      res <- run_step("Running enrichment", {
        need_pkg("clusterProfiler"); need_pkg("AnnotationDbi")
        orgdb <- orgdb_for(store$organism)
        i <- store$de_info; d <- store$de
        ez_all <- to_entrez(d$gene, orgdb)
        if (mean(!is.na(ez_all)) < 0.3) stop("Fewer than 30% of gene IDs mapped to Entrez. Check the organism setting.")
        d$entrez <- ez_all
        sig <- !is.na(d$padj) & d$padj < store$fdr & abs(d$log2FoldChange) > store$lfc
        if (input$direction == "Up") sig <- sig & d$log2FoldChange > 0
        if (input$direction == "Down") sig <- sig & d$log2FoldChange < 0
        genes <- unique(na.omit(d$entrez[sig]))
        universe <- unique(na.omit(to_entrez(rownames(store$counts), orgdb)))

        if (input$db %in% c("BP", "MF", "CC")) {
          if (length(genes) < 5) stop("Fewer than 5 significant genes mapped; relax the DE thresholds.")
          r <- clusterProfiler::enrichGO(genes, OrgDb = orgdb, keyType = "ENTREZID", ont = input$db,
                                         universe = universe, pvalueCutoff = input$pcut, readable = TRUE)
        } else if (input$db == "KEGG") {
          if (length(genes) < 5) stop("Fewer than 5 significant genes mapped; relax the DE thresholds.")
          org <- if (store$organism == "Mm") "mmu" else "hsa"
          r <- clusterProfiler::enrichKEGG(genes, organism = org, universe = universe, pvalueCutoff = input$pcut)
          if (!is.null(r)) r <- clusterProfiler::setReadable(r, orgdb, keyType = "ENTREZID")
        } else {
          # GSEA: rank all tested genes by signed -log10(p)
          rk <- d %>% filter(!is.na(entrez), !is.na(pvalue)) %>%
            mutate(score = sign(log2FoldChange) * -log10(pmax(pvalue, 1e-300))) %>%
            group_by(entrez) %>% slice_max(abs(score), n = 1, with_ties = FALSE) %>% ungroup()
          ranks <- sort(setNames(rk$score, rk$entrez), decreasing = TRUE)
          r <- clusterProfiler::gseGO(ranks, OrgDb = orgdb, keyType = "ENTREZID", ont = "BP",
                                      pvalueCutoff = input$pcut, verbose = FALSE)
          if (!is.null(r)) r <- clusterProfiler::setReadable(r, orgdb, keyType = "ENTREZID")
        }
        list(res = r, n_genes = length(genes), n_univ = length(universe), db = input$db)
      })
      req(res)
      store$enrich <- res
    })

    has_rows <- reactive({
      need_step(store$enrich, "Choose an analysis and click 'Run'.")
      validate(need(!is.null(store$enrich$res) && nrow(as.data.frame(store$enrich$res)) > 0,
                    "No significant terms at this cutoff."))
      TRUE
    })

    output$note <- renderUI({
      req(store$enrich); e <- store$enrich
      HTML(if (startsWith(e$db, "GSEA")) sprintf("<p>GSEA on %d ranked genes.</p>", e$n_univ)
           else sprintf("<p>%d query genes against a background of %d tested genes.</p>", e$n_genes, e$n_univ))
    })
    output$dot <- renderPlot({ has_rows(); need_pkg("enrichplot")
      enrichplot::dotplot(store$enrich$res, showCategory = 20) })
    output$table <- renderDT({ has_rows()
      datatable(as.data.frame(store$enrich$res), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15)) })
    output$cnet <- renderPlot({ has_rows(); need_pkg("enrichplot")
      enrichplot::cnetplot(store$enrich$res, showCategory = 5) })

    output$download <- downloadHandler(
      filename = function() paste0("enrichment_", store$enrich$db, ".csv"),
      content = function(file) write.csv(as.data.frame(store$enrich$res), file, row.names = FALSE))
  })
}
