# app.R ------------------------------------------------------------------
# BioShiny: modular bulk + single-cell RNA-seq analysis app.
# Run with:  shiny::runApp("path/to/BioShiny")
#
# How the pieces fit:
#   global.R      packages + options
#   R/utils.R     shared helpers
#   R/mod_*.R     one file per analysis (UI function + server function)
#   store         a single reactiveValues object created below and passed to
#                 every module. Modules read their inputs from it and write
#                 their results back, which is how the tabs talk to each other.

source("global.R")

ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = "BioShiny"),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("1. Data", tabName = "data", icon = icon("database")),
      menuItem("2. QC", tabName = "qc", icon = icon("vial")),
      menuItem("3. Normalisation", tabName = "norm", icon = icon("wave-square")),
      menuItem("4. Batch", tabName = "batch", icon = icon("broom")),
      menuItem("5. Differential expression", tabName = "de", icon = icon("chart-line")),
      menuItem("6. Enrichment", tabName = "enrich", icon = icon("project-diagram")),
      menuItem("7. Report", tabName = "report", icon = icon("file-alt")),
      menuItem("Single-cell", icon = icon("braille"), startExpanded = FALSE,
        menuSubItem("8. Seurat: QC & clusters", tabName = "sc"),
        menuSubItem("9. Trajectory", tabName = "traj"),
        menuSubItem("10. Cell-cell comm.", tabName = "cellchat"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML(".content-wrapper{background:#f7f8fa} .box{border-radius:8px}"))),
    tabItems(
      tabItem("data",   mod_data_ui("data")),
      tabItem("qc",     mod_qc_ui("qc")),
      tabItem("norm",   mod_norm_ui("norm")),
      tabItem("batch",  mod_batch_ui("batch")),
      tabItem("de",     mod_de_ui("de")),
      tabItem("enrich", mod_enrich_ui("enrich")),
      tabItem("report", mod_report_ui("report")),
      tabItem("sc",       mod_singlecell_ui("sc")),
      tabItem("traj",     mod_trajectory_ui("traj")),
      tabItem("cellchat", mod_cellchat_ui("cellchat"))
    )
  )
)

server <- function(input, output, session) {
  # Shared state. Every field starts NULL; a module fills its own fields.
  store <- reactiveValues(
    counts = NULL, meta = NULL, organism = NULL, n_genes_raw = NULL,  # mod_data
    norm = NULL, norm_method = NULL,                                  # mod_norm
    batch_corrected = NULL, batch_info = NULL,                        # mod_batch
    de = NULL, de_info = NULL, fdr = 0.05, lfc = 1,                   # mod_de
    enrich = NULL,                                                    # mod_enrich
    sc = NULL, sc_markers = NULL,                                     # mod_singlecell
    traj = NULL,                                                      # mod_trajectory
    cellchat = NULL                                                   # mod_cellchat
  )

  mod_data_server("data", store)
  mod_qc_server("qc", store)
  mod_norm_server("norm", store)
  mod_batch_server("batch", store)
  mod_de_server("de", store)
  mod_enrich_server("enrich", store)
  mod_report_server("report", store)
  mod_singlecell_server("sc", store)
  mod_trajectory_server("traj", store)
  mod_cellchat_server("cellchat", store)
}

shinyApp(ui, server)
