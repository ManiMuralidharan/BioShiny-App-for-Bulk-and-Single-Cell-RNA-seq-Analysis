# install.R --------------------------------------------------------------
# Run once:  source("install.R")
# Tested with R 4.3+ and Bioconductor 3.18+.

cran <- c("shiny", "shinydashboard", "DT", "ggplot2", "plotly", "dplyr",
          "pheatmap", "rmarkdown", "knitr", "BiocManager",
          "Seurat", "SeuratObject", "hdf5r", "remotes",
          "NMF", "circlize", "ggalluvial")   # last three: CellChat dependencies
bioc <- c("SummarizedExperiment", "DESeq2", "edgeR", "limma", "sva",
          "AnnotationDbi", "clusterProfiler", "enrichplot",
          "org.Hs.eg.db", "org.Mm.eg.db", "airway",
          "SingleCellExperiment", "slingshot", "ComplexHeatmap")

new_cran <- setdiff(cran, rownames(installed.packages()))
if (length(new_cran)) install.packages(new_cran)

new_bioc <- setdiff(bioc, rownames(installed.packages()))
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

# CellChat is GitHub-only (needed for tab 10 only)
if (!requireNamespace("CellChat", quietly = TRUE)) {
  remotes::install_github("jinworks/CellChat", upgrade = "never")
}
# Optional: presto makes Seurat's FindAllMarkers much faster
# remotes::install_github("immunogenomics/presto")

missing <- setdiff(c(cran, bioc, "CellChat"), rownames(installed.packages()))
if (length(missing)) {
  message("Still missing: ", paste(missing, collapse = ", "))
} else {
  message("All packages installed. Start the app with shiny::runApp()")
}
