# utils.R ----------------------------------------------------------------
# Small helpers shared by all modules.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Show a friendly placeholder instead of an error when an upstream step
# has not been run yet.
need_step <- function(x, msg) {
  validate(need(!is.null(x), msg))
}

# Check a package is installed before running an analysis.
need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Run install.R first.", pkg), call. = FALSE)
  }
  invisible(TRUE)
}

# Run an analysis step with a progress bar and a readable error notification.
run_step <- function(label, expr) {
  withProgress(message = label, value = 0.3, {
    tryCatch(
      {
        out <- force(expr)
        setProgress(1)
        out
      },
      error = function(e) {
        showNotification(paste(label, "failed:", conditionMessage(e)),
                         type = "error", duration = NULL)
        NULL
      }
    )
  })
}

# Guess the gene ID type from the row names.
guess_keytype <- function(ids) {
  ids <- head(ids[!is.na(ids)], 200)
  if (mean(grepl("^ENS[A-Z]*G[0-9]+", ids)) > 0.5) return("ENSEMBL")
  if (mean(grepl("^[0-9]+$", ids)) > 0.5) return("ENTREZID")
  "SYMBOL"
}

# Strip Ensembl version suffixes (ENSG000001.12 -> ENSG000001).
strip_version <- function(ids) sub("\\.[0-9]+$", "", ids)

orgdb_for <- function(organism) {
  pkg <- switch(organism, Hs = "org.Hs.eg.db", Mm = "org.Mm.eg.db", NULL)
  if (is.null(pkg)) stop("Enrichment supports Human and Mouse only.", call. = FALSE)
  need_pkg(pkg)
  getExportedValue(pkg, pkg)
}

# PCA on a samples-in-columns matrix, using the most variable genes.
pca_df <- function(mat, meta, ntop = 500) {
  mat <- as.matrix(mat)
  v <- apply(mat, 1, var)
  keep <- head(order(v, decreasing = TRUE), min(ntop, sum(v > 0)))
  p <- prcomp(t(mat[keep, , drop = FALSE]), scale. = FALSE)
  ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  df <- cbind(meta, PC1 = p$x[, 1], PC2 = p$x[, 2])
  attr(df, "var") <- ve[1:2]
  df
}
