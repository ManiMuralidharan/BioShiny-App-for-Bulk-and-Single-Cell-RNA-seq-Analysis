# BioShiny

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%E2%89%A54.3-276DC3?logo=r)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18%2B-1a81c2)](https://bioconductor.org/)
[![Shiny](https://img.shields.io/badge/Shiny-app-blue?logo=rstudio)](https://shiny.posit.co/)

**A modular R Shiny application for interactive bulk and single-cell RNA-seq analysis.**

BioShiny takes you from a count matrix to differential expression, pathway enrichment and a rendered report, and from a 10X matrix to clusters, trajectories and cell-cell communication, all in the browser with no coding required.

<!-- Add a screenshot: save it as docs/screenshot.png and uncomment the next line -->
<!-- ![BioShiny screenshot](docs/screenshot.png) -->

---

## Features

### Bulk RNA-seq
- **Data import and validation:** CSV/TSV counts and metadata, sample matching, low-count filtering, automatic rounding of non-integer estimates
- **Quality control:** library size, detected genes, PCA, sample correlation
- **Normalisation:** TMM logCPM (edgeR), RLE and VST (DESeq2)
- **Batch correction:** ComBat-Seq with before/after PCA (raw counts are never overwritten)
- **Differential expression:** DESeq2, edgeR quasi-likelihood or limma-voom, with user-chosen reference and test groups and optional covariates
- **Visualisation:** interactive volcano and MA plots, top-gene heatmap, per-gene plots
- **Functional enrichment:** GO (BP/MF/CC) and KEGG over-representation, GSEA on GO BP (clusterProfiler)
- **Reporting:** one-click rendered HTML report plus an `.rds` of all results

### Single-cell RNA-seq
- **Seurat workflow:** 10X (`.zip` or `.h5`) or Seurat `.rds` input, QC with live filter preview, normalisation, PCA, clustering, UMAP and marker genes
- **Trajectory inference:** Slingshot lineages and pseudotime, root cluster chosen from your real clusters, gene trends along pseudotime
- **Cell-cell communication:** CellChat networks, heatmaps, ligand-receptor bubble plots and interaction tables

---

## Installation

Requires **R ≥ 4.3**.

```bash
git clone https://github.com/ManiMuralidharan/BioShiny.git
cd BioShiny
```

```r
source("install.R")   # installs CRAN, Bioconductor and GitHub dependencies (first run only)
shiny::runApp()
```

`install.R` only installs what is missing. CellChat is installed from GitHub and is needed only for the cell-cell communication tab; the rest of the app works without it.

## Quick start

| Try | Steps |
|---|---|
| Bulk demo | **1. Data** → *Demo: airway* → **Load & validate**, then work through tabs 2 to 7 |
| Single-cell demo | **Single-cell → 8. Seurat** → *Demo: pbmc_small* → **Load** → **Process** → **Find markers** |

---

## Input formats

**Bulk**
- **Counts:** CSV or TSV, genes in rows and samples in columns; the first column is the gene ID (Ensembl, Entrez or symbol). Raw integer counts from featureCounts, HTSeq or STAR are expected; non-integer estimates (Salmon/RSEM) are rounded with a warning.
- **Metadata:** CSV with one row per sample, one column of sample IDs matching the count-matrix column names, plus grouping columns (e.g. `condition`, `batch`).

**Single-cell**
- A `.zip` of a 10X folder (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`), a 10X `.h5` file, or a Seurat / SingleCellExperiment `.rds`.

---

## Workflow

| Tab | Does | Uses | Produces |
|---|---|---|---|
| 1. Data | Upload, validate, filter low counts | files | `counts`, `meta` |
| 2. QC | Library size, detected genes, PCA, correlation | `counts` | plots |
| 3. Normalisation | TMM logCPM, RLE or VST | `counts` (or ComBat-Seq counts) | `norm` |
| 4. Batch | ComBat-Seq, before/after PCA | `counts` | `batch_corrected` |
| 5. DE | DESeq2, edgeR QL or limma-voom; reference, test group and covariates | `counts`, `meta` | `de` |
| 6. Enrichment | GO/KEGG over-representation, GSEA on GO BP | `de` | `enrich` |
| 7. Report | Rendered HTML report of everything above | all | `.html`, `.rds` |
| 8. Seurat | QC, filter, normalise, PCA, cluster, UMAP, markers | 10X / `.rds` | `sc`, `sc_markers` |
| 9. Trajectory | Slingshot lineages, pseudotime, gene trends | `sc` | `traj` |
| 10. Cell-cell comm. | CellChat networks, heatmap, L-R bubble plot | `sc` | `cellchat` |

**Design notes**
- DE always runs on **raw counts** with each method's own normalisation. To account for batch in DE, add it as a covariate in the DE tab rather than using ComBat-Seq counts.
- Re-processing the Seurat object clears trajectory and CellChat results, so stale results are never shown.
- Tabs whose inputs are not ready show a "run X first" message instead of an error or placeholder plot.

---

## Architecture

```
BioShiny/
├── app.R              # UI layout, shared `store`, wires modules together
├── global.R           # packages and options
├── install.R          # dependencies
├── R/
│   ├── utils.R        # shared helpers
│   └── mod_*.R        # one analysis per file (UI + server)
└── inst/report.Rmd    # report template
```

Every module receives the same `store` (a `reactiveValues` object), reads its inputs from it and writes its results back. This is how results flow between tabs.

**Adding an analysis:** create `R/mod_new.R` with `mod_new_ui(id)` and `mod_new_server(id, store)`, then add one `menuItem`, one `tabItem` and one server call in `app.R`.

---

## Limitations

- Two-group comparisons per DE run (one reference vs one test level).
- Enrichment supports human and mouse; KEGG requires internet access.
- Single-cell analysis handles one sample per run (no multi-sample integration yet).
- The trajectory "genes changing with pseudotime" table is an exploratory Spearman ranking, not a formal test (e.g. tradeSeq).

## Roadmap

- [ ] Multi-sample single-cell integration (Harmony)
- [ ] ChIP-seq / ATAC-seq differential binding
- [ ] Variant annotation
- [ ] Multi-omics integration (MOFA2)

---

## Citation

If you use BioShiny in your research, please cite it using the **"Cite this repository"** button on this page (generated from [`CITATION.cff`](CITATION.cff)), or:

> Mani, M. (2026). *BioShiny: a modular R Shiny application for bulk and single-cell RNA-seq analysis* (Version 1.0.0) [Computer software]. https://github.com/ManiMuralidharan/BioShiny

Please also cite the underlying methods you use (e.g. DESeq2, edgeR, limma, sva, clusterProfiler, Seurat, Slingshot, CellChat).

## Author

**Muralidharan Mani, Ph.D.**
Honorary Fellow, Department of Biochemistry, University of Wisconsin-Madison

[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--1194--3616-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0000-0002-1194-3616)
[![GitHub](https://img.shields.io/badge/GitHub-ManiMuralidharan-181717?logo=github)](https://github.com/ManiMuralidharan)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Muralidharan%20Mani-0A66C2?logo=linkedin)](https://www.linkedin.com/in/muralidharan-mani-b4b32093/)

## License

Released under the [MIT License](LICENSE).
