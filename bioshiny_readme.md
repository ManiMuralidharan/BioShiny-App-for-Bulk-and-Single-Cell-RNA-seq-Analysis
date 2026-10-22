# BioShiny

**BioShiny** is a modular R Shiny application built for comprehensive, end-to-end RNA-sequencing data analysis and visualization. It provides interactive workflows for both **Bulk RNA-seq** and **Single-Cell RNA-seq (scRNA-seq)**, allowing users to move from raw count matrices to differential expression, pathway enrichment, trajectory inference, and cell-cell communication networks—all without writing code.

---

## Key Features

- **Bulk RNA-seq Workflow**:
  - **Data Import & Validation**: Load raw count matrices and sample metadata with automatic filtering of low-count genes.
  - **Quality Control (QC)**: Library size distribution, detected gene counts, Principal Component Analysis (PCA), and sample correlation heatmaps.
  - **Normalization & Batch Correction**: Supports TMM logCPM, RLE, VST, and batch effect adjustment via ComBat-Seq.
  - **Differential Expression (DE)**: Integrated DESeq2, edgeR (QL), and limma-voom engines with flexible covariate and reference group selection.
  - **Functional Enrichment**: Gene Ontology (GO) and KEGG over-representation analysis, along with Gene Set Enrichment Analysis (GSEA) on GO Biological Processes.
  - **Reporting**: One-click generation of fully rendered HTML summaries and downloadable RDS objects.

- **Single-Cell RNA-seq (scRNA-seq) Workflow**:
  - **Seurat Pipeline**: Import 10X Genomics format (`.zip` matrix/barcodes/features or `.h5`) or pre-processed `.rds` files; execute QC filtering, normalization, PCA, UMAP dimensionality reduction, clustering, and marker identification.
  - **Trajectory Inference**: Reconstruct cell lineages and pseudotime trends using Slingshot directly on Seurat objects.
  - **Cell-Cell Communication**: Analyze ligand-receptor interaction networks, interaction heatmaps, and bubble plots powered by CellChat.

---

## Installation & Running

### Prerequisites
Ensure you have **R (>= 4.0)** installed on your machine.

### Quick Start

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/your-username/BioShiny.git
   cd BioShiny
   ```

2. **Install Dependencies**:
   Run the installation script once to fetch required CRAN, Bioconductor, and GitHub packages:
   ```r
   source("install.R")
   ```

3. **Launch the Application**:
   Launch the app directly from inside the `BioShiny` directory:
   ```r
   shiny::runApp()
   ```

---

## Input Formats

### 1. Bulk RNA-seq
* **Count Matrix**: CSV or TSV format where rows represent genes and columns represent samples.
  * The first column must contain Gene IDs (Ensembl, Entrez, or Gene Symbols).
  * Inputs must be raw integer counts (e.g., from `featureCounts`, `HTSeq`, or `STAR`). Non-integer estimates (e.g., Salmon/RSEM) will be automatically rounded with a warning.
* **Metadata**: CSV format containing one row per sample. Must include a column with Sample IDs matching the column names of the count matrix, alongside experimental covariates (e.g., `condition`, `batch`).

### 2. Single-Cell RNA-seq
* **10X Data**: `.zip` archive containing `matrix.mtx`, `barcodes.tsv`, and `features.tsv` (or `genes.tsv`), or a standard 10X `.h5` file.
* **Seurat Object**: Pre-saved `.rds` file containing a valid Seurat object.

> **Tip**: You can explore the app immediately using built-in demonstration datasets:
> - Bulk RNA-seq Demo: **Data > Demo: airway**
> - Single-Cell Demo: **Single-cell > Demo: pbmc_small** (80 cells)

---

## Analysis Workflows

### Bulk RNA-seq Modules

| Tab | Functionality | Input | Output |
|---|---|---|---|
| **1. Data** | File upload, format validation, low-count gene filtering | Raw files | `counts`, `meta` |
| **2. QC** | Library size, gene detection rates, PCA, sample correlation heatmaps | `counts` | Diagnostic plots |
| **3. Normalisation** | Calculates TMM logCPM, RLE, or VST normalized counts | `counts` (or ComBat-Seq) | `norm` |
| **4. Batch** | Batch effect correction via ComBat-Seq with pre/post PCA | `counts` | `batch_corrected` |
| **5. DE** | Differential expression via DESeq2, edgeR QL, or limma-voom | `counts`, `meta` | `de` |
| **6. Enrichment** | GO/KEGG over-representation analysis, GSEA on GO BP | `de` | `enrich` |
| **7. Report** | Generates a comprehensive HTML summary report & exports objects | All analysis steps | `.html`, `.rds` |

*Note: Differential Expression (DE) tools perform their own internal normalization on raw counts. To account for batch effects in DE, include batch as a covariate in the DE tab rather than passing ComBat-Seq adjusted counts.*

### Single-Cell RNA-seq Modules

| Tab | Functionality | Input | Output |
|---|---|---|---|
| **8. Seurat** | QC filtering, normalization, PCA, UMAP, clustering, marker discovery | Raw file / `.rds` | `sc`, `sc_markers` |
| **9. Trajectory** | Slingshot lineage inference, pseudotime calculation, gene trends | `sc` | `traj` |
| **10. Cell-Cell Comm.** | CellChat communication networks, heatmaps, L-R bubble plots | `sc` | `cellchat` |

*Note: Re-processing or re-clustering the Seurat object in Tab 8 will automatically clear downstream Trajectory and CellChat results to ensure data consistency.*

---

## Code Architecture

The app is built around a modular architecture where a central reactive store manages data flow across independent analysis modules:

```text
BioShiny/
├── app.R            # UI layout, global reactive `store`, module wiring
├── global.R         # Package loading and app configurations
├── R/
│   ├── utils.R      # Shared helper functions
│   └── mod_*.R      # Individual Shiny modules (UI + server logic per tab)
├── inst/
│   └── report.Rmd   # RMarkdown template for HTML report generation
└── install.R        # Dependency setup script
```

### Extending BioShiny
Every module shares the central `store` (`reactiveValues` object). To add a new analysis module:
1. Create `R/mod_new.R` containing `mod_new_ui()` and `mod_new_server(id, store)`.
2. Register the UI tab in `app.R` under `menuItem` and `tabItem`.
3. Call `mod_new_server("new_id", store)` in the main server function in `app.R`.

---

## Limitations

- **Two-group Comparisons**: Differential expression is currently restricted to two-level comparisons (one reference vs. one test group per run).
- **Species & Enrichment**: Functional enrichment natively supports Human and Mouse gene identifiers. Active internet access is required for KEGG queries.
- **Single-Cell Processing**: Currently supports one sample per scRNA-seq run (multi-sample integration/Harmony is not yet supported). CellChat is installed from GitHub; the rest of the app works without it.
- **Trajectory Gene Ranking**: The "genes changing along pseudotime" table uses Spearman correlation ranking rather than dynamic statistical testing (e.g., tradeSeq).

---

## Author / Contact

**Muralidharan Mani, Ph.D.**  
Honorary Fellow  
Department of Biochemistry  
University of Wisconsin–Madison  
Email: [mmani3@wisc.edu](mailto:mmani3@wisc.edu)