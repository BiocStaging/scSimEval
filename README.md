# scSimEval: Unified Benchmarking and Accuracy Evaluation for Single-Cell Multiomics Simulation Methods <a href="https://kabilanbio.github.io/scSimEval/"><img src="man/figures/logo.png" align="right" height="135" alt="scSimEval logo" /></a>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://kabilanbio.github.io/scSimEval/)
[![Bioconductor submission](https://img.shields.io/badge/Bioconductor-submission%20v0.99.10-brightgreen.svg)](https://github.com/BiocStaging/scSimEval)

---

## Overview

**`scSimEval`** is a unified evaluation and benchmarking toolkit designed to systematically quantify how faithfully single-cell and multiomics simulation methods reproduce empirical biological experiments.

The package provides a standardized, reference-based benchmarking framework supporting four primary operational modalities:
* **`scRNA-seq`** (Single-cell transcriptomic count matrices)
* **`scATAC-seq`** (Single-cell chromatin accessibility peak matrices)
* **`Paired Multiomics`** (Simultaneous co-assays across identical cells, e.g., 10x Chromium Multiome, SHARE-seq)
* **`Unpaired Multiomics`** (Independent scRNA-seq and scATAC-seq profiled from separate cells of the same biological system)

---

## Installation

### Prerequisites & System Requirements

Before installing `scSimEval`, please ensure your system meets the following requirements:

1. **R Version:** R (>= 4.0.0, recommended R >= 4.3.0).
2. **C++ Compilation Tools:**
   - **Windows:** Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version (required for compiling C++ extensions such as `Rcpp` and `RANN`).
   - **macOS:** Install Xcode Command Line Tools (`xcode-select --install`) and `gfortran`.
   - **Linux:** Ensure standard development packages (`build-essential`, `g++`, `libcurl4-openssl-dev`, `libxml2-dev`) are installed.
3. **Bioconductor Dependencies:**
   `scSimEval` relies on core Bioconductor packages for single-cell data structures and pseudobulk differential testing. We recommend installing them first via `BiocManager`:
   ```R
   if (!requireNamespace("BiocManager", quietly = TRUE)) {
     install.packages("BiocManager")
   }

   BiocManager::install(c(
     "SingleCellExperiment",
     "SummarizedExperiment",
     "edgeR",
     "BiocNeighbors",
     "BiocParallel"
   ), update = FALSE)
   ```

### Installing `scSimEval` from GitHub

Install the development release of **`scSimEval`** from GitHub using `remotes`:

```R
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("kabilanbio/scSimEval", dependencies = TRUE)
```

### Verification

Verify that `scSimEval` has installed successfully:

```R
library(scSimEval)
packageVersion("scSimEval")
```

---

## Execution Options: Interactive Shiny Studio & Reproducible R Script

`scSimEval` supports two distinct, flexible execution workflows tailored for both interactive exploratory analysis and large-scale computational pipelines:

### Option 1: Interactive Shiny App Studio (`scSimEvalApp`)
For interactive graphical exploration, `scSimEval` provides an all-in-one, no-code web studio directly within R:
* **Interactive Upload & Validation:** Automated file format checking, metadata alignment, and assay detection.
* **Comparative Visualizations:** Live generation of 600 DPI publication-grade plots including Comparative Bubble Matrices, Category Radar Charts, Heatmaps, Boxplots, and Multi-Simulator UMAP/PCA Embeddings.
* **Filter & Leaderboard Exploration:** Real-time simulator rank ordering, metric selection, and Excel/CSV summary report exports.

Launch the Shiny studio directly from your R console:
```R
library(scSimEval)
launch_scSimEval_app()
```

### Option 2: Standalone & Headless R Script (Workstations & HPC Clusters)
For large-scale, automated, or compute-intensive benchmarking on multi-core workstations and High-Performance Computing (HPC) clusters:
* **Download Pre-Configured R Script:** In the Shiny app **Data Hub**, clicking **"Download R Script"** automatically generates a complete, self-contained, and procedural R workflow script pre-populated with your exact dataset paths, metadata columns, and core settings.
* **Headless Command-Line Execution:** Run directly in batch scripts or via R without launching a browser or consuming graphical display resources:
  ```R
  library(scSimEval)

  # Run full benchmark evaluation headlessly
  benchmark_results <- evaluate_simulation_accuracy(
    ref_data   = reference_counts,
    sim_data   = simulated_counts,
    cell_types = cell_type_annotations,
    batch_info = batch_annotations,
    threads    = 4,
    verbose    = TRUE
  )

  # Save results object
  saveRDS(benchmark_results, "scSimEval_benchmark_results.rds")
  ```

---

## Authors

`scSimEval` is developed and maintained at the **Division of Agricultural Bioinformatics, ICAR-Indian Agricultural Statistics Research Institute (ICAR-IASRI)**, New Delhi, India:

- **Kabilan Sakthivel** (Maintainer & Ph.D Scholar) — `kabilan151414@gmail.com`
- **Dr. Dwijesh Chandra Mishra** (Principal Scientist & Thesis Guide) — `dwij.mishra@gmail.com`
- **Dr. Shashi Bhushan Lal** (Principal Scientist) — `sblall16@gmail.com`
- **Dr. Ulavappa Angadi** (Principal Scientist) — `ubaangadi.iasri@icar.org.in`
- **Dr. Sudhir Srivastava** (Senior Scientist) — `sudhir0401bm@gmail.com`
- **Dr. Neeraj Budhlakoti** (Senior Scientist) — `neeraj35669@gmail.com`
- **Dr. Krishna Kumar Chaturvedi** (Principal Scientist) — `kkcchaturvedi@gmail.com`
- **Dr. Ranjeet Ranjan Kumar** (Senior Scientist) — `ranjeetranjaniari@gmail.com`
- **Dr. Navin Chandra Gupta** (Senior Scientist) — `navinbtc@gmail.com`
- **Dr. Sharanbasappa** (Scientist) — `smadival509@gmail.com`

- **Bug Reports & Issues:** [https://github.com/kabilanbio/scSimEval/issues](https://github.com/kabilanbio/scSimEval/issues)
