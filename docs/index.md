# scSimEval: Unified Benchmarking for Single-Cell Multiomics Simulations

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://kabilanbio.github.io/scSimEval/)

------------------------------------------------------------------------

## Overview

**`scSimEval`** is a unified evaluation and benchmarking toolkit
designed to systematically quantify how faithfully single-cell and
multiomics simulation methods reproduce empirical biological data across
**56 curated evaluation measures** spanning **8 foundational
categories**.

### Supported Modalities

- **`scRNA-seq`** (Single-cell transcriptomics)
- **`scATAC-seq`** (Single-cell chromatin accessibility)
- **`Paired Multiomics`** (Joint scRNA-seq + scATAC-seq across identical
  cells)
- **`Unpaired Multiomics`** (Independent scRNA-seq and scATAC-seq from
  the same biological system)

------------------------------------------------------------------------

## Installation

Install the package from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("scSimEval")
```

To install the development version from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("kabilanbio/scSimEval")
```

------------------------------------------------------------------------

## Quick Start

### 1. Reproducible R Workflow (SingleCellExperiment)

Run benchmarking directly on `SingleCellExperiment` objects with
cell-type and batch annotations stored in `colData`:

``` r
library(scSimEval)

# Users starting from Seurat objects should convert first:
# ref_sce <- Seurat::as.SingleCellExperiment(ref_seurat)

# Run benchmark evaluation (cell_type & batch read from colData)
benchmark_results <- evaluate_simulation_accuracy(
  ref_sce       = ref_sce,
  sim_sce       = sim_sce,
  cell_type_col = "cell_type",
  batch_col     = "batch",
  threads       = 4
)

# Visualize summary
plot_evaluation_summary(benchmark_results)
```

### 2. Interactive Shiny Studio

Interactive exploration, point-and-click benchmarking, and report export
are available through the companion **`scSimEvalStudio`** package:

``` r
library(scSimEvalStudio)
launch_scSimEval_app()
```

------------------------------------------------------------------------

## Documentation

Full documentation, tutorials, and demonstration workflows are available
at:  
👉 **<https://kabilanbio.github.io/scSimEval/>**

------------------------------------------------------------------------

## Development & Support

`scSimEval` is developed and maintained at the **Division of
Agricultural Bioinformatics, ICAR-IASRI, New Delhi**.

- **Maintainer:** Kabilan Sakthivel (`kabilan151414@gmail.com`)
- **Bug Reports & Issues:**
  <https://github.com/kabilanbio/scSimEval/issues>
