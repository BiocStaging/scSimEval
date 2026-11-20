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

\
`if`` ``(``!`[`requireNamespace`](https://rdrr.io/r/base/ns-load.html)`(``"BiocManager"``, quietly ``=`` ``TRUE``)``)`` ``{`\
`    `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"BiocManager"``)`\
`}`\
`BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"scSimEval"``)`

To install the development version from GitHub:

\
`if`` ``(``!`[`requireNamespace`](https://rdrr.io/r/base/ns-load.html)`(``"remotes"``, quietly ``=`` ``TRUE``)``)`` ``{`\
`    `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"remotes"``)`\
`}`\
`remotes``::`[`install_github`](https://remotes.r-lib.org/reference/install_github.html)`(``"kabilanbio/scSimEval"``)`

------------------------------------------------------------------------

## Quick Start

### 1. Reproducible R Workflow (SingleCellExperiment)

Run benchmarking directly on `SingleCellExperiment` objects with
cell-type and batch annotations stored in `colData`:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`scSimEval`](https://kabilanbio.github.io/scSimEval/)`)`\
\
`# Users starting from Seurat objects should convert first:`\
`# ref_sce <- Seurat::as.SingleCellExperiment(ref_seurat)`\
\
`# Run benchmark evaluation (cell_type & batch read from colData)`\
`benchmark_results`` ``<-`` `[`evaluate_simulation_accuracy`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md)`(`\
`  ref_sce       ``=`` ``ref_sce``,`\
`  sim_sce       ``=`` ``sim_sce``,`\
`  cell_type_col ``=`` ``"cell_type"``,`\
`  batch_col     ``=`` ``"batch"``,`\
`  threads       ``=`` ``4`\
`)`\
\
`# Visualize summary`\
[`plot_evaluation_summary`](https://kabilanbio.github.io/scSimEval/reference/plot_evaluation_summary.md)`(``benchmark_results``)`

### 2. Interactive Shiny Studio

Interactive exploration, point-and-click benchmarking, and report export
are available through the companion **`scSimEvalStudio`** package:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`scSimEvalStudio`](https://kabilanbio.github.io/scSimEvalStudio/)`)`\
[`launch_scSimEval_app`](https://kabilanbio.github.io/scSimEvalStudio/reference/launch_scSimEval_app.html)`(``)`

------------------------------------------------------------------------

## Documentation

Full documentation, tutorials, and demonstration workflows are available
at:\
👉 **<https://kabilanbio.github.io/scSimEval/>**

------------------------------------------------------------------------

## Development & Support

`scSimEval` is developed and maintained at the **Division of
Agricultural Bioinformatics, ICAR-IASRI, New Delhi**.

- **Maintainer:** Kabilan Sakthivel (`kabilan151414@gmail.com`)
- **Bug Reports & Issues:**
  <https://github.com/kabilanbio/scSimEval/issues>
