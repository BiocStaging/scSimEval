# scSimEval: Unified Benchmarking for Single-Cell Multiomics Simulations <a href="https://kabilanbio.github.io/scSimEval/"><img src="man/figures/logo.png" align="right" height="135" alt="scSimEval logo" /></a>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://kabilanbio.github.io/scSimEval/)

---

## Overview

**`scSimEval`** is a unified evaluation and benchmarking toolkit designed to systematically quantify how faithfully single-cell and multiomics simulation methods reproduce empirical biological data across **56 curated evaluation measures** spanning **8 foundational categories**.

### Supported Modalities
* **`scRNA-seq`** (Single-cell transcriptomics)
* **`scATAC-seq`** (Single-cell chromatin accessibility)
* **`Paired Multiomics`** (Joint scRNA-seq + scATAC-seq across identical cells)
* **`Unpaired Multiomics`** (Independent scRNA-seq and scATAC-seq from the same biological system)

---

## Installation

Install the development version of **`scSimEval`** directly from GitHub:

```R
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("kabilanbio/scSimEval")
```

---

## Quick Start

### 1. Interactive Shiny Studio
Launch the built-in, interactive graphical studio directly from R:

```R
library(scSimEval)
launch_scSimEval_app()
```

The Shiny app provides automated dataset validation, real-time metric calculation, interactive diagnostic plots with zoom/scroll, and summary report exports.

### 2. Command-Line R Workflow
Run benchmarking directly in reproducible R scripts or HPC pipelines:

```R
library(scSimEval)

# Run benchmark evaluation
benchmark_results <- evaluate_simulation_accuracy(
  ref_data   = reference_counts,
  sim_data   = simulated_counts,
  cell_types = cell_type_annotations,
  batch_info = batch_annotations,
  threads    = 4
)

# Visualize summary
plot_evaluation_summary(benchmark_results)
```

---

## Documentation

Full documentation, tutorials, and demonstration workflows are available at:  
👉 **[https://kabilanbio.github.io/scSimEval/](https://kabilanbio.github.io/scSimEval/)**

---

## Development & Support

`scSimEval` is developed and maintained at the **Division of Agricultural Bioinformatics, ICAR-IASRI, New Delhi**.

- **Maintainer:** Kabilan Sakthivel (`kabilan151414@gmail.com`)
- **Bug Reports & Issues:** [https://github.com/kabilanbio/scSimEval/issues](https://github.com/kabilanbio/scSimEval/issues)
