# Changelog

## scSimEval 0.99.13

- **Curated User-Focused Reference Documentation**: Streamlined the
  `pkgdown` Reference tab from 14 sprawling sections (~120 topics) down
  to 5 clear sections featuring the 36 essential user-facing workflows,
  scoring/embedding functions, publication visualizations, category
  evaluation suites, and bundled datasets (low-level kernels and
  individual `calc_*` helpers remain exported in `NAMESPACE` and
  documented via `@keywords internal`).
- **Complete, Uncropped Reference Examples & Figures**: Completed all
  truncated visualization examples (`plot_scalability_benchmark`,
  `plot_metric_heatmap`, `plot_metric_boxplots`,
  `plot_category_metric_bars`, `plot_individual_metric_bar`,
  `plot_metric_pca`, `plot_metric_mds`), added concise printed outputs
  across workflow/dataset examples, widened the `pkgdown` figure canvas
  (`11 in` at `150 DPI`), expanded Panel B label margins and added
  `ggrepel` label separation in
  [`plot_scalability_benchmark()`](https://kabilanbio.github.io/scSimEval/reference/plot_scalability_benchmark.md),
  and stacked bottom legends in
  [`plot_benchmark_bubble_matrix()`](https://kabilanbio.github.io/scSimEval/reference/plot_benchmark_bubble_matrix.md)
  examples so zero text or legends are cropped.
- **Modern `edgeR` Normalization (`normLibSizes` + `TMMwsp`)**: Replaced
  deprecated
  [`edgeR::calcNormFactors()`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)
  calls with `edgeR::normLibSizes(..., method = "TMMwsp")` for sparse
  single-cell count matrices, eliminating all deprecation messages and
  sparse-cell [`max()`](https://rdrr.io/r/base/Extremes.html) warnings
  without using
  [`suppressWarnings()`](https://rdrr.io/r/base/warning.html).
- **Proper Mathematical Notation Rendering (KaTeX)**: Enabled
  `math-rendering: katex` in `_pkgdown.yml` and converted remaining raw
  `$...$` roxygen strings to valid `\eqn{}` / `\deqn{}` notation so all
  formulas render cleanly as typeset mathematics across the
  documentation website.

## scSimEval 0.99.12

- **Bioconductor Installation Documentation**: Added standard
  `BiocManager::install("scSimEval")` installation instructions across
  `README.md`, package vignettes, and `pkgdown` documentation.

## scSimEval 0.99.11

- **Package Tarball Size Reduction**: Excluded web-only `pkgdown`
  demonstration articles (`vignettes/demo-*.Rmd`), pre-rendered figure
  assets (`vignettes/figures/`), and pre-saved benchmark `.rds` archives
  (`vignettes/*.rds`) via `.Rbuildignore`, reducing the built
  Bioconductor source tarball size to `< 1 MB`.

## scSimEval 0.99.10

- **Standalone Interactive GUI (`scSimEvalStudio`)**: Decoupled the
  Shiny application into a dedicated companion package
  (`scSimEvalStudio`), maintaining `scSimEval` as a focused, lightweight
  computational benchmarking engine per Bioconductor design principles.
- **Standardized SingleCellExperiment API**: Core evaluation workflows
  (`evaluate_simulation_accuracy`, `evaluate_multiomics_accuracy`,
  `extract_dataset_summary`) now natively require and operate directly
  on standard Bioconductor `SingleCellExperiment` objects with
  standardized `colData` annotations.
- **Bioconductor Multiomics via `altExp`**: Upgraded multiomics
  evaluation pipelines to represent paired multimodal data (scRNA-seq +
  scATAC-seq) using Bioconductor’s standard `altExp()` assay
  representation.
- **Streamlined Dependencies**: Replaced direct dependencies on heavy
  external single-cell ecosystems (`Seurat`, `Signac`) by moving
  optional adapters to `Suggests`, drastically reducing dependency
  overhead.
- **Accessor Method Compliance**: Eliminated all direct S4 `@` slot
  accesses across metric computation, trajectories, and dimension
  reductions in favor of formal accessor methods (`assay()`, `assays()`,
  `colData()`, `reducedDim()`, `altExp()`).
- **Clean Condition & Error Handling**: Removed all instances of
  [`suppressWarnings()`](https://rdrr.io/r/base/warning.html),
  [`suppressMessages()`](https://rdrr.io/r/base/message.html),
  [`unlockBinding()`](https://rdrr.io/r/base/bindenv.html), and global
  variable assignments (`<<-`) to adhere to strict Bioconductor code
  hygiene standards.
- **Modular Function Architecture**: Refactored monolithic visualization
  and benchmarking functions into cohesive, single-purpose helper
  functions for improved maintainability and extensibility.
- **Comprehensive Return Documentation**: Upgraded `@return`
  documentation across all exported functions to explicit `\describe{}`
  format detailing output slots, dimensions, and object structures.
- **Fully Evaluated & Reproducible Vignettes**: Reworked vignettes
  (`scSimEval-workflow.Rmd` and `metrics-catalog.Rmd`) into fully
  executable, reproducible workflows incorporating sample datasets and
  concluding with
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).
- **High-DPI Visual Assets**: Upgraded workflow architecture diagrams
  (`workflow_diagram.png`) to crisp 350 DPI high-resolution figures
  across documentation, vignettes, and web resources.

## scSimEval 0.99.9

- **Category VII Trajectory Standardization**: Upgraded Category VII to
  the `dynverse` trajectory inference protocol (`him`, `edge_flip`,
  `geodesic_correlation`, `F1_branches`, and `wcor_features`).
- **Category VIII Scalability Normalization**: Standardized
  computational efficiency to core-normalized metrics (`core_seconds`
  and `memory_per_core_mib`).
- **S4 Object Standardization**: Standardized input handling to require
  `Seurat` or `SingleCellExperiment` objects with precomputed
  dimensional reductions.
- **Category VI Integration Metrics**: Standardized multiomics
  integration to 5 core metrics (`batch_asw`, `cell_type_asw`,
  `graph_connectivity`, `ari`, and `nmi`).
- **Category V Pseudobulk DE/DA**: Structured Category V into scRNA-seq
  DE and scATAC-seq DA sub-categories using biological replicate
  pseudobulk testing via `edgeR-QLF`.
- **Multi-Dataset Benchmarking**: Added
  [`combine_benchmark_results()`](https://kabilanbio.github.io/scSimEval/reference/combine_benchmark_results.md)
  to aggregate and average performance across multiple experimental
  datasets.
- **Interactive Diagnostic Logs**: Added real-time log streaming and
  diagnostic console in the Shiny application.

## scSimEval 0.99.8

- **Bubble Matrix Visual Refinements**: High-contrast bold circular
  outlines highlight best-performing simulators.
- **Modality Selectors**: Independent assay modality filters (scRNA-seq,
  scATAC-seq, Combined) across visualization tabs.
- **Export Enhancements**: Standalone reproducible R script generators
  for all diagnostic plots.

## scSimEval 0.99.7

- **Input Capacity Guidance**: Context-aware guidance detailing
  available metrics based on uploaded matrix, cell-type, and batch
  annotations.
- **Reproducible Script Downloads**: 1-click generation of
  self-contained ggplot2 R scripts for each visualization panel.

## scSimEval 0.99.6

- **Cross-Dataset Averaging**: Seamless toggling between cross-dataset
  averages and individual datasets across all visualization panels.
- **Memory Standardization**: Standardized computational memory
  reporting to binary MiB RAM (1\text{ MiB} = 1,048,576\text{ bytes}).

## scSimEval 0.99.5

- **Multi-Dataset Benchmark Suite**: Added
  [`benchmark_simulation_suite()`](https://kabilanbio.github.io/scSimEval/reference/benchmark_simulation_suite.md)
  for high-throughput multi-dataset evaluations.
- **Redesigned Data Hub**: Four streamlined workflow modes: Example
  Datasets, Single-Cell Upload, Multiomics Upload, and Saved Results.

## scSimEval 0.99.4

- **Large Matrix Uploads**: Increased Shiny upload limit to 10 GB
  (`shiny.maxRequestSize`).
- **Organism-Agnostic Metrics**: Removed species-specific priors to
  ensure universal applicability to plant and animal datasets.

## scSimEval 0.99.3

- **Paired & Unpaired Multiomics**: Dedicated evaluation workflows for
  paired (co-assay) and unpaired single-cell multiomics.
- **Cell Embeddings**: Added UMAP, t-SNE, and PCA visualization panels
  with quantitative embedding quality metrics.
- **Complete Benchmark Archive**: Multi-format ZIP export containing CSV
  tables, Excel workbooks, RDS objects, vector PDFs, and 600 DPI
  figures.

## scSimEval 0.99.2

- **Multiomics Methodology**: Documented metric applicability across
  paired, unpaired, and mosaic multiomics experimental designs.
- **Package Check Fixes**: Corrected documentation examples and vignette
  chunk settings for Bioconductor compliance.

## scSimEval 0.99.1

- **Documentation & Build Fixes**: Cleaned up documentation files and
  build ignore paths for R-universe / Bioconductor staging.

## scSimEval 0.99.0

- **Initial Submission**: Initial submission of `scSimEval` to
  Bioconductor for benchmarking single-cell and multiomics simulations
  across 8 biological and computational categories.
