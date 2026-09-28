# Changelog

## scSimEval 0.99.9

### Documentation

- **Global search bar added to GitHub Pages site:**
  - Enabled pkgdown’s built-in Fuse.js client-side full-text search
    across all reference pages, vignettes, and news entries via
    `search: exclude: ['news/index.html']` in `_pkgdown.yml`.
  - Users can now instantly search any function name, parameter, or
    concept from the documentation navbar.
- **Shiny Studio vignette completely rewritten
  (`vignettes/shiny-app.Rmd`):**
  - Previous version was 72 lines with minimal coverage; new guide is
    comprehensive.
  - Documents all **6 navigation tabs** (Home, Data Hub, Comparative
    Bubble Matrix, Visualizations, Download Results, Help & Getting
    Started) and the **8 Visualization sub-panels** (Evaluation Summary,
    Distribution QC, Scalability Benchmark, Metric Boxplots, Metric
    Heatmap, PCA Ordination, MDS Metric Space, Cell Embeddings).
  - Every UI control is documented with its options, range, and default
    value.
  - All 4 Data Hub evaluation modes (Demo, Single-Cell, Multiomics,
    Upload Saved) are explained in detail including Paired vs. Unpaired
    multiomics distinction.
  - Every export button (JPEG 600 DPI, PDF, Excel, CSV, RDS, zip
    archive) is listed.
  - Cross-reference table maps all Shiny controls to their R API
    equivalents.
- **Shiny Studio added to navbar:**
  - `articles/shiny-app.html` now appears as a top-level **Shiny
    Studio** link in the navbar.
  - Added `articles:` section to `_pkgdown.yml` to explicitly index and
    organize all vignettes.

## scSimEval 0.99.8

### Bug Fixes

- **[`plot_metric_pca()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_pca.md)
  and
  [`plot_metric_mds()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_mds.md)
  — `category = "all"` support:**
  - Fixed a category filter edge case where passing `category = "all"`
    would return an empty data frame (no category named “all” exists)
    instead of showing all available categories.
  - Both functions now treat `category = "all"` or `NULL` identically —
    no category filter is applied, so all metrics from every category
    contribute to the PCA / MDS ordination.
  - This matches the expected behavior when users select “All” in the
    interactive Shiny Studio.

## scSimEval 0.99.7

### Documentation & Online Integration

- **Interactive Shiny Studio Documentation Portal:**
  - Integrated official documentation links
    (`https://kabilanbio.github.io/scSimEval`) into Tab 6 (“Help &
    Getting Started”) of the Shiny application via a top hero callout
    badge and a dedicated gradient footer card with quick-launch
    actions.
  - Linked users directly to online vignettes, tutorials, paired
    vs. unpaired multiomics protocols, and function reference manuals.
- **GitHub Pages Site Synchronization (`pkgdown`):**
  - Updated `_pkgdown.yml` and rebuilt documentation website to align
    with all new features, including
    [`compute_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/compute_dataset_embeddings.md),
    [`plot_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/plot_dataset_embeddings.md),
    [`compute_embedding_quality_metrics()`](https://kabilanbio.github.io/scSimEval/reference/compute_embedding_quality_metrics.md),
    and
    [`extract_dataset_summary()`](https://kabilanbio.github.io/scSimEval/reference/extract_dataset_summary.md).
  - Re-rendered full documentation site in `docs/` reflecting all 62
    metrics, multiomics compatibility guides, and low-dimensional cell
    embedding workflows.

## scSimEval 0.99.6

### Bug Fixes & Standalone Robustness

- **Self-Contained Shiny Studio Runtime:**
  - Resolved `could not find function "compute_dataset_embeddings"`
    error when launching the Shiny app under pre-existing library
    installations by embedding fully self-contained runtime definitions
    for
    [`compute_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/compute_dataset_embeddings.md),
    [`plot_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/plot_dataset_embeddings.md),
    and
    [`compute_embedding_quality_metrics()`](https://kabilanbio.github.io/scSimEval/reference/compute_embedding_quality_metrics.md)
    directly within `inst/shiny/scSimEvalApp/app.R`.
  - Guaranteed seamless rendering of Sub-panel 8 Cell Embeddings for
    demo datasets and custom uploaded matrices alike, without requiring
    manual package re-installation.
  - Updated Tab 6 (“Help & Getting Started”) in the Shiny app studio to
    fully document Sub-panel 8: Cell Embeddings (t-SNE & UMAP), layout
    modes, and quantitative quality metrics table.

### Documentation Across Package & App

- **Comprehensive Package-Wide Documentation Overhaul:**
  - Synchronized `DESCRIPTION`, `README.md`, `NEWS.md`,
    `vignettes/scSimEval-workflow.Rmd`, and
    `vignettes/demo-multiomics.Rmd` across all recent additions.
  - Documented low-dimensional cell embedding workflows, multi-simulator
    comparison layouts, and embedding quality metrics (Silhouette, ARI,
    library size deviation).
  - Documented Dataset Properties Summary extraction
    ([`extract_dataset_summary()`](https://kabilanbio.github.io/scSimEval/reference/extract_dataset_summary.md))
    and multiomics pairing mode selection (paired vs. unpaired).

## scSimEval 0.99.5

### New Features & Enhancements

- **Sub-panel 8: Cell Embeddings (t-SNE & UMAP) & Quantitative Quality
  Metrics:**
  - Implemented high-level dimensionality reduction and visualization
    functions:
    [`compute_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/compute_dataset_embeddings.md),
    [`plot_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/plot_dataset_embeddings.md),
    and
    [`compute_embedding_quality_metrics()`](https://kabilanbio.github.io/scSimEval/reference/compute_embedding_quality_metrics.md).
  - Supports UMAP, t-SNE, and PCA reductions computed simultaneously
    across the empirical reference and all evaluated simulators.
  - Interactive layout options: Multi-Simulator Faceted Grid (comparing
    reference with all simulators side by side) and Direct 1-to-1
    Comparison (Reference vs. Selected Simulator).
  - Rich aesthetic coloring options: Cell-type labels, Unsupervised
    cluster recovery (k-means on principal components), Library size
    gradient, and Detected features gradient.
  - Integrated an interactive summary table of quantitative quality
    metrics directly beneath the embedding plots (Mean Silhouette Score,
    Adjusted Rand Index cluster concordance, Mean Library Size, Mean
    Detected Features, and Library Size Discrepancy %).
  - Full export support: 600 DPI publication JPEG, vectorized PDF,
    inclusion in the multi-page PDF report, and automatic bundling in
    the complete benchmark ZIP download.

## scSimEval 0.99.4

### Documentation & Methodological Guidance

- **Unpaired & Mosaic Multiomics Compatibility Framework:**
  - Expanded `DESCRIPTION`, `vignettes/scSimEval-workflow.Rmd`, and
    `vignettes/demo-multiomics.Rmd` with deep-dive documentation on
    experimental design compatibility.
  - Clarified support for (i) paired multiomics (e.g. 10x Multiome,
    SHARE-seq), (ii) unpaired multiomics (separate cells from the same
    biological system), and (iii) mosaic multiomics (partially
    overlapping cell cohorts or modality subsets).
  - Provided practical code recipes and metric breakdown matrices for
    decomposing mosaic multiomics datasets into paired and unpaired
    evaluation blocks.

## scSimEval 0.99.3

### New Features & Enhancements

- **Paired vs. Unpaired Multiomics Selection:**
  - Added `pairing = c("paired", "unpaired")` parameter (with
    `is_paired` alias) to
    [`evaluate_multiomics_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_accuracy.md)
    and
    [`evaluate_multiple_datasets()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiple_datasets.md).
  - In **Paired mode** (default; e.g. 10x Chromium Multiome, SHARE-seq),
    all 62 measures across all 8 canonical categories are evaluated,
    including cell-level pairing metrics (FOSCTTM, <Match@1>, in silico
    cross-modal generation, direct peak-to-gene regulatory coupling,
    cross-modality correlation).
  - In **Unpaired mode** (separate cells from the same tissue),
    cell-level pairing metrics that require 1-to-1 matching cell
    barcodes are safely omitted with clear informational logging.
    Population-level cross-modal metrics (cross-modal cell-type label
    transfer accuracy and macro-F1, chromatin peak co-accessibility RV
    coefficient, gene co-expression module correlation, and
    accessibility profile concordance) and unimodal metrics for both
    layers are fully computed and standardized.
- **Interactive Shiny Studio App:**
  - Added a dynamic **Multiomics Dataset Type** selector (radio buttons)
    in Tab 2 (Data Hub: Mode 3) allowing users to switch between Paired
    and Unpaired multiomics evaluation with live contextual warning/info
    banners.
  - Automatically routes datasets to the appropriate evaluation pipeline
    and updates benchmark headers and method names accordingly.
  - Updated Help & Manual (Tab 5) with full metric compatibility tables
    and guidelines for paired vs. unpaired multiomics.

## scSimEval 0.99.2

### Bug Fixes

- - Fixed the Bioconductor staging build failure (\`Error in
    x\[clustering == i, \]:

    (subscript) logical subscript too
    long`): corrected all roxygen`@examples`blocks in the clustering/batch/DE/trajectory/multiomics metric functions so cluster labels align with matrix dimensions and all required arguments are supplied. All examples now pass`R
    CMD check\`.

- Removed the example from the internal helper `calc_expected_mi`, which
  errored in an installed-namespace context (“could not find function”).

- Cleaned up Rd generation: `\figure` widths now declared in pixels,
  explicit `@title` tags for distribution metrics, ASCII-only
  documentation text, and regenerated `man/` pages.

- Labeled all vignette code chunks and added
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) blocks to
  every vignette as required by BiocCheck.

### Documentation

- **Paired / Unpaired / Mosaic Multiomics:** Added comprehensive,
  scientifically accurate documentation throughout the package
  clarifying which evaluation metrics are valid for each multiomics
  experimental design:
  - **Paired** (10x Multiome, SHARE-seq, SNARE-seq): Full 62-metric
    evaluation including all Category 7 cross-modal coupling metrics
    (FOSCTTM, <Match@1>, cross-modal generation fidelity, peak-to-gene
    linkage).
  - **Unpaired** (independent scRNA-seq + scATAC-seq from separate
    cells): Full unimodal evaluation (Categories 1-6, 8) plus
    population-level cross-modal metrics (co-expression module fidelity,
    network Jaccard, cross-modal label transfer). Cell-pairing metrics
    (FOSCTTM, cross-modal generation) are not applicable and excluded.
  - **Mosaic** (partial co-measurement, e.g. DOGMA-seq): Per-modality
    unimodal evaluation on the full matrix; paired coupling metrics
    applied only to the co-assayed cell subset.
- Updated `DESCRIPTION`, `README.md`, `vignettes/demo-multiomics.Rmd`,
  and `vignettes/scSimEval-workflow.Rmd` to reflect these distinctions.

## scSimEval 0.99.1

### Bug Fixes

- Fixed R-universe / Bioconductor build failure: removed `^doc$` and
  `^inst/doc$` from `.Rbuildignore` so that pre-rendered vignette HTML
  files are correctly included in the source tarball and
  `R CMD INSTALL --html` succeeds.
- Cleaned up `docs/` folder by removing scratch files (`git_help.Rmd`),
  redundant source Markdown files, and the empty `tutorials/` directory.

## scSimEval 0.99.0

### New Features

- Initial submission to Bioconductor.
- Unified benchmarking and fidelity assessment for single-cell and
  multiomics simulation algorithms (scRNA-seq, scATAC-seq, and paired
  multiomics).
- 62 quantitative fidelity metrics spanning 8 biological and
  computational categories:
  - 1.  Distributional Properties
  - 2.  Correlations & Zero-Inflation
  - 3.  Cellular Structure & Concordance
  - 4.  Batch Effects & Confounder Mixing
  - 22. Biological Signal & Downstream Fidelity
  - 6.  Trajectory & Lineage Dynamics
  - 7.  Cross-Modal Coupling & Modularity
  - 8.  Computational Scalability
- Comprehensive visualization suite with publication-ready comparative
  bubble matrices, individual metric barplots, heatmaps, PCA, and MDS
  projections.
- Interactive Shiny Studio
  ([`launch_scSimEval_app()`](https://kabilanbio.github.io/scSimEval/reference/launch_scSimEval_app.md))
  for real-time benchmark exploration, custom data upload, and 600 DPI
  publication exports.
