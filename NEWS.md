# scSimEval 0.99.4

## Bug Fixes & Improvements

* **`compute_method_leaderboard()` & `.ingest_bubble_data()` Robustness:**
  - Added automatic fallback to `Method = "Simulation"` and `Category = "Uncategorized"` when evaluating single simulation accuracy results (`evaluate_simulation_accuracy()`) or pre-formed benchmark data frames lacking explicit method labels.
  - Resolved `R CMD check` example failure in `man/compute_method_leaderboard.Rd`.

* **Dependency & Package Best Practices:**
  - Added `irlba` to `Suggests:` in `DESCRIPTION` to properly declare conditional SVD acceleration in `compute_dataset_embeddings()`.
  - Replaced direct `set.seed()` invocation with RNG-safe execution that preserves and restores `.Random.seed` on function exit in accordance with Bioconductor guidelines.
  - Removed global `eval = FALSE` in `vignettes/shiny-app.Rmd` to satisfy Bioconductor vignette check requirements.


# scSimEval 0.99.3

## New Features & Enhancements

* **Paired, Unpaired & Mosaic Multiomics Framework:**
  - Added explicit `pairing = c("paired", "unpaired")` parameter (with `is_paired` alias) to `evaluate_multiomics_accuracy()` and `evaluate_multiple_datasets()`.
  - In **Paired mode** (default; e.g. 10x Chromium Multiome, SHARE-seq), all 62 measures across all 8 canonical categories are evaluated, including cell-level pairing metrics (FOSCTTM, Match@1, in silico cross-modal generation, direct peak-to-gene linkage).
  - In **Unpaired mode** (separate cells from the same tissue/condition), cell-level pairing metrics that require 1-to-1 matching cell barcodes are safely omitted with clear informational logging, while population-level cross-modal metrics and unimodal evaluations are fully computed and standardized.
  - Added comprehensive methodological guidance for decomposing mosaic multiomics datasets into paired and unpaired evaluation blocks.
  - Integrated dynamic **Multiomics Dataset Type** selector (radio buttons) in Shiny Studio Tab 2 (Data Hub: Mode 3) with live contextual warning/info banners.

* **Sub-panel 8: Cell Embeddings (UMAP, t-SNE, PCA) & Quality Metrics:**
  - Implemented high-level dimensionality reduction and visualization functions: `compute_dataset_embeddings()`, `plot_dataset_embeddings()`, and `compute_embedding_quality_metrics()`.
  - Interactive layout options: Multi-Simulator Faceted Grid and Direct 1-to-1 Comparison (Reference vs. Selected Simulator).
  - Flexible cell coloring: Cell Type labels, Unsupervised cluster recovery, Library size gradient, and Detected features gradient.
  - Integrated interactive summary table of quantitative embedding fidelity metrics (Mean Silhouette width, Adjusted Rand Index cluster concordance, and library size discrepancy %).

* **Comprehensive All-in-One Benchmark Archive (.zip) & Export Suite:**
  - Overhauled the complete benchmark bundle (`download_complete_zip`) in Shiny Studio Tab 5 ("Download Results") to deliver all quantitative results, spreadsheets, reports, and figures in an organized hierarchy.
  - **All Metrics in Multiple Formats:** Bundles complete evaluation results in Excel (`.xlsx`), CSV (`.csv`), tab-delimited text (`.txt`), and native R object (`.rds`) formats in both root and `metrics/` subfolders.
  - **Rankings & Metadata Tables:** Automatically generates and packages `method_rankings_leaderboard` (`.csv` & `.txt`) and `dataset_properties_summary` (`.csv` & `.txt`).
  - **Multi-Page Compiled PDF Report:** Compiles all 9 diagnostic and comparative figures into a single publication-quality vector PDF (`scSimEval_all_plots_report.pdf`).
  - **All Kinds of Figures (Grouped & Individual):**
    - `figures/grouped/`: Both 600 DPI publication-grade JPEGs and vectorized PDFs for comparative bubble matrix, overall evaluation summary, scalability benchmark, metric boxplots by category, performance heatmap, PCA simulator ordination, MDS metric space, comparative distribution QC, and cell embeddings comparison grid.
    - `figures/individual/`: Individual publication-grade barplots for every evaluated metric (`metric_bar_<MetricName>.jpeg` & `.pdf`) and simulator-specific 1-to-1 comparison cell embeddings against empirical reference (`cell_embeddings_compare_<SimulatorName>.jpeg` & `.pdf`).
  - Added dedicated **Download TXT Table (.txt)** action button to Card 2 in Tab 5 for one-click tab-delimited exports.

* **Documentation, Search & Site Integration:**
  - Enabled client-side Fuse.js full-text search across the GitHub Pages documentation site, allowing instant search of all functions, parameters, vignettes, and metric descriptions.
  - Completely rewritten `vignettes/shiny-app.Rmd` documenting all 6 navigation tabs, 8 visualization panels, UI controls, export options, and an R-API cross-reference table.
  - Integrated official documentation links and quick-launch actions into Tab 6 ("Help & Getting Started").

* **Dedicated Individual Figure Export (Shiny Tab 5 - Card 4):**
  - Added on-demand figure export suite in Tab 5 ("Download Results") allowing researchers to download any single diagnostic or comparative visualization in either publication-grade vector PDF (`.pdf`) or ultra-high-resolution 600 DPI JPEG (`.jpeg`).
  - Covers all 14 figure types: Comparative Bubble Matrix, Overall Evaluation Summary, Scalability Benchmark, Metric Boxplots, Performance Heatmap, PCA Ordination, MDS Metric Space, Comparative Distribution QC, UMAP Grid, t-SNE Grid, PCA Grid, and 1-to-1 simulator vs reference comparisons.
  - Dynamically populated simulator selector for custom 1-to-1 side-by-side comparison figure downloads.

* **Expanded Embeddings in Complete Benchmark Archive (`download_complete_zip`):**
  - Upgraded the comprehensive `.zip` bundle to systematically generate and organize cell embedding figures for all three supported dimensionality reduction techniques (UMAP, t-SNE, and PCA) in `figures/grouped/` (`09_cell_embeddings_umap_grid`, `10_cell_embeddings_tsne_grid`, `11_cell_embeddings_pca_grid`).
  - Systematically exports 1-to-1 comparison plots for every simulator against the empirical reference across UMAP, t-SNE, and PCA to `figures/individual/`.

* **Documentation & Vignette Synchronization:**
  - Fully updated `vignettes/shiny-app.Rmd` to document the new Card 4 export controls and the expanded zip bundle directory structure.

## Bug Fixes

* **`plot_metric_pca()` and `plot_metric_mds()` — `category = "all"` support:**
  - Fixed category filter edge case where passing `category = "all"` or `NULL` now correctly retains all metrics across all categories without filtering out rows.
* **Self-Contained Shiny Studio Runtime:**
  - Embedded runtime definitions directly within `inst/shiny/scSimEvalApp/app.R` to ensure seamless execution across any R environment without requiring package re-installation.


# scSimEval 0.99.2

## Bug Fixes
* Fixed the Bioconductor staging build failure (`Error in x[clustering == i, ]
  : (subscript) logical subscript too long`): corrected all roxygen `@examples`
  blocks in the clustering/batch/DE/trajectory/multiomics metric functions so
  cluster labels align with matrix dimensions and all required arguments are
  supplied. All examples now pass `R CMD check`.
* Removed the example from the internal helper `calc_expected_mi`, which errored
  in an installed-namespace context ("could not find function").
* Cleaned up Rd generation: `\figure` widths now declared in pixels, explicit
  `@title` tags for distribution metrics, ASCII-only documentation text, and
  regenerated `man/` pages.
* Labeled all vignette code chunks and added `sessionInfo()` blocks to every
  vignette as required by BiocCheck.

## Documentation

* **Paired / Unpaired / Mosaic Multiomics:** Added comprehensive, scientifically
  accurate documentation throughout the package clarifying which evaluation metrics
  are valid for each multiomics experimental design:
  - **Paired** (10x Multiome, SHARE-seq, SNARE-seq): Full 62-metric evaluation
    including all Category 7 cross-modal coupling metrics (FOSCTTM, Match@1,
    cross-modal generation fidelity, peak-to-gene linkage).
  - **Unpaired** (independent scRNA-seq + scATAC-seq from separate cells): Full
    unimodal evaluation (Categories 1-6, 8) plus population-level cross-modal
    metrics (co-expression module fidelity, network Jaccard, cross-modal label
    transfer). Cell-pairing metrics (FOSCTTM, cross-modal generation) are not
    applicable and excluded.
  - **Mosaic** (partial co-measurement, e.g. DOGMA-seq): Per-modality unimodal
    evaluation on the full matrix; paired coupling metrics applied only to the
    co-assayed cell subset.
* Updated `DESCRIPTION`, `README.md`, `vignettes/demo-multiomics.Rmd`, and
  `vignettes/scSimEval-workflow.Rmd` to reflect these distinctions.

# scSimEval 0.99.1

## Bug Fixes
* Fixed R-universe / Bioconductor build failure: removed `^doc$` and `^inst/doc$`
  from `.Rbuildignore` so that pre-rendered vignette HTML files are correctly
  included in the source tarball and `R CMD INSTALL --html` succeeds.
* Cleaned up `docs/` folder by removing scratch files (`git_help.Rmd`),
  redundant source Markdown files, and the empty `tutorials/` directory.

# scSimEval 0.99.0

## New Features
* Initial submission to Bioconductor.
* Unified benchmarking and fidelity assessment for single-cell and multiomics simulation algorithms (scRNA-seq, scATAC-seq, and paired multiomics).
* 62 quantitative fidelity metrics spanning 8 biological and computational categories:
  - (I) Distributional Properties
  - (II) Correlations & Zero-Inflation
  - (III) Cellular Structure & Concordance
  - (IV) Batch Effects & Confounder Mixing
  - (V) Biological Signal & Downstream Fidelity
  - (VI) Trajectory & Lineage Dynamics
  - (VII) Cross-Modal Coupling & Modularity
  - (VIII) Computational Scalability
* Comprehensive visualization suite with publication-ready comparative bubble matrices, individual metric barplots, heatmaps, PCA, and MDS projections.
* Interactive Shiny Studio (`launch_scSimEval_app()`) for real-time benchmark exploration, custom data upload, and 600 DPI publication exports.
