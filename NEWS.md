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
