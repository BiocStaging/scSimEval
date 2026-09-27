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
