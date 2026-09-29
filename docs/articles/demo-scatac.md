# Demonstration 2: scATAC-seq Simulation Evaluation

## Introduction

Single-cell assay for transposase-accessible chromatin sequencing
(scATAC-seq) measures chromatin accessibility at single-cell resolution.
Unlike scRNA-seq expression counts, scATAC-seq data have unique
properties: \* **High Sparsity:** Over 90% to 98% of values in
peak-by-cell matrices are zeroes because each locus has at most 2 copies
per cell. \* **Near-Binary Nature:** Values primarily represent open
(accessible) vs. closed chromatin states. \* **Peak Co-Accessibility:**
Genomic loci regulated together exhibit correlated accessibility across
cells.

This tutorial shows how to evaluate synthetic scATAC-seq data using
**`scSimEval`**.

------------------------------------------------------------------------

## 1. Loading the Packaged scATAC-seq Dataset

`scSimEval` includes an empirical reference and simulated scATAC-seq
peak accessibility dataset (`example_scatac`):

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`scSimEval`](https://kabilanbio.github.io/scSimEval)`)`\
\
`# Load packaged scATAC-seq example data`\
[`data`](https://rdrr.io/r/utils/data.html)`(``example_scatac``)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Reference scATAC-seq dimensions:"``, `[`dim`](https://rdrr.io/r/base/dim.html)`(``example_scatac``$``ref``)``, ``"\n"``)`\
`#> Reference scATAC-seq dimensions: 60 80`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Simulated scATAC-seq dimensions:"``, `[`dim`](https://rdrr.io/r/base/dim.html)`(``example_scatac``$``sim``)``, ``"\n"``)`\
`#> Simulated scATAC-seq dimensions: 60 80`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Cell types:"``, `[`levels`](https://rdrr.io/r/base/levels.html)`(``example_scatac``$``cell_types``)``, ``"\n"``)`\
`#> Cell types: TypeA TypeB`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Batches:"``, `[`levels`](https://rdrr.io/r/base/levels.html)`(``example_scatac``$``batch_info``)``, ``"\n"``)`\
`#> Batches: Batch1 Batch2`

The matrices contain integer insertion counts across genomic peaks
($`60\text{ peaks} \times 80\text{ cells}`$).

------------------------------------------------------------------------

## 2. Epigenomic Distribution & Summary Properties

We evaluate whether simulated peak insertion frequencies and cell
library sizes match real distributions using
[`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md):

\
`atac_res`` ``<-`` `[`evaluate_simulation_accuracy`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md)`(`\
`  ref_data          ``=`` ``example_scatac``$``ref``,`\
`  sim_data          ``=`` ``example_scatac``$``sim``,`\
`  compute_bivariate ``=`` ``FALSE``,`\
`  verbose           ``=`` ``FALSE`\
`)`\
\
`# Inspect peak accessibility and library size statistical distances`\
[`head`](https://rdrr.io/r/utils/head.html)`(``atac_res``$``metrics_summary_table``, ``8``)`\
`#>                    Category     Property        Metric        Value`\
`#> 1 Distributional Properties library_size           MAD 2.0000000000`\
`#> 2 Distributional Properties library_size            KS 0.1750000000`\
`#> 3 Distributional Properties library_size           MAE 2.3375000000`\
`#> 4 Distributional Properties library_size          RMSE 2.4874685928`\
`#> 5 Distributional Properties library_size            OV 0.8346458491`\
`#> 6 Distributional Properties library_size Bhattacharyya 0.0001055641`\
`#> 7 Distributional Properties library_size   Wasserstein 2.3375000000`\
`#> 8 Distributional Properties library_size ECDF_DiffArea 0.0863425926`

------------------------------------------------------------------------

## 3. Epigenomic Cluster & Cell-Type Concordance

In scATAC-seq, identifying cell types requires clustering cells based on
chromatin accessibility profiles. We test whether the simulator
preserves genuine biological separation using
[`evaluate_clustering_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_clustering_metrics.md):

\
`atac_clust`` ``<-`` `[`evaluate_clustering_metrics`](https://kabilanbio.github.io/scSimEval/reference/evaluate_clustering_metrics.md)`(`\
`  data       ``=`` ``example_scatac``$``sim``,`\
`  cell_types ``=`` ``example_scatac``$``cell_types``,`\
`  ref_data   ``=`` ``example_scatac``$``ref`\
`)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Average Silhouette Width (Simulated):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_clust``$``silhouette``, ``4``)``, ``"\n"``)`\
`#> Average Silhouette Width (Simulated): 0.0465`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Davies-Bouldin Index:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_clust``$``davies_bouldin``, ``4``)``, ``"\n"``)`\
`#> Davies-Bouldin Index: 3.9852`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Adjusted Rand Index (ARI):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_clust``$``ARI``, ``4``)``, ``"\n"``)`\
`#> Adjusted Rand Index (ARI): 0.7627`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Normalized Mutual Information (NMI):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_clust``$``NMI``, ``4``)``, ``"\n"``)`\
`#> Normalized Mutual Information (NMI): 0.721`

------------------------------------------------------------------------

## 4. Epigenomic Batch Integration & Confounder Handling

When simulating multi-sample experiments, batch variations can confound
biological signals. We evaluate batch mixing metrics using
[`evaluate_batch_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_batch_metrics.md):

\
`atac_batch`` ``<-`` `[`evaluate_batch_metrics`](https://kabilanbio.github.io/scSimEval/reference/evaluate_batch_metrics.md)`(`\
`  data       ``=`` ``example_scatac``$``sim``,`\
`  batch_info ``=`` ``example_scatac``$``batch_info``,`\
`  cell_types ``=`` ``example_scatac``$``cell_types``,`\
`  verbose    ``=`` ``FALSE`\
`)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Batch Shannon Entropy:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_batch``$``shannon_entropy``, ``4``)``, ``"\n"``)`\
`#> Batch Shannon Entropy: 0.9671`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Principal Component Regression (PCR R2):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_batch``$``pcr_r2``, ``4``)``, ``"\n"``)`\
`#> Principal Component Regression (PCR R2): 0.0125`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Cross-Batch Transfer Accuracy:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_batch``$``cross_batch_accuracy``, ``4``)``, ``"\n"``)`\
`#> Cross-Batch Transfer Accuracy: 0.825`

------------------------------------------------------------------------

## 5. Chromatin Peak Co-Accessibility & Regulatory Coupling

In real biological cells, peaks that are close to each other or
regulated by the same transcription factor complexes show correlated
accessibility:

\
`# Evaluate peak co-accessibility fidelity`\
`coacc_res`` ``<-`` `[`calc_peak_coaccessibility_fidelity`](https://kabilanbio.github.io/scSimEval/reference/calc_peak_coaccessibility_fidelity.md)`(`\
`  ref_atac ``=`` ``example_scatac``$``ref``,`\
`  sim_atac ``=`` ``example_scatac``$``sim`\
`)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Co-Accessibility RV Coefficient:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``coacc_res``$``rv_coefficient``, ``4``)``, ``"\n"``)`\
`#> Co-Accessibility RV Coefficient: 0.6014`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Peak Co-Accessibility Pearson Correlation:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``coacc_res``$``coaccessibility_pearson``, ``4``)``, ``"\n"``)`\
`#> Peak Co-Accessibility Pearson Correlation: 0.1128`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Peak Co-Accessibility Spearman Correlation:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``coacc_res``$``coaccessibility_spearman``, ``4``)``, ``"\n"``)`\
`#> Peak Co-Accessibility Spearman Correlation: 0.0902`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Frobenius Matrix Distance:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``coacc_res``$``frobenius_distance``, ``4``)``, ``"\n"``)`\
`#> Frobenius Matrix Distance: 0.0053`

------------------------------------------------------------------------

## 6. Master Unimodal scATAC-seq Pipeline

Execute the comprehensive evaluation workflow in a single coordinated
command:

\
`master_atac`` ``<-`` `[`evaluate_simulation_accuracy`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md)`(`\
`  ref_data          ``=`` ``example_scatac``$``ref``,`\
`  sim_data          ``=`` ``example_scatac``$``sim``,`\
`  elapsed_time      ``=`` ``18.2``,`\
`  memory_mb         ``=`` ``512.4``,`\
`  compute_bivariate ``=`` ``FALSE``,`\
`  verbose           ``=`` ``FALSE`\
`)`\
\
`# Inspect tidy summary table`\
[`head`](https://rdrr.io/r/utils/head.html)`(``master_atac``$``metrics_summary_table``)`\
`#>                    Category     Property        Metric        Value`\
`#> 1 Distributional Properties library_size           MAD 2.0000000000`\
`#> 2 Distributional Properties library_size            KS 0.1750000000`\
`#> 3 Distributional Properties library_size           MAE 2.3375000000`\
`#> 4 Distributional Properties library_size          RMSE 2.4874685928`\
`#> 5 Distributional Properties library_size            OV 0.8346458491`\
`#> 6 Distributional Properties library_size Bhattacharyya 0.0001055641`

------------------------------------------------------------------------

## 7. Conclusion

`scSimEval` effectively benchmarks scATAC-seq simulations, validating
peak insertion distributions, extreme sparsity, cell type clustering,
batch mixing, and chromatin peak co-accessibility without requiring
synthetic ground-truth annotations.

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`\
`#> R version 4.6.1 (2026-06-24 ucrt)`\
`#> Platform: x86_64-w64-mingw32/x64`\
`#> Running under: Windows 11 x64 (build 26200)`\
`#> `\
`#> Matrix products: default`\
`#>   LAPACK version 3.12.1`\
`#> `\
`#> locale:`\
`#> [1] LC_COLLATE=English_India.utf8  LC_CTYPE=English_India.utf8   `\
`#> [3] LC_MONETARY=English_India.utf8 LC_NUMERIC=C                  `\
`#> [5] LC_TIME=English_India.utf8    `\
`#> `\
`#> time zone: Asia/Calcutta`\
`#> tzcode source: internal`\
`#> `\
`#> attached base packages:`\
`#> [1] stats     graphics  grDevices utils     datasets  methods   base     `\
`#> `\
`#> other attached packages:`\
`#> [1] scSimEval_0.99.4`\
`#> `\
`#> loaded via a namespace (and not attached):`\
`#>  [1] SummarizedExperiment_1.42.0 xfun_0.60                  `\
`#>  [3] bslib_0.12.0                htmlwidgets_1.6.4          `\
`#>  [5] Biobase_2.73.2              lattice_0.22-9             `\
`#>  [7] tools_4.6.1                 generics_0.1.4             `\
`#>  [9] parallel_4.6.1              clValid_0.7                `\
`#> [11] stats4_4.6.1                flexmix_2.3-21             `\
`#> [13] proxy_0.4-29                DEoptimR_1.2-0             `\
`#> [15] cluster_2.1.8.3             pkgconfig_2.0.3            `\
`#> [17] BiocNeighbors_2.6.0         Matrix_1.7-6               `\
`#> [19] RColorBrewer_1.1-3          desc_1.4.3                 `\
`#> [21] S4Vectors_0.50.1            lifecycle_1.0.5            `\
`#> [23] clusterSim_0.51-6           compiler_4.6.1             `\
`#> [25] textshaping_1.0.5           statmod_1.5.2              `\
`#> [27] bluster_1.22.0              codetools_0.2-20           `\
`#> [29] Seqinfo_1.2.0               clue_0.3-68                `\
`#> [31] htmltools_0.5.9             class_7.3-24               `\
`#> [33] sass_0.4.10                 yaml_2.3.12                `\
`#> [35] pkgdown_2.2.1               jquerylib_0.1.4            `\
`#> [37] prabclus_2.3-5              MASS_7.3-66                `\
`#> [39] BiocParallel_1.47.0         diptest_0.77-2             `\
`#> [41] SingleCellExperiment_1.34.0 DelayedArray_0.38.2        `\
`#> [43] cachem_1.1.0                limma_3.68.4               `\
`#> [45] fpc_2.2-15                  abind_1.4-8                `\
`#> [47] mclust_6.1.3                robustbase_0.99-7          `\
`#> [49] locfit_1.5-9.12             digest_0.6.39              `\
`#> [51] kernlab_0.9-33              splines_4.6.1              `\
`#> [53] ade4_1.7-24                 fastmap_1.2.0              `\
`#> [55] grid_4.6.1                  cli_3.6.6                  `\
`#> [57] SparseArray_1.13.2          magrittr_2.0.5             `\
`#> [59] S4Arrays_1.13.0             e1071_1.7-17               `\
`#> [61] edgeR_4.10.1                rmarkdown_2.31             `\
`#> [63] XVector_0.52.0              matrixStats_1.5.0          `\
`#> [65] igraph_2.3.3                otel_0.2.0                 `\
`#> [67] nnet_7.3-21                 RANN_2.6.2                 `\
`#> [69] ragg_1.5.2                  modeltools_0.2-24          `\
`#> [71] evaluate_1.0.5              knitr_1.51                 `\
`#> [73] GenomicRanges_1.64.0        IRanges_2.46.0             `\
`#> [75] rlang_1.3.0                 Rcpp_1.1.2                 `\
`#> [77] BiocGenerics_0.58.1         jsonlite_2.0.0             `\
`#> [79] R6_2.6.1                    MatrixGenerics_1.24.0      `\
`#> [81] systemfonts_1.3.2           fs_2.1.0`
