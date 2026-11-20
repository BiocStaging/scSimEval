# Comprehensive Catalog of the 56 Evaluation Measures

## Overview of the scSimEval Evaluation Suite

Single-cell and multiomics simulation algorithms are central to
computational biology. They let researchers benchmark new bioinformatic
methods, optimize experimental designs, and test analytical workflows
without the expense of wet-lab sequencing. However, for a simulator to
be truly useful, the synthetic data it generates must faithfully mirror
real biology across multiple dimensions—from basic gene expression
distributions to complex cell-cell interactions, batch variations, and
developmental trajectories.

**`scSimEval`** provides a rigorous, standardized, reference-based
benchmarking framework to evaluate simulation quality. Rather than
relying on a single overall score or a handful of summary statistics, it
assesses synthetic datasets across **56 curated evaluation measures
organized into 8 foundational categories**.

Every measure is calculated by directly comparing an in silico simulated
dataset against an empirical biological reference dataset across
single-cell RNA sequencing (**scRNA-seq**), single-cell ATAC sequencing
(**scATAC-seq**), or joint **multiomics** (paired and unpaired RNA +
ATAC) modalities.

------------------------------------------------------------------------

## Summary of the 8 Canonical Categories

| Category | Domain | Measures | Supported Modalities | Primary Evaluation Focus |
|:--:|:---|:--:|:---|:---|
| **I** | **Distributional Properties** | **14** | scRNA-seq, scATAC-seq, Multiomics | 10 univariate (1D) and 4 bivariate (2D) statistical distances comparing fundamental cell-level, feature-level, and joint property distributions. |
| **II** | **Sparsity & Zero-Inflation Dynamics** | **5** | scRNA-seq, scATAC-seq, Multiomics | Evaluates whether dropout rates, zero fractions, and non-linear sparsity curves match real experimental technologies. |
| **III** | **Cellular Manifold & Clustering Concordance** | **8** | scRNA-seq, scATAC-seq, Multiomics | Checks whether simulated cells preserve high-dimensional geometry, cell-type clustering boundaries, and avoid artificial cell states. |
| **IV** | **Technical Variation & Batch Mixing Dynamics** | **5** | scRNA-seq, scATAC-seq, Multiomics | Quantifies how realistically technical batch effects and donor-to-donor variations intermix across local neighborhoods without obscuring biological signal. |
| **V** | **Biological Signal & Differential Fidelity** | **12** | scRNA-seq (6 DE) & scATAC-seq (6 DA) | Tests whether real marker genes (DEGs) and accessible peaks (DARs) can be reliably recovered via pseudobulk testing without false-positive inflation. |
| **VI** | **Multiomics Integration & Alignment** | **5** | Paired & Unpaired Multiomics | Evaluates cross-modal alignment in shared latent space, assessing whether RNA and ATAC layers coherently integrate while preserving cell identities. |
| **VII** | **Trajectory & Lineage Topology** | **5** | scRNA-seq | Evaluates continuous developmental paths, branching tree topologies, cell pseudotime progression, and dynamic driver genes using the dynverse framework. |
| **VIII** | **Computational Scalability & Efficiency** | **2** | All Modalities | Hardware-normalized assessment of compute runtime (Core-Seconds) and resident memory usage (Peak RAM per Allocated Core). |
| **Total** |  | **56** | **Full Benchmark Suite** | **Comprehensive, multi-faceted validation of single-cell and multiomics simulators.** |

------------------------------------------------------------------------

## Complete Master Reference Table (All 56 Measures)

The table below lists all 56 evaluation measures implemented in
`scSimEval`. For each metric, it highlights the internal metric
identifier, category, the ideal scoring direction, and what the metric
captures.

| \# | Metric ID | Display Name | Category | Best Direction | Plain-English Evaluation Focus |
|:--:|:---|:---|:---|:--:|:---|
| 1 | `Wasserstein` | 1D Wasserstein Distance | I. Distributional | Lower is better | Physical transport work needed to transform simulated distributions into reference distributions across 14 univariate properties. |
| 2 | `KS` | Kolmogorov–Smirnov Distance | I. Distributional | Lower is better | Maximum vertical divergence between the cumulative distribution curves of real and simulated properties. |
| 3 | `MAD` | Median Absolute Deviation | I. Distributional | Lower is better | Median absolute discrepancy across aligned quantile points; highly resistant to extreme outliers. |
| 4 | `MAE` | Mean Absolute Error | I. Distributional | Lower is better | Average absolute discrepancy across the entire range of cell and feature properties. |
| 5 | `RMSE` | Root Mean Squared Error | I. Distributional | Lower is better | Average squared discrepancy across quantiles, placing greater penalty on severe distribution tail mismatches. |
| 6 | `Bhattacharyya` | Bhattacharyya Distance | I. Distributional | Lower is better | Measures the geometric overlap and statistical divergence between reference and simulated probability histograms. |
| 7 | `KDE_ISE` | Kernel Density Error | I. Distributional | Lower is better | Integrated squared difference between smoothed continuous probability density curves. |
| 8 | `Runs_Statistic` | Wald–Wolfowitz Runs | I. Distributional | Closer to 0 is better | Tests whether real and simulated values interleave randomly without clustering by dataset source. |
| 9 | `NN_Mismatch` | 1D Nearest-Neighbor Mismatch | I. Distributional | Lower is better | Fraction of pooled data points whose nearest neighbor comes from the wrong dataset. |
| 10 | `Between_Dataset_Silh_Global` | 1D Inter-Dataset Silhouette | I. Distributional | Closer to 0 is better | Verifies that reference and simulated data do not artificially form separate clusters along a single property axis. |
| 11 | `Fasano_Franceschini_2D_KS` | 2D Kolmogorov–Smirnov | I. Distributional | Lower is better | Two-dimensional extension of the KS test measuring divergence across joint bivariate quadrant distributions. |
| 12 | `EMD_2D` | 2D Earth Mover’s Distance | I. Distributional | Lower is better | Optimal transport work required to morph 2D bivariate property density surfaces into one another. |
| 13 | `NN_Mismatch_2D` | 2D Nearest-Neighbor Mismatch | I. Distributional | Lower is better | Evaluates local joint neighborhood mixing across bivariate property pairs. |
| 14 | `Between_Dataset_Silh_2D_Global` | 2D Inter-Dataset Silhouette | I. Distributional | Closer to 0 is better | Ensures real and simulated cells do not artificially segregate into distinct clusters in 2D property space. |
| 15 | `midpoint_error` | Dropout Midpoint Error | II\. Sparsity Dynamics | Lower is better | Absolute shift in the expression level at which genes have a 50% probability of dropout (scRNA-seq). |
| 16 | `slope_error` | Dropout Slope Error | II\. Sparsity Dynamics | Lower is better | Discrepancy in the rate of transition from silent to detected as expression increases (scRNA-seq). |
| 17 | `r2_discrepancy` | Sparsity Fit Discrepancy | II\. Sparsity Dynamics | Lower is better | Difference in how tightly feature sparsity follows empirical trend lines compared to real biological data. |
| 18 | `curve_rmse` | Accessibility Sparsity RMSE | II\. Sparsity Dynamics | Lower is better | Root mean squared error between peak accessibility and non-zero proportion curves (scATAC-seq). |
| 19 | `excess_zero_weight_discrepancy` | Excess Zero Discrepancy | II\. Sparsity Dynamics | Lower is better | Quantifies differences in the proportion of zeroes attributable to technical dropout beyond biological expectation. |
| 20 | `MMD` | Maximum Mean Discrepancy | III\. Cellular Manifold | Lower is better | Multi-scale kernel distance in high-dimensional PCA space, evaluating overall multivariate point-cloud alignment. |
| 21 | `Generative_Precision` | Generative Precision | III\. Cellular Manifold | Higher is better | Proportion of simulated cells that fall within the valid support of real cell states (penalizes hallucinated cells). |
| 22 | `Generative_Recall` | Generative Recall | III\. Cellular Manifold | Higher is better | Proportion of real biological cell diversity successfully captured by the simulator (penalizes mode collapse). |
| 23 | `silhouette_sim` | Cluster Silhouette (ASW) | III\. Cellular Manifold | Higher is better | Measures how tightly cells of the same cell type cluster together compared to neighboring cell types. |
| 24 | `cdi` | Clustering Deviation Index | III\. Cellular Manifold | Lower is better | Penalized statistical likelihood measuring whether simulated cell populations are over-clustered or under-clustered. |
| 25 | `neighborhood_purity` | Neighborhood Purity | III\. Cellular Manifold | Higher is better | Average percentage of a cell’s nearest neighbors that belong to the exact same biological cell type. |
| 26 | `ari` | Adjusted Rand Index | III\. Cellular Manifold | Higher is better | Agreement between unsupervised cell clusters and true biological cell types, corrected for chance. |
| 27 | `nmi` | Normalized Mutual Info | III\. Cellular Manifold | Higher is better | Information-theoretic measure of shared information between unsupervised clusters and known cell types. |
| 28 | `batch_silhouette` | Batch Silhouette | IV\. Technical Batch | Closer to 0 is better | Measures whether cells from different technical batches intermix naturally rather than separating into artificial groups. |
| 29 | `pcr_r2` | PC Regression R2 | IV\. Technical Batch | Lower is better | Percentage of total expression variance in top principal components driven purely by technical batch artifacts. |
| 30 | `cms` | Cell-Specific Mixing Score | IV\. Technical Batch | Higher is better | Statistical test evaluating whether local cell neighborhoods reflect the overall batch proportions of the experiment. |
| 31 | `isi` | Inverse Simpson Index (LISI) | IV\. Technical Batch | Higher is better | Effective number of technical batches represented in each cell’s immediate neighborhood. |
| 32 | `cross_batch_F1` | Cross-Batch Classification F1 | IV\. Technical Batch | Higher is better | Tests whether cell-type classifiers trained on one batch can accurately identify cell types in another batch. |
| 33 | `pval_dist_score` (scRNA) | DE Null P-Value Uniformity | V. Biological Signal | Higher is better | Checks whether non-differential genes produce flat, uniform p-values under the null hypothesis (prevents false positives). |
| 34 | `pval_dist_score` (scATAC) | DA Null P-Value Uniformity | V. Biological Signal | Higher is better | Checks whether non-differential chromatin peaks produce properly calibrated, uniform null p-values. |
| 35 | `accuracy` (scRNA) | DE Detection Accuracy | V. Biological Signal | Higher is better | Overall binary accuracy in distinguishing differential genes from non-differential genes. |
| 36 | `accuracy` (scATAC) | DA Detection Accuracy | V. Biological Signal | Higher is better | Overall binary accuracy in distinguishing differential chromatin peaks from background regions. |
| 37 | `precision` (scRNA) | DE Detection Precision | V. Biological Signal | Higher is better | Proportion of declared marker genes that correspond to genuine biological markers in the reference data. |
| 38 | `precision` (scATAC) | DA Detection Precision | V. Biological Signal | Higher is better | Proportion of declared accessible peaks that represent genuine regulatory differences in reference data. |
| 39 | `recall` (scRNA) | DE Detection Recall | V. Biological Signal | Higher is better | Sensitivity: what percentage of true biological marker genes are successfully detected in the simulated data. |
| 40 | `recall` (scATAC) | DA Detection Recall | V. Biological Signal | Higher is better | Sensitivity: what percentage of true regulatory chromatin peaks are successfully identified. |
| 41 | `f1_score` (scRNA) | DE Detection F1-Score | V. Biological Signal | Higher is better | Harmonic balance between marker detection precision and recall for single-cell gene expression. |
| 42 | `f1_score` (scATAC) | DA Detection F1-Score | V. Biological Signal | Higher is better | Harmonic balance between peak detection precision and recall for single-cell chromatin accessibility. |
| 43 | `auc` (scRNA) | DE Ranking ROC-AUC | V. Biological Signal | Higher is better | Threshold-independent ranking test: verifies that true marker genes receive stronger statistical significance than noise. |
| 44 | `auc` (scATAC) | DA Ranking ROC-AUC | V. Biological Signal | Higher is better | Threshold-independent ranking test: verifies that true differential chromatin peaks are ranked at the top. |
| 45 | `batch_asw` | Cross-Modal Batch ASW | VI\. Multiomics Alignment | Higher is better | Evaluates whether RNA and ATAC modalities intermix smoothly within matching cell types in joint latent space. |
| 46 | `cell_type_asw` | Joint Cell-Type ASW | VI\. Multiomics Alignment | Higher is better | Verifies that biological cell types remain distinctly separated from one another in the joint multiomics space. |
| 47 | `graph_connectivity` | Graph Connectivity | VI\. Multiomics Alignment | Higher is better | Checks whether cells of the same biological type form an unbroken, interconnected network across modalities. |
| 48 | `ari` | Integration ARI | VI\. Multiomics Alignment | Higher is better | Accuracy of unsupervised clustering on joint multiomics embeddings compared to ground-truth cell labels. |
| 49 | `nmi` | Integration NMI | VI\. Multiomics Alignment | Higher is better | Shared information between joint multiomics clusters and true biological cell identities. |
| 50 | `him` | HIM Topology Similarity | VII\. Lineage Trajectory | Higher is better | Combines local edge structure and global graph spectra to test whether trajectory branching trees match real biology. |
| 51 | `edge_flip` | Edge Flip Similarity | VII\. Lineage Trajectory | Higher is better | Measures the minimum number of graph edits needed to morph simulated differentiation paths into real paths. |
| 52 | `geodesic_correlation` | Geodesic Distance Corr. | VII\. Lineage Trajectory | Higher is better | Tests whether continuous pseudotime progression distances between cells are preserved along developmental paths. |
| 53 | `F1_branches` | Branch Assignment F1 | VII\. Lineage Trajectory | Higher is better | Assesses whether cells are assigned to the correct developmental lineages and terminal fate branches. |
| 54 | `wcor_features` | Driver Gene Importance Corr. | VII\. Lineage Trajectory | Higher is better | Verifies that the exact same dynamic regulatory genes govern trajectory progression in simulated and real cells. |
| 55 | `core_seconds` | Cumulative Core-Seconds | VIII\. Scalability | Lower is better | Total computational effort: elapsed wall-clock seconds multiplied by the number of allocated processor cores. |
| 56 | `memory_per_core_mib` | Peak RAM per Core (MiB) | VIII\. Scalability | Lower is better | Hardware-normalized resident memory footprint: peak RAM consumption divided by allocated cores. |

------------------------------------------------------------------------

## Session Information

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`\
`#> R version 4.6.2 (2026-10-29)`\
`#> Platform: x86_64-pc-linux-gnu`\
`#> Running under: Ubuntu 24.04.5 LTS`\
`#> `\
`#> Matrix products: default`\
`#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 `\
`#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0`\
`#> `\
`#> locale:`\
`#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              `\
`#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    `\
`#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   `\
`#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 `\
`#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            `\
`#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       `\
`#> `\
`#> time zone: Asia/Kolkata`\
`#> tzcode source: system (glibc)`\
`#> `\
`#> attached base packages:`\
`#> [1] stats     graphics  grDevices utils     datasets  methods   base     `\
`#> `\
`#> loaded via a namespace (and not attached):`\
`#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    `\
`#>  [5] xfun_0.62         cachem_1.1.0      knitr_1.52        htmltools_0.5.9  `\
`#>  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      `\
`#> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.3`\
`#> [17] compiler_4.6.2    tools_4.6.2       ragg_1.5.2        bslib_0.12.0     `\
`#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.3.0        jsonlite_2.0.0   `\
`#> [25] rlang_1.3.0       fs_2.1.0          htmlwidgets_1.6.4`
