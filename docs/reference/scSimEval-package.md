# scSimEval: Unified Benchmarking and Accuracy Evaluation for Single-Cell Multiomics Simulation Methods

A simple and unified evaluation toolkit for benchmarking single-cell and
multiomics data simulation methods. Evaluates how realistic simulated
data are compared to real biological datasets across three multiomics
experimental designs: (i) paired multiomics (simultaneous co-assay from
identical cells, e.g. 10x Chromium Multiome, SHARE-seq, SNARE-seq),
which unlocks the full Category 6 multiomics integration and cross-modal
alignment suite (joint-embedding Batch ASW, Cell-Type ASW, Graph
Connectivity, Integration ARI, and Integration NMI); (ii) unpaired
multiomics (independent scRNA-seq and scATAC-seq profiled from separate
cells of the same tissue), supported through per-modality unimodal
evaluation (Categories 1-5, 7, 8) and joint-embedding integration
fidelity; and (iii) mosaic multiomics (partial co-measurement where only
a subset of cells has both modalities, e.g. DOGMA-seq or mosaic pooling
designs), handled by evaluating the paired and unimodal cell fractions
across all 8 canonical categories. The toolkit quantifies simulation
realism using 56 evaluation measures organized into 8 canonical
categories: (1) Distributional Properties: 10 univariate metrics
(Wasserstein, KS, MAD, MAE, RMSE, Bhattacharyya, KDE-ISE, Runs
statistic, NN mismatch, Inter-dataset silhouette) and 4 bivariate
metrics (Fasano-Franceschini 2D KS, 2D EMD, 2D NN mismatch, 2D
Inter-dataset silhouette) evaluated across 7 cell-level, 7
feature-level, and 4 bivariate property pairs, (2) Sparsity &
Zero-Inflation Dynamics: logistic dropout curve midpoint/slope/R2
(scRNA-seq), quadratic accessibility-sparsity curve RMSE/R2
(scATAC-seq), and posterior excess zero weight discrepancy, (3) Cellular
Manifold & Clustering Concordance: MMD, Generative Precision/Recall,
Silhouette ASW, Clustering Deviation Index (CDI), Neighborhood Purity,
ARI, and NMI (with auxiliary Davies-Bouldin Index), (4) Technical
Variation & Batch Mixing Dynamics: Batch Silhouette, PCR R2, CMS,
Inverse Simpson Index (ISI), and Cross-Batch F1, (5) Biological Signal &
Differential Fidelity: pseudobulk edgeR-QLF null P-value uniformity,
accuracy, precision, recall, F1-score, and ROC-AUC across parallel
scRNA-seq DE and scATAC-seq DA sub-categories, (6) Multiomics
Integration & Alignment: joint-embedding Batch ASW, Cell-Type ASW, Graph
Connectivity, Integration ARI, and Integration NMI, (7) Trajectory &
Lineage Topology: HVG-filtered dynverse benchmark (Saelens et al., 2019)
comprising HIM topology similarity, edge flip similarity, geodesic
distance correlation, branch assignment F1, and Random Forest feature
importance weighted correlation (wcor_features), (8) Computational
Scalability & Efficiency: core-normalized execution runtime
(core_seconds) and memory footprint per allocated core
(memory_per_core_mib), with automated multi-core scaling and unified
memory_mb resource benchmarking. Includes non-linear low-dimensional
cell embedding projections (UMAP, t-SNE, PCA) with quantitative cluster
separability (Silhouette) and fidelity (ARI) metrics, automated score
normalization, leaderboards, and publication visualizations. Interactive
GUI exploration is available via the companion `scSimEvalStudio`
package.

## See also

Useful links:

- <https://kabilanbio.github.io/scSimEval/>

- <https://github.com/kabilanbio/scSimEval>

- Report bugs at <https://github.com/kabilanbio/scSimEval/issues>

## Author

**Maintainer**: Kabilan Sakthivel <kabilan151414@gmail.com>
([ORCID](https://orcid.org/0009-0002-9422-0938))

Authors:

- Kabilan Sakthivel <kabilan151414@gmail.com>
  ([ORCID](https://orcid.org/0009-0002-9422-0938))

- Dr Dwijesh Chandra Mishra <dwij.mishra@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-3327-6156)) \[thesis advisor\]

- Dr Shashi Bhushan Lal <sblall16@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-2958-6178))

- Dr Ulavappa Angadi <ubaangadi.iasri@icar.org.in>
  ([ORCID](https://orcid.org/0000-0003-0012-0179))

- Dr Sudhir Srivastava <sudhir0401bm@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-4990-7693))

- Dr Neeraj Budhlakoti <neeraj35669@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-3380-7831))

- Dr Krishna Kumar Chaturvedi <kkcchaturvedi@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-6529-5466))

- Dr Ranjeet Ranjan Kumar <ranjeetranjaniari@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-2279-6340))

- Dr Navin Chandra Gupta <navinbtc@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-6502-448X))

Other contributors:

- Dr Sharanbasappa <smadival509@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-6571-9070)) \[contributor\]
