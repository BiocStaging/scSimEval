# Package index

## Master Benchmarking Workflows

Core end-to-end evaluation pipelines executing the 56 curated evaluation
metrics across unimodal (scRNA-seq, scATAC-seq) and multiomics (paired,
unpaired, mosaic) simulations.

- [`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md)
  : Evaluate Comprehensive Simulation Accuracy (Unimodal
  SingleCellExperiment)
- [`evaluate_multiomics_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_accuracy.md)
  : Master Multiomics Simulation Benchmarking Pipeline
- [`evaluate_multiple_datasets()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiple_datasets.md)
  : Evaluate Multiple Single-Cell Datasets Simultaneously
- [`benchmark_simulation_suite()`](https://kabilanbio.github.io/scSimEval/reference/benchmark_simulation_suite.md)
  : Benchmark Simulation Methods Across Multiple Single-Cell Datasets
- [`combine_benchmark_results()`](https://kabilanbio.github.io/scSimEval/reference/combine_benchmark_results.md)
  [`append_benchmark_results()`](https://kabilanbio.github.io/scSimEval/reference/combine_benchmark_results.md)
  : Combine and Append Multiple Single-Cell Benchmark Results

## Scoring, Ranking & Cell Embeddings

Direction-aware score normalization, overall simulator leaderboards, and
low-dimensional cell embedding projections (UMAP, t-SNE, PCA) with
quantitative fidelity metrics.

- [`compute_method_leaderboard()`](https://kabilanbio.github.io/scSimEval/reference/compute_method_leaderboard.md)
  : Compute Simulator Method Performance Leaderboard
- [`normalize_bubble_scores()`](https://kabilanbio.github.io/scSimEval/reference/normalize_bubble_scores.md)
  : Normalize Bubble Matrix Scores Across Benchmark Datasets and
  Simulators
- [`compute_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/compute_dataset_embeddings.md)
  : Compute Low-Dimensional Embeddings for Reference and Simulated
  Single-Cell Datasets
- [`compute_embedding_quality_metrics()`](https://kabilanbio.github.io/scSimEval/reference/compute_embedding_quality_metrics.md)
  : Compute Quantitative Quality Metrics for Low-Dimensional Cell
  Embeddings
- [`plot_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/plot_dataset_embeddings.md)
  : Plot High-Dimensional Dataset Embeddings for Reference and Simulated
  Single-Cell Data

## Publication Visualizations

High-resolution comparative figures for multi-metric, multi-simulator,
distribution QC, ordination, and computational scalability benchmarking.

- [`plot_benchmark_bubble_matrix()`](https://kabilanbio.github.io/scSimEval/reference/plot_benchmark_bubble_matrix.md)
  [`plot_bubble_matrix()`](https://kabilanbio.github.io/scSimEval/reference/plot_benchmark_bubble_matrix.md)
  : Plot Multi-Dimensional Benchmarking Bubble Matrix
- [`plot_evaluation_summary()`](https://kabilanbio.github.io/scSimEval/reference/plot_evaluation_summary.md)
  [`plot_benchmark_summary_bars()`](https://kabilanbio.github.io/scSimEval/reference/plot_evaluation_summary.md)
  [`plot_summary_bars()`](https://kabilanbio.github.io/scSimEval/reference/plot_evaluation_summary.md)
  : Plot Single-Cell Simulator Evaluation Summary
- [`plot_distribution_qc()`](https://kabilanbio.github.io/scSimEval/reference/plot_distribution_qc.md)
  : Plot Single-Cell Summary Distribution Quality (Comparative Overlays)
- [`plot_metric_heatmap()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_heatmap.md)
  : Plot Multi-Simulator Comparative Metric Heatmap Across Canonical
  Categories
- [`plot_metric_boxplots()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_boxplots.md)
  : Plot Benchmark Metric Score Distributions Across Simulators
- [`plot_category_metric_bars()`](https://kabilanbio.github.io/scSimEval/reference/plot_category_metric_bars.md)
  : Plot Benchmark Metric Performance Barplots by Category
- [`plot_individual_metric_bar()`](https://kabilanbio.github.io/scSimEval/reference/plot_individual_metric_bar.md)
  : Plot Individual Metric Performance Barplot
- [`plot_metric_pca()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_pca.md)
  : Principal Component Analysis (PCA) Dashboard of Benchmark Metrics
  and Methods
- [`plot_metric_mds()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_mds.md)
  : Multi-Dimensional Scaling (MDS) Ordination of Evaluation Metrics or
  Simulators
- [`plot_scalability_benchmark()`](https://kabilanbio.github.io/scSimEval/reference/plot_scalability_benchmark.md)
  : Plot Multi-Dimensional Computational Scalability Benchmark

## Category-Specific Evaluation Suites

Focused entry points to run individual evaluation categories directly
(for metric formulas and conceptual definitions across all 56 measures,
see the Metric Catalog tab).

- [`evaluate_zero_probability_curve()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_zero_probability_curve.md)
  : Evaluate Zero-Probability Dropout Curve Fidelity (Splatter &
  ZINB-WaVE)
- [`evaluate_accessibility_sparsity_curve()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_accessibility_sparsity_curve.md)
  : Evaluate Accessibility-Sparsity Curve Concordance (simATAC)
- [`evaluate_clustering_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_clustering_metrics.md)
  : Full Clustering Performance Evaluation
- [`evaluate_batch_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_batch_metrics.md)
  : Calculate Batch Integration Metrics
- [`evaluate_deg_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_deg_fidelity.md)
  [`evaluate_rna_de_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_deg_fidelity.md)
  : Evaluate scRNA-seq Differential Expression (DE) Fidelity
  (Sub-category V-A)
- [`evaluate_atac_da_fidelity()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_atac_da_fidelity.md)
  : Evaluate scATAC-seq Differential Accessibility (DA) Fidelity
  (Sub-category V-B)
- [`evaluate_multiomics_integration()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_integration.md)
  : Evaluate Multiomics Integration & Cross-Modal Alignment Fidelity
  (Category VI)
- [`evaluate_trajectory_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_trajectory_metrics.md)
  : Full Category VII Trajectory & Lineage Topology Evaluation (dynverse
  Workflow)
- [`benchmark_resource_usage()`](https://kabilanbio.github.io/scSimEval/reference/benchmark_resource_usage.md)
  : Benchmark Execution Runtime and Memory Usage
- [`compute_core_seconds()`](https://kabilanbio.github.io/scSimEval/reference/compute_core_seconds.md)
  : Compute Core-Seconds Runtime Metric
- [`compute_memory_per_core()`](https://kabilanbio.github.io/scSimEval/reference/compute_memory_per_core.md)
  : Compute Peak RAM per Allocated CPU Core

## Bundled Example Datasets

Pre-packaged synthetic and empirical single-cell and multiomics
SingleCellExperiment benchmark datasets for immediate experimentation.

- [`example_scrna`](https://kabilanbio.github.io/scSimEval/reference/example_scrna.md)
  : Example scRNA-seq Benchmarking Dataset (Rice Bud)
- [`example_scatac`](https://kabilanbio.github.io/scSimEval/reference/example_scatac.md)
  : Example scATAC-seq Benchmarking Dataset (Rice Bud)
- [`example_multiomics`](https://kabilanbio.github.io/scSimEval/reference/example_multiomics.md)
  : Example Paired Multiomics Benchmarking Dataset (Rice Bud)
- [`example_rice_bud`](https://kabilanbio.github.io/scSimEval/reference/example_rice_bud.md)
  : Real Paired Single-Cell Multiomics Benchmark Dataset (Rice Bud,
  SingleCellExperiment)
- [`example_soybean_root`](https://kabilanbio.github.io/scSimEval/reference/example_soybean_root.md)
  : Real Unpaired Single-Cell Multiomics Benchmark Dataset (Soybean
  Root, SingleCellExperiment)
