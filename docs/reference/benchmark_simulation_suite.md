# Benchmark Simulation Methods Across Multiple Single-Cell Datasets

Unified benchmarking runner for multi-dataset and multi-simulator
studies. Manages sparse matrix coercion, `SingleCellExperiment`
conversion, parallel metric computation, and produces both granular and
cross-dataset aggregated summary tables ready for visualization.

## Usage

``` r
benchmark_simulation_suite(
  datasets,
  simulators,
  compute_bivariate = FALSE,
  threads = 1L,
  auto_sparse = TRUE,
  run_gc = TRUE,
  save_path = NULL,
  verbose = TRUE
)
```

## Arguments

- datasets:

  Named list of empirical reference datasets (`SingleCellExperiment`,
  count matrix, file path, or list containing `ref`, `cell_types`,
  `batch_info`).

- simulators:

  Named list of simulator outputs corresponding to `datasets`.

- compute_bivariate:

  Logical; compute 2D bivariate tests. Default `FALSE`.

- threads:

  Integer number of parallel workers. Default `1L`.

- auto_sparse:

  Logical; coerce dense matrices to `dgCMatrix`. Default `TRUE`.

- run_gc:

  Logical; run [`gc()`](https://rdrr.io/r/base/gc.html) between
  datasets. Default `TRUE`.

- save_path:

  Optional `.rds` path to save the output bundle.

- verbose:

  Logical; print progress logging. Default `TRUE`.

## Value

An S3 object of class
`c("scSimEval_multidataset_benchmark", "scSimEval_consolidated", "list")`
containing:

- granular_benchmark_table:

  Tidy data.frame across all dataset-simulator evaluations.

- benchmark_summary_table:

  Cross-dataset averaged fidelity scores formatted for
  [`plot_benchmark_bubble_matrix`](https://kabilanbio.github.io/scSimEval/reference/plot_benchmark_bubble_matrix.md).

- cross_dataset_stats:

  Mean, SD, SE, Min, Max, and `N_Datasets` for each simulator and
  metric.

- dataset_summary:

  Summary table of reference dimensions, sparsity, and cell-type counts.

- dataset_summary_df:

  Alias of `dataset_summary`.

- methods:

  Character vector of evaluated simulator names.

- datasets:

  Character vector of evaluated dataset names.

- is_multidataset:

  Logical `TRUE`.

- n_datasets:

  Integer count of datasets.

- n_simulators:

  Integer count of simulators.

## Examples

``` r
data(example_scrna, package = "scSimEval")
refs <- list("CohortA" = example_scrna$ref)
sims <- list("Splatter" = list("CohortA" = example_scrna$sim))
res <- benchmark_simulation_suite(refs, sims, verbose = FALSE)
print(res)
#> =====================================================================
#>   scSimEval Multi-Dataset Benchmark Results
#> =====================================================================
#> Datasets Evaluated     : 1 (CohortA)
#> Simulators Evaluated   : 1 (Splatter)
#> Total Evaluated Scores : 168 across all datasets & metrics
#> Cross-Dataset Average  : 38 averaged metric entries in summary table
#> =====================================================================
head(res$benchmark_summary_table, 8)
#>     Method                      Category                      Metric Score_Norm
#> 1 Splatter (I) Distributional Properties Between_Dataset_Silh_Global     0.9872
#> 2 Splatter (I) Distributional Properties               Bhattacharyya     0.9907
#> 3 Splatter (I) Distributional Properties                     KDE_ISE     0.8996
#> 4 Splatter (I) Distributional Properties                          KS     0.8768
#> 5 Splatter (I) Distributional Properties                         MAD     0.6454
#> 6 Splatter (I) Distributional Properties                         MAE     0.6194
#> 7 Splatter (I) Distributional Properties                 NN_Mismatch     0.7832
#> 8 Splatter (I) Distributional Properties                        RMSE     0.5627
#>   Score_Raw  Score Simulator
#> 1    0.0069 0.9872  Splatter
#> 2    0.0093 0.9907  Splatter
#> 3    0.1004 0.8996  Splatter
#> 4    0.1232 0.8768  Splatter
#> 5    1.8449 0.6454  Splatter
#> 6    3.7466 0.6194  Splatter
#> 7    0.2168 0.7832  Splatter
#> 8    7.2937 0.5627  Splatter
```
