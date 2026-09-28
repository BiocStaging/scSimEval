# Evaluate Comprehensive Simulation Accuracy (Unimodal Omics Layer)

Computes distribution distances, cellular and feature-level properties,
zero-inflation curves, and manifold distances between reference and
simulated datasets.

## Usage

``` r
evaluate_simulation_accuracy(
  ref_data,
  sim_data,
  compute_bivariate = TRUE,
  threads = 1,
  memory_mb = NULL,
  elapsed_time = NULL,
  peak_memory_mb = NULL,
  verbose = TRUE
)
```

## Arguments

- ref_data:

  Reference count matrix (features x cells).

- sim_data:

  Simulated count matrix (features x cells).

- compute_bivariate:

  Logical, whether to compute 2D bivariate tests. Default is TRUE.

- threads:

  Number of threads for parallel computation. Default is 1.

- memory_mb:

  Optional numeric value of peak memory allocation in MB (or
  peak_memory_mb).

- elapsed_time:

  Optional numeric value of wall-clock elapsed time in seconds.

- peak_memory_mb:

  Optional alias for `memory_mb`.

- verbose:

  Logical, whether to print execution messages. Default is TRUE.

## Value

A list containing tidy summary tables and detailed metric lists.

## Details

![scSimEval Unified Benchmarking Workflow
Diagram](figures/workflow_diagram.png)

## Examples

``` r
data(example_scrna, package = "scSimEval")
res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#> [1/5] Extracting cell-level properties...
#> [2/5] Extracting feature-level properties...
#> [3/5] Computing univariate accuracy metrics...
#> [4/5] Computing bivariate metrics (Fasano-Franceschini, Peacock, KDE zstat, 2D EMD)...
#> [5/5] Computing zero-probability and manifold distances...
#> Unimodal accuracy evaluation complete.
```
