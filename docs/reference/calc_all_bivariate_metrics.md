# Compute All Bivariate (2D) Distance & Accuracy Metrics for Joint Distributions

Compute All Bivariate (2D) Distance & Accuracy Metrics for Joint
Distributions

## Usage

``` r
calc_all_bivariate_metrics(ref_mat, sim_mat, metric_prefix = "", threads = 1)
```

## Arguments

- ref_mat:

  2-column matrix for reference.

- sim_mat:

  2-column matrix for simulation.

- metric_prefix:

  Optional prefix string for metric names.

- threads:

  CPU threads.

## Value

A named list of 4 numeric 2D bivariate distributional metrics (names
prefixed with `metric_prefix_` when provided):

- Fasano_Franceschini_2D_KS:

  Fasano-Franceschini 2D Kolmogorov-Smirnov test statistic.

- EMD_2D:

  Normalized 2D Earth Mover's Distance between kernel density surfaces.

- NN_Mismatch_2D:

  Fraction of 2D k-NN neighborhoods with dataset composition mismatch.

- Between_Dataset_Silh_2D_Global:

  Global between-dataset silhouette width in standardized 2D space.

## Examples

``` r
ref <- matrix(stats::rnorm(100), 50, 2)
sim <- matrix(stats::rnorm(100), 50, 2)
calc_all_bivariate_metrics(ref, sim)
#> $Fasano_Franceschini_2D_KS
#> [1] 1450
#> 
#> $EMD_2D
#> [1] 0.08436649
#> 
#> $NN_Mismatch_2D
#> [1] 0
#> 
#> $Between_Dataset_Silh_2D_Global
#> [1] 0.01767296
#> 
```
