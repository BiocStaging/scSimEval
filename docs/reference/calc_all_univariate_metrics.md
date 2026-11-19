# Compute All Univariate Distance & Accuracy Metrics for a Feature

Compute All Univariate Distance & Accuracy Metrics for a Feature

## Usage

``` r
calc_all_univariate_metrics(ref, sim, metric_prefix = "")
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

- metric_prefix:

  Optional prefix string for metric names.

## Value

A named list of 11 numeric univariate distance and similarity metrics
(names prefixed with `metric_prefix_` when provided):

- MAD:

  Median absolute difference between quantile-aligned vectors.

- KS:

  1D Kolmogorov-Smirnov maximum vertical distance between empirical
  CDFs.

- MAE:

  Mean absolute error between quantile-aligned vectors.

- RMSE:

  Root mean squared error between quantile-aligned vectors.

- Bhattacharyya:

  Bhattacharyya divergence between normalized probability vectors.

- Wasserstein:

  1D Wasserstein-1 distance between distributions.

- Runs_Statistic:

  Wald-Wolfowitz runs test z-statistic measuring distribution
  interleaving.

- Runs_PValue:

  Left-sided p-value from the Wald-Wolfowitz runs test.

- NN_Mismatch:

  Fraction of local k-NN neighborhoods showing dataset composition
  mismatch.

- Between_Dataset_Silh_Global:

  Global between-dataset silhouette width.

- KDE_ISE:

  Integrated squared error between kernel density estimates.

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_all_univariate_metrics(ref, sim)
#> $MAD
#> [1] 0.06811227
#> 
#> $KS
#> [1] 0.08
#> 
#> $MAE
#> [1] 0.1339628
#> 
#> $RMSE
#> [1] 0.2398392
#> 
#> $Bhattacharyya
#> [1] 0.00184261
#> 
#> $Wasserstein
#> [1] 0.1339628
#> 
#> $Runs_Statistic
#> [1] 1.005089
#> 
#> $Runs_PValue
#> [1] 0.842573
#> 
#> $NN_Mismatch
#> [1] 0
#> 
#> $Between_Dataset_Silh_Global
#> [1] -0.01544037
#> 
#> $KDE_ISE
#> [1] 0.008094561
#> 
```
