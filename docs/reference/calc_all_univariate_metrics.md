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

Named list of univariate distance metrics.

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_all_univariate_metrics(ref, sim)
#> $MAD
#> [1] 0.2246476
#> 
#> $KS
#> [1] 0.2
#> 
#> $MAE
#> [1] 0.2742963
#> 
#> $RMSE
#> [1] 0.3427931
#> 
#> $OV
#> [1] 0.8053545
#> 
#> $Bhattacharyya
#> [1] 0.005302867
#> 
#> $Wasserstein
#> [1] 0.2742963
#> 
#> $ECDF_DiffArea
#> [1] 0.05858361
#> 
#> $Runs_Statistic
#> [1] -1.005089
#> 
#> $Runs_PValue
#> [1] 0.157427
#> 
#> $NN_Mismatch
#> [1] 0.06
#> 
#> $Between_Dataset_Silh_Global
#> [1] 0.007346281
#> 
#> $Between_Dataset_Silh_Local
#> [1] -0.0212288
#> 
```
