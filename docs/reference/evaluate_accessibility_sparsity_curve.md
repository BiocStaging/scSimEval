# Evaluate Accessibility-Sparsity Curve Concordance (simATAC)

Compares the non-linear relationship between peak mean accessibility and
non-zero proportion (NZP) between reference and simulated scATAC-seq
datasets.

## Usage

``` r
evaluate_accessibility_sparsity_curve(ref_data, sim_data, poly_degree = 2)
```

## Arguments

- ref_data:

  Reference count matrix (peaks x cells).

- sim_data:

  Simulated count matrix (peaks x cells).

- poly_degree:

  Degree of polynomial (default: 2).

## Value

A list of reference and simulated curve parameters, absolute
discrepancies, and curve prediction RMSE.

## Examples

``` r
data(example_scrna, package = "scSimEval")
evaluate_accessibility_sparsity_curve(example_scrna$ref, example_scrna$sim)
#> $ref_c0
#> [1] 0.5117571
#> 
#> $ref_c1
#> [1] 0.1032392
#> 
#> $ref_c2
#> [1] -0.006616585
#> 
#> $ref_r_squared
#> [1] 0.7560448
#> 
#> $sim_c0
#> [1] 0.5033413
#> 
#> $sim_c1
#> [1] 0.1126697
#> 
#> $sim_c2
#> [1] -0.008578903
#> 
#> $sim_r_squared
#> [1] 0.6265665
#> 
#> $delta_c0
#> [1] 0.008415788
#> 
#> $delta_c1
#> [1] 0.009430472
#> 
#> $delta_c2
#> [1] 0.001962318
#> 
#> $delta_r_squared
#> [1] 0.1294783
#> 
#> $delta_spearman
#> [1] 0.04562409
#> 
#> $curve_rmse
#> [1] 0.0177535
#> 
```
