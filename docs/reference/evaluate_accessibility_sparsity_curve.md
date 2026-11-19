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

A named list quantifying concordance between reference and simulated
accessibility-sparsity curves:

- `curve_rmse`:

  Root mean squared error between predicted NZP curves over a 100-point
  grid spanning the reference mean range.

- `r2_discrepancy`:

  Absolute difference between simulated and reference polynomial \\R^2\\
  values.

- `ref_r_squared`:

  Polynomial \\R^2\\ of the reference dataset.

- `sim_r_squared`:

  Polynomial \\R^2\\ of the simulated dataset.

## Examples

``` r
data(example_scatac, package = "scSimEval")
evaluate_accessibility_sparsity_curve(example_scatac$ref, example_scatac$sim)
#> $curve_rmse
#> [1] 0.03827511
#> 
#> $r2_discrepancy
#> [1] 0.2758995
#> 
#> $ref_r_squared
#> [1] 0.9101528
#> 
#> $sim_r_squared
#> [1] 0.6342533
#> 
```
