# Inverse Simpson Index (ISI) for Batch Mixing

Evaluates the diversity of batch labels within each cell's k-nearest
neighborhood using the Inverse Simpson Index (1 / sum(p_b^2)).
Implements a pure-R, lightweight formulation of LISI with optional
distance weighting.

## Usage

``` r
calc_isi(coords, batch_info, k = 30, weighted = TRUE)
```

## Arguments

- coords:

  Matrix of cell coordinates (cells x dimensions).

- batch_info:

  Factor or character vector of batch assignments.

- k:

  Number of nearest neighbors. Default is 30.

- weighted:

  Logical, whether to weight neighbor contributions by inverse distance.
  Default is TRUE.

## Value

A named list:

- isi_scores:

  Numeric vector of cell-level ISI values

- mean_isi:

  Mean ISI score (ranges from 1 to number of batches; higher indicates
  better mixing)

- median_isi:

  Median ISI score

## Examples

``` r
coords <- matrix(stats::rnorm(100), 50, 2)
batch <- factor(rep(c("B1", "B2"), length.out = 50))
calc_isi(coords, batch)
#> $isi_scores
#>  [1] 1.926602 1.936097 1.996715 1.978073 1.926113 1.969539 1.938800 1.947332
#>  [9] 1.931143 1.941027 1.919895 1.889481 1.969783 1.993213 1.991806 1.980674
#> [17] 1.756755 1.922820 1.954500 1.958572 1.999999 1.982085 1.988130 1.887140
#> [25] 1.919450 1.945164 1.744224 1.899145 1.858070 1.997629 1.960305 1.940331
#> [33] 1.998231 1.974506 1.749880 1.987254 1.765746 1.979309 1.983751 1.913723
#> [41] 1.994643 1.965708 1.933262 1.991392 1.977415 1.955302 1.784964 1.994640
#> [49] 1.990485 1.895894
#> 
#> $mean_isi
#> [1] 1.935734
#> 
#> $median_isi
#> [1] 1.954901
#> 
```
