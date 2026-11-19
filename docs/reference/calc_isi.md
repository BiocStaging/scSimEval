# Inverse Simpson Index (ISI) for Batch Mixing

Evaluates the diversity of batch labels within each cell's k-nearest
neighborhood using the Inverse Simpson Index (1 / sum(p_b^2)).
Implements a pure-R, lightweight formulation of LISI with optional
distance weighting.

## Usage

``` r
calc_isi(
  coords,
  batch_info,
  k = 30,
  weighted = TRUE,
  threads = 1,
  knn_res = NULL
)
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

- threads:

  Number of CPU threads (default 1).

- knn_res:

  Optional precomputed k-NN list (\`idx\`, \`dists\`) from
  \`.knn_coords\`.

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
#>  [1] 1.959196 1.984330 1.997477 1.990234 1.984297 1.978011 1.963880 1.999958
#>  [9] 1.936712 1.996648 1.990049 1.980124 1.984232 1.955490 1.997168 1.991017
#> [17] 1.999136 1.977667 1.990014 1.998647 1.999498 1.997099 1.996292 1.999698
#> [25] 1.966916 1.999930 1.917644 1.989364 1.935741 1.962303 1.949662 1.993992
#> [33] 1.969402 1.999474 1.998409 1.986547 1.987219 1.997075 1.982406 1.995580
#> [41] 1.994262 1.992530 1.995701 1.965695 1.994997 1.997942 1.957628 1.994393
#> [49] 1.968758 1.999202
#> 
#> $mean_isi
#> [1] 1.982793
#> 
#> $median_isi
#> [1] 1.990141
#> 
```
