# Local Structure Preservation Metric

Calculates the proportion of overlapping k-nearest neighbors between an
original/reference embedding and an integrated/simulated embedding
(adapted from Seurat LocalStruct and CellMixS locStructure).

## Usage

``` r
calc_local_structure_metric(coords_pre, coords_post, k = 30)
```

## Arguments

- coords_pre:

  Matrix of coordinates before integration / reference (cells x
  dimensions).

- coords_post:

  Matrix of coordinates after integration / simulated (cells x
  dimensions).

- k:

  Number of nearest neighbors. Default is 30.

## Value

A named list:

- cell_overlaps:

  Numeric vector of overlap fractions per cell

- mean_local_structure:

  Mean overlap fraction across cells (higher indicates better
  preservation)

- median_local_structure:

  Median overlap fraction

## Examples

``` r
coords <- matrix(stats::rnorm(100), 50, 2)
batch <- factor(rep(c("B1", "B2"), length.out = 50))
calc_local_structure_metric(coords, coords)
#> $cell_overlaps
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [39] 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> $mean_local_structure
#> [1] 1
#> 
#> $median_local_structure
#> [1] 1
#> 
```
