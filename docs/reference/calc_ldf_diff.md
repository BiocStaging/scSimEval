# Local Density Differences (ldfDiff)

Quantifies cell-specific changes in the Local Density Factor (LDF)
before and after data integration or between reference and simulated
datasets.

## Usage

``` r
calc_ldf_diff(coords_pre, coords_post, k = 15, h = 1, c = 1)
```

## Arguments

- coords_pre:

  Matrix of coordinates before integration / reference (cells x
  dimensions).

- coords_post:

  Matrix of coordinates after integration / simulated (cells x
  dimensions).

- k:

  Number of nearest neighbors. Default is 15.

- h:

  Bandwidth parameter for Gaussian kernel. Default is 1.

- c:

  Scaling constant. Default is 1.

## Value

A named list:

- diff:

  Numeric vector of absolute LDF differences per cell

- mean_ldf_diff:

  Mean absolute LDF difference (lower indicates better structure
  preservation)

- median_ldf_diff:

  Median absolute LDF difference

- ldf_pre:

  LDF in pre-integration / reference space

- ldf_post:

  LDF in post-integration / simulated space

## Examples

``` r
coords <- matrix(stats::rnorm(100), 50, 2)
batch <- factor(rep(c("B1", "B2"), length.out = 50))
calc_ldf_diff(coords, coords)
#> $diff
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [39] 0 0 0 0 0 0 0 0 0 0 0 0
#> 
#> $mean_ldf_diff
#> [1] 0
#> 
#> $median_ldf_diff
#> [1] 0
#> 
#> $ldf_pre
#>  [1] 0.5234678 0.4952380 0.5647165 0.4859537 0.5171036 0.4831887 0.5016189
#>  [8] 0.7046081 0.7335257 0.5775315 0.4987396 0.6985569 0.5652276 0.5290416
#> [15] 0.5971797 0.5168758 0.4871057 0.4795465 0.5324614 0.4796246 0.5375272
#> [22] 0.4851946 0.4994395 0.8866663 0.5122622 0.5433988 0.6205420 0.4795286
#> [29] 0.4934460 0.6149252 0.4836573 0.6141218 0.7934550 0.6945098 0.6099527
#> [36] 0.6918709 0.6374834 0.6478343 0.5094234 0.4849180 0.4874206 0.7205556
#> [43] 0.5259607 0.7819178 0.5416860 0.4888105 0.4950431 0.6156297 0.5397597
#> [50] 0.5611875
#> 
#> $ldf_post
#>  [1] 0.5234678 0.4952380 0.5647165 0.4859537 0.5171036 0.4831887 0.5016189
#>  [8] 0.7046081 0.7335257 0.5775315 0.4987396 0.6985569 0.5652276 0.5290416
#> [15] 0.5971797 0.5168758 0.4871057 0.4795465 0.5324614 0.4796246 0.5375272
#> [22] 0.4851946 0.4994395 0.8866663 0.5122622 0.5433988 0.6205420 0.4795286
#> [29] 0.4934460 0.6149252 0.4836573 0.6141218 0.7934550 0.6945098 0.6099527
#> [36] 0.6918709 0.6374834 0.6478343 0.5094234 0.4849180 0.4874206 0.7205556
#> [43] 0.5259607 0.7819178 0.5416860 0.4888105 0.4950431 0.6156297 0.5397597
#> [50] 0.5611875
#> 
```
