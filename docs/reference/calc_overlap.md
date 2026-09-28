# Calculate Distribution Overlapping Index (OV)

Calculate Distribution Overlapping Index (OV)

## Usage

``` r
calc_overlap(ref, sim)
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

## Value

Area of overlap under probability densities (0 to 1).

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_overlap(ref, sim)
#> [1] 0.8805348
```
