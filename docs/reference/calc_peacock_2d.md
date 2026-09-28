# Calculate Peacock 2D Kolmogorov-Smirnov Test Statistic

Integrated from HelenaLC/simulation-comparison.

## Usage

``` r
calc_peacock_2d(ref_mat, sim_mat)
```

## Arguments

- ref_mat:

  2-column numeric matrix for reference.

- sim_mat:

  2-column numeric matrix for simulation.

## Value

Peacock test statistic.

## Examples

``` r
ref <- matrix(stats::rnorm(100), 50, 2)
sim <- matrix(stats::rnorm(100), 50, 2)
calc_peacock_2d(ref, sim)
#> [1] NA
```
