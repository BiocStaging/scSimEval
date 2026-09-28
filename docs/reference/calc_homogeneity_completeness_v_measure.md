# Calculate Homogeneity, Completeness, and V-Measure

Integrated from scCluBench (Xu et al., AAAI 2026).

## Usage

``` r
calc_homogeneity_completeness_v_measure(pred, truth, beta = 1)
```

## Arguments

- pred:

  Vector of predicted cluster labels.

- truth:

  Vector of ground truth cell type labels.

- beta:

  Weight of completeness vs. homogeneity (default 1).

## Value

A named list containing homogeneity, completeness, and v_measure.

## Examples

``` r
pred <- factor(rep(c("A", "B"), each = 20))
truth <- factor(rep(c("A", "B"), each = 20))
calc_homogeneity_completeness_v_measure(pred, truth)
#> $homogeneity
#> [1] 1
#> 
#> $completeness
#> [1] 1
#> 
#> $v_measure
#> [1] 1
#> 
```
