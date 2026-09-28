# Calculate Fowlkes-Mallows Index (FMI)

Integrated from scCluBench (Xu et al., AAAI 2026). Computes the
geometric mean of pairwise cluster precision and recall.

## Usage

``` r
calc_fmi(pred, truth)
```

## Arguments

- pred:

  Vector of predicted cluster labels.

- truth:

  Vector of ground truth cell type labels.

## Value

FMI score between 0 and 1.

## Examples

``` r
pred <- factor(rep(c("A", "B"), each = 20))
truth <- factor(rep(c("A", "B"), each = 20))
calc_fmi(pred, truth)
#> [1] 1
```
