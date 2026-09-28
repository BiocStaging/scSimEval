# Calculate Adjusted Rand Index (ARI)

Calculate Adjusted Rand Index (ARI)

## Usage

``` r
calc_ari(pred, truth)
```

## Arguments

- pred:

  Predicted cluster labels.

- truth:

  Ground truth cell type labels.

## Value

ARI value between -1 and 1.

## Examples

``` r
pred <- factor(rep(c("A", "B"), each = 20))
truth <- factor(rep(c("A", "B"), each = 20))
calc_ari(pred, truth)
#> [1] 1
```
