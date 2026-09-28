# Calculate Normalized Mutual Information (NMI)

Calculate Normalized Mutual Information (NMI)

## Usage

``` r
calc_nmi(pred, truth)
```

## Arguments

- pred:

  Predicted cluster labels.

- truth:

  Ground truth cell type labels.

## Value

NMI value between 0 and 1.

## Examples

``` r
pred <- factor(rep(c("A", "B"), each = 20))
truth <- factor(rep(c("A", "B"), each = 20))
calc_nmi(pred, truth)
#> A 
#> 1 
```
