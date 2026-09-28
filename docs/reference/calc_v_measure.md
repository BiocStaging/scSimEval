# Calculate V-Measure

Calculate V-Measure

## Usage

``` r
calc_v_measure(pred, truth, beta = 1)
```

## Arguments

- pred:

  Vector of predicted cluster labels.

- truth:

  Vector of ground truth cell type labels.

- beta:

  Weight parameter (default 1).

## Value

V-measure score between 0 and 1.

## Examples

``` r
pred <- factor(rep(c("A", "B"), each = 20))
truth <- factor(rep(c("A", "B"), each = 20))
calc_v_measure(pred, truth)
#> [1] 1
```
