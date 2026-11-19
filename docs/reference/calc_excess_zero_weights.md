# Calculate Posterior Excess Zero Weights (zingeR & ZINB-WaVE)

Calculates the posterior probability of zero counts being technical
dropouts (excess zeros) versus biological sampling zeros under a
Negative Binomial model.

## Usage

``` r
calc_excess_zero_weights(counts)
```

## Arguments

- counts:

  Count matrix (features x cells).

## Value

A named list summarizing posterior excess-zero probabilities under a
Negative Binomial model:

- `mean_excess_zero_weight`:

  Mean posterior excess-zero weight across features.

- `median_excess_zero_weight`:

  Median posterior excess-zero weight across features.

- `mean_zero_inflation`:

  Mean estimated excess-zero fraction (`pi`) beyond Negative Binomial
  expectations.

- `gene_excess_weights`:

  Numeric vector of per-feature posterior excess-zero weights.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_excess_zero_weights(example_scrna$ref)
#> $mean_excess_zero_weight
#> [1] 0.05275921
#> 
#> $median_excess_zero_weight
#> [1] 0
#> 
#> $mean_zero_inflation
#> [1] 0.05051367
#> 
#> $gene_excess_weights
#>   [1] 0.008837151 0.000000000 0.000000000 0.061916799 0.000000000 0.000000000
#>   [7] 0.000000000 0.000000000 0.091204951 0.000000000 0.156452306 0.000000000
#>  [13] 0.000000000 0.000000000 0.000000000 0.022312528 0.000000000 0.000000000
#>  [19] 0.168793802 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#>  [25] 0.054780653 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#>  [31] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.003288746
#>  [37] 0.100556844 0.276404977 0.047334945 0.185797187 0.032953446 0.000000000
#>  [43] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#>  [49] 0.418098470 0.000000000 0.000000000 0.000000000 0.221434696 0.013160768
#>  [55] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#>  [61] 0.127330401 0.000000000 0.000000000 0.000000000 0.481131063 0.000000000
#>  [67] 0.000000000 0.276404977 0.043869620 0.003711566 0.000000000 0.000000000
#>  [73] 0.129820499 0.000000000 0.000000000 0.276404977 0.000000000 0.211976240
#>  [79] 0.165902841 0.013202566 0.053966737 0.000000000 0.070846551 0.084939589
#>  [85] 0.211976240 0.000000000 0.000000000 0.000000000 0.000000000 0.128425169
#>  [91] 0.194007676 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#>  [97] 0.126253481 0.278237588 0.126253481 0.407931896
#> 
```
