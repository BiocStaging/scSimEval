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

  Count matrix (genes x cells).

## Value

A list with mean excess zero weight, gene-level weights, and estimated
zero-inflation rate.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_excess_zero_weights(example_scrna$ref)
#> $mean_excess_zero_weight
#> [1] 0.125175
#> 
#> $median_excess_zero_weight
#> [1] 0.05941929
#> 
#> $mean_zero_inflation
#> [1] 0.02386441
#> 
#> $gene_excess_weights
#>    Gene_01    Gene_02    Gene_03    Gene_04    Gene_05    Gene_06    Gene_07 
#> 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.09777775 0.40423218 
#>    Gene_08    Gene_09    Gene_10    Gene_11    Gene_12    Gene_13    Gene_14 
#> 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.44714706 
#>    Gene_15    Gene_16    Gene_17    Gene_18    Gene_19    Gene_20    Gene_21 
#> 0.08562163 0.18806978 0.27267149 0.06694220 0.00000000 0.00000000 0.01533320 
#>    Gene_22    Gene_23    Gene_24    Gene_25    Gene_26    Gene_27    Gene_28 
#> 0.00000000 0.52388338 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 
#>    Gene_29    Gene_30    Gene_31    Gene_32    Gene_33    Gene_34    Gene_35 
#> 0.49862464 0.28463070 0.00000000 0.02307153 0.31073986 0.00000000 0.06367513 
#>    Gene_36    Gene_37    Gene_38    Gene_39    Gene_40    Gene_41    Gene_42 
#> 0.15088815 0.24689382 0.00000000 0.09061268 0.31942074 0.21672462 0.18998773 
#>    Gene_43    Gene_44    Gene_45    Gene_46    Gene_47    Gene_48    Gene_49 
#> 0.07329427 0.00000000 0.22786326 0.00000000 0.34310829 0.06448393 0.05516346 
#>    Gene_50    Gene_51    Gene_52    Gene_53    Gene_54    Gene_55    Gene_56 
#> 0.00000000 0.30971601 0.19499308 0.28160864 0.37442912 0.00000000 0.27711280 
#>    Gene_57    Gene_58    Gene_59    Gene_60 
#> 0.26124999 0.45666487 0.09386449 0.00000000 
#> 
```
