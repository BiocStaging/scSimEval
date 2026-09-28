# Calculate Delta Variance (Pseudoreplication Bias Metric)

Measures the difference in gene expression variance between true
biological replicates and cell-shuffled pseudoreplicates (Squair et al.,
Nature Communications 2021). Positive delta variance indicates that
biological replicates have higher inter-individual variability than
randomly pooled pseudoreplicates.

## Usage

``` r
calc_delta_variance(
  counts,
  replicates,
  conditions,
  cell_types = NULL,
  n_perm = 3
)
```

## Arguments

- counts:

  Count matrix (features x cells).

- replicates:

  Factor of biological replicate / batch assignments for each cell.

- conditions:

  Factor of condition or cell-type labels.

- cell_types:

  Optional cell type vector to evaluate per cell type or overall.

- n_perm:

  Number of permutation shuffles to compute mean pseudoreplicate
  variance. Default is 3.

## Value

A named list:

- delta_variance:

  Numeric vector of gene-level delta variance values

- mean_delta_variance:

  Mean delta variance across all genes

- median_delta_variance:

  Median delta variance across all genes

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_delta_variance(example_scrna$ref,
                    replicates = example_scrna$batch_info,
                    conditions = example_scrna$cell_types)
#> $delta_variance
#>     Gene_01     Gene_02     Gene_03     Gene_04     Gene_05     Gene_06 
#> -10940591.7   2432260.2  14123334.8  -2854588.6  47604863.8   9901503.6 
#>     Gene_07     Gene_08     Gene_09     Gene_10     Gene_11     Gene_12 
#>   3733923.6   1769357.2  26565363.6 -28362940.4   -758448.1 -22377436.7 
#>     Gene_13     Gene_14     Gene_15     Gene_16     Gene_17     Gene_18 
#>   6199826.6  -2993490.4   5802911.0 -14398700.1   2892876.6   4871112.8 
#>     Gene_19     Gene_20     Gene_21     Gene_22     Gene_23     Gene_24 
#>   2326013.7 -29697716.8   4309671.9   1710500.4   8554148.6  14142670.2 
#>     Gene_25     Gene_26     Gene_27     Gene_28     Gene_29     Gene_30 
#>  35878706.6  -2952790.4 -21213909.5   5529417.1   3267786.0    927488.4 
#>     Gene_31     Gene_32     Gene_33     Gene_34     Gene_35     Gene_36 
#>   2840249.6  -3172399.6  12226394.1  -4849385.7   6961863.5    914659.6 
#>     Gene_37     Gene_38     Gene_39     Gene_40     Gene_41     Gene_42 
#>   2765727.5    847562.4    841130.7   4828235.3 -19787527.1   1919476.1 
#>     Gene_43     Gene_44     Gene_45     Gene_46     Gene_47     Gene_48 
#>   2176420.7   4904627.5   -550599.1 -32764966.7   2805171.1    944679.1 
#>     Gene_49     Gene_50     Gene_51     Gene_52     Gene_53     Gene_54 
#> -22811344.1   9479473.6    315241.3 -43707187.1   7946736.2   1483283.6 
#>     Gene_55     Gene_56     Gene_57     Gene_58     Gene_59     Gene_60 
#>   8810718.9 -19070147.8 -46345376.8    940847.2   -231595.9   5585329.3 
#> 
#> $mean_delta_variance
#> [1] -795993
#> 
#> $median_delta_variance
#> [1] 1844417
#> 
```
