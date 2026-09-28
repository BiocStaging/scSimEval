# Detect Differential Proportion / Zero-Inflation (DP) Genes (Chisq Test)

Detect Differential Proportion / Zero-Inflation (DP) Genes (Chisq Test)

## Usage

``` r
calc_signal_dp(exprs_mat, cell_types, threshold = 0)
```

## Arguments

- exprs_mat:

  Count or expression matrix (genes x cells).

- cell_types:

  Factor of 2 cell types.

- threshold:

  Value threshold for zero/detection status (default 0).

## Value

Vector of adjusted p-values.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_signal_dp(example_scrna$ref, example_scrna$cell_types)
#> Gene_01 Gene_02 Gene_03 Gene_04 Gene_05 Gene_06 Gene_07 Gene_08 Gene_09 Gene_10 
#>       1       1       1       1       1       1       1       1       1       1 
#> Gene_11 Gene_12 Gene_13 Gene_14 Gene_15 Gene_16 Gene_17 Gene_18 Gene_19 Gene_20 
#>       1       1       1       1       1       1       1       1       1       1 
#> Gene_21 Gene_22 Gene_23 Gene_24 Gene_25 Gene_26 Gene_27 Gene_28 Gene_29 Gene_30 
#>       1       1       1       1       1       1       1       1       1       1 
#> Gene_31 Gene_32 Gene_33 Gene_34 Gene_35 Gene_36 Gene_37 Gene_38 Gene_39 Gene_40 
#>       1       1       1       1       1       1       1       1       1       1 
#> Gene_41 Gene_42 Gene_43 Gene_44 Gene_45 Gene_46 Gene_47 Gene_48 Gene_49 Gene_50 
#>       1       1       1       1       1       1       1       1       1       1 
#> Gene_51 Gene_52 Gene_53 Gene_54 Gene_55 Gene_56 Gene_57 Gene_58 Gene_59 Gene_60 
#>       1       1       1       1       1       1       1       1       1       1 
```
