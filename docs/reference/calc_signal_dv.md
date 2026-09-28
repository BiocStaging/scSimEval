# Detect Differential Variability (DV) Genes (Bartlett Test)

Detect Differential Variability (DV) Genes (Bartlett Test)

## Usage

``` r
calc_signal_dv(exprs_mat, cell_types)
```

## Arguments

- exprs_mat:

  Log-normalized matrix (genes x cells).

- cell_types:

  Factor of 2 cell types.

## Value

Vector of adjusted p-values.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_signal_dv(example_scrna$ref, example_scrna$cell_types)
#>      Gene_01      Gene_02      Gene_03      Gene_04      Gene_05      Gene_06 
#> 7.657754e-02 2.758638e-02 1.768190e-05 2.258438e-01 2.010465e-01 8.016458e-01 
#>      Gene_07      Gene_08      Gene_09      Gene_10      Gene_11      Gene_12 
#> 4.553091e-01 9.328136e-01 7.171864e-01 7.657754e-02 2.574963e-02 2.709945e-02 
#>      Gene_13      Gene_14      Gene_15      Gene_16      Gene_17      Gene_18 
#> 1.360342e-07 3.912340e-02 7.038018e-03 9.534963e-01 9.534963e-01 9.534963e-01 
#>      Gene_19      Gene_20      Gene_21      Gene_22      Gene_23      Gene_24 
#> 9.534963e-01 2.937155e-02 7.523990e-02 2.292232e-01 7.657754e-02 3.644729e-01 
#>      Gene_25      Gene_26      Gene_27      Gene_28      Gene_29      Gene_30 
#> 9.359562e-01 8.240076e-01 7.171864e-01 2.665657e-01 8.240076e-01 8.240076e-01 
#>      Gene_31      Gene_32      Gene_33      Gene_34      Gene_35      Gene_36 
#> 3.644729e-01 7.657754e-02 8.638967e-01 3.851983e-02 8.815885e-03 4.759076e-01 
#>      Gene_37      Gene_38      Gene_39      Gene_40      Gene_41      Gene_42 
#> 4.344485e-01 8.815885e-03 2.980833e-02 8.240076e-01 9.909542e-01 1.636242e-01 
#>      Gene_43      Gene_44      Gene_45      Gene_46      Gene_47      Gene_48 
#> 4.759076e-01 9.433325e-01 1.636242e-01 3.851983e-02 8.600326e-01 8.240076e-01 
#>      Gene_49      Gene_50      Gene_51      Gene_52      Gene_53      Gene_54 
#> 8.240076e-01 4.238035e-01 6.827933e-01 6.028167e-01 4.344485e-01 6.827933e-01 
#>      Gene_55      Gene_56      Gene_57      Gene_58      Gene_59      Gene_60 
#> 1.326235e-03 7.858222e-01 3.905269e-01 4.759076e-01 1.308364e-02 2.258438e-01 
```
