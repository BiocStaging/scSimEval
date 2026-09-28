# Detect Differential Distribution (DD) Genes (Kolmogorov-Smirnov Test)

Detect Differential Distribution (DD) Genes (Kolmogorov-Smirnov Test)

## Usage

``` r
calc_signal_dd(exprs_mat, cell_types)
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
calc_signal_dd(example_scrna$ref, example_scrna$cell_types)
#>   Gene_01   Gene_02   Gene_03   Gene_04   Gene_05   Gene_06   Gene_07   Gene_08 
#> 0.9895096 0.5642625 0.5642625 0.8388210 0.4794990 0.9895096 0.6921413 0.9895096 
#>   Gene_09   Gene_10   Gene_11   Gene_12   Gene_13   Gene_14   Gene_15   Gene_16 
#> 0.8388210 0.6921413 0.5642625 0.9895096 0.4794990 0.2070355 0.5642625 0.9895096 
#>   Gene_17   Gene_18   Gene_19   Gene_20   Gene_21   Gene_22   Gene_23   Gene_24 
#> 0.9895096 0.9895096 0.9895096 0.9176670 0.9176670 0.2096603 0.8388210 0.8388210 
#>   Gene_25   Gene_26   Gene_27   Gene_28   Gene_29   Gene_30   Gene_31   Gene_32 
#> 0.9895096 0.9895096 0.9895096 0.9895096 0.9895096 0.4794990 0.9895096 0.9895096 
#>   Gene_33   Gene_34   Gene_35   Gene_36   Gene_37   Gene_38   Gene_39   Gene_40 
#> 0.9895096 0.9176670 0.9895096 0.9895096 0.9895096 0.9895096 0.9895096 0.9895096 
#>   Gene_41   Gene_42   Gene_43   Gene_44   Gene_45   Gene_46   Gene_47   Gene_48 
#> 0.9895096 0.9176670 0.9895096 0.9895096 0.8388210 0.9895096 0.9895096 0.9895096 
#>   Gene_49   Gene_50   Gene_51   Gene_52   Gene_53   Gene_54   Gene_55   Gene_56 
#> 0.9895096 0.9895096 0.8388210 0.6921413 0.5642625 0.9895096 0.9176670 0.9895096 
#>   Gene_57   Gene_58   Gene_59   Gene_60 
#> 0.9895096 0.2095607 0.5642625 0.9895096 
```
