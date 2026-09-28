# Detect Bimodal Distribution (BD) Genes (Bimodal Separation Index)

Detect Bimodal Distribution (BD) Genes (Bimodal Separation Index)

## Usage

``` r
calc_signal_bd(exprs_mat, cell_types)
```

## Arguments

- exprs_mat:

  Matrix (genes x cells).

- cell_types:

  Factor of 2 cell types.

## Value

Numeric vector of Bimodality Index values.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_signal_bd(example_scrna$ref, example_scrna$cell_types)
#>    Gene_01    Gene_02    Gene_03    Gene_04    Gene_05    Gene_06    Gene_07 
#> 0.15353894 0.57529828 0.45592597 0.27766187 0.52927974 0.04833610 0.47021143 
#>    Gene_08    Gene_09    Gene_10    Gene_11    Gene_12    Gene_13    Gene_14 
#> 0.12805879 0.19109691 0.50898100 0.60771870 0.21182865 0.54809053 0.68211347 
#>    Gene_15    Gene_16    Gene_17    Gene_18    Gene_19    Gene_20    Gene_21 
#> 0.51942805 0.03626745 0.13049770 0.17832521 0.11657537 0.36575190 0.23370295 
#>    Gene_22    Gene_23    Gene_24    Gene_25    Gene_26    Gene_27    Gene_28 
#> 0.38613408 0.21391309 0.29013226 0.25822364 0.01256994 0.06623510 0.07172948 
#>    Gene_29    Gene_30    Gene_31    Gene_32    Gene_33    Gene_34    Gene_35 
#> 0.12867062 0.45456925 0.13620441 0.19546751 0.01243096 0.33224383 0.28943024 
#>    Gene_36    Gene_37    Gene_38    Gene_39    Gene_40    Gene_41    Gene_42 
#> 0.23230470 0.27366321 0.18048584 0.23140746 0.07293277 0.09747107 0.29419980 
#>    Gene_43    Gene_44    Gene_45    Gene_46    Gene_47    Gene_48    Gene_49 
#> 0.18896420 0.08713736 0.13104329 0.32947316 0.10525201 0.06234958 0.06758319 
#>    Gene_50    Gene_51    Gene_52    Gene_53    Gene_54    Gene_55    Gene_56 
#> 0.01749721 0.36271040 0.43127717 0.31431313 0.08717148 0.41562115 0.02319334 
#>    Gene_57    Gene_58    Gene_59    Gene_60 
#> 0.19236317 0.47002627 0.52742974 0.14504665 
```
