# Fraction of Samples Closer Than The True Match (FOSCTTM)

Computes the Fraction of Samples Closer Than The True Match (FOSCTTM) to
quantify cell alignment error across single-cell multi-omics modalities
in a shared latent or integrated space (Zhai et al., Genome Biology
2024; Liu et al., Nat Biotechnol 2023). A score of 0 represents perfect
alignment (where the true paired cell is the nearest neighbor), while
0.5 corresponds to random chance.

## Usage

``` r
calc_foscttm(x, y, metric = c("euclidean", "cosine"))
```

## Arguments

- x:

  Coordinates matrix for Modality 1 (cells x dimensions).

- y:

  Coordinates matrix for Modality 2 (cells x dimensions), where row i of
  y corresponds to the true match of row i of x.

- metric:

  Distance metric to evaluate: "euclidean" (default) or "cosine".

## Value

A list containing:

- foscttm:

  Bidirectional mean FOSCTTM score across all cells (lower is better, 0
  to 0.5).

- foscttm_xy:

  Directional FOSCTTM from Modality 1 to Modality 2.

- foscttm_yx:

  Directional FOSCTTM from Modality 2 to Modality 1.

- match_at_1:

  Top-1 match rate (proportion of cells where the true match is rank 1).

- match_at_5:

  Top-5 match rate (proportion of cells where the true match is in top
  5).

- cell_foscttm:

  Vector of bidirectional FOSCTTM scores for individual cells.

## Examples

``` r
data(example_multiomics, package = "scSimEval")
r_rna <- example_multiomics$ref_multi$rna
r_atac <- example_multiomics$ref_multi$atac
calc_foscttm(r_rna, r_atac)
#> $foscttm
#> [1] 0.4636111
#> 
#> $foscttm_xy
#> [1] 0.4411111
#> 
#> $foscttm_yx
#> [1] 0.4861111
#> 
#> $match_at_1
#> [1] 0.01666667
#> 
#> $match_at_5
#> [1] 0.1
#> 
#> $cell_foscttm
#>    Gene_01    Gene_02    Gene_03    Gene_04    Gene_05    Gene_06    Gene_07 
#> 0.70833333 0.32500000 0.64166667 0.42500000 0.51666667 0.30833333 0.60833333 
#>    Gene_08    Gene_09    Gene_10    Gene_11    Gene_12    Gene_13    Gene_14 
#> 0.36666667 0.53333333 0.42500000 0.46666667 0.46666667 0.44166667 0.57500000 
#>    Gene_15    Gene_16    Gene_17    Gene_18    Gene_19    Gene_20    Gene_21 
#> 0.36666667 0.41666667 0.29166667 0.43333333 0.41666667 0.64166667 0.57500000 
#>    Gene_22    Gene_23    Gene_24    Gene_25    Gene_26    Gene_27    Gene_28 
#> 0.13333333 0.36666667 0.77500000 0.90833333 0.50000000 0.53333333 0.20833333 
#>    Gene_29    Gene_30    Gene_31    Gene_32    Gene_33    Gene_34    Gene_35 
#> 0.24166667 0.55000000 0.58333333 0.49166667 0.55833333 0.48333333 0.40000000 
#>    Gene_36    Gene_37    Gene_38    Gene_39    Gene_40    Gene_41    Gene_42 
#> 0.24166667 0.59166667 0.50833333 0.28333333 0.41666667 0.55833333 0.10833333 
#>    Gene_43    Gene_44    Gene_45    Gene_46    Gene_47    Gene_48    Gene_49 
#> 0.32500000 0.41666667 0.19166667 0.55000000 0.25833333 0.44166667 0.66666667 
#>    Gene_50    Gene_51    Gene_52    Gene_53    Gene_54    Gene_55    Gene_56 
#> 0.77500000 0.78333333 0.55833333 0.90000000 0.45000000 0.35833333 0.45833333 
#>    Gene_57    Gene_58    Gene_59    Gene_60 
#> 0.36666667 0.45833333 0.05833333 0.40833333 
#> 
```
