# Detect Differential Expression (DE) Genes (Mean Shift via Limma or T-Test)

Detect Differential Expression (DE) Genes (Mean Shift via Limma or
T-Test)

## Usage

``` r
calc_signal_de(exprs_mat, cell_types, p_sig = 0.05)
```

## Arguments

- exprs_mat:

  Log-normalized matrix (genes x cells).

- cell_types:

  Factor or binary vector of 2 cell types.

- p_sig:

  Significance threshold (default 0.05).

## Value

Vector of adjusted p-values.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_signal_de(example_scrna$ref, example_scrna$cell_types)
#>  [1] 0.7659767 0.2011222 0.2011222 0.5012172 0.2041681 0.9124780 0.2390949
#>  [8] 0.7926715 0.6940252 0.2046114 0.2011222 0.6347146 0.1367366 0.1367366
#> [15] 0.2011222 0.9475222 0.7926715 0.6940252 0.7926715 0.3992100 0.5974588
#> [22] 0.3422234 0.6940252 0.4909899 0.5313905 0.9569914 0.8547439 0.8547439
#> [29] 0.7926715 0.2390949 0.7876928 0.6940252 0.9569914 0.4037275 0.4909899
#> [36] 0.6347146 0.5904349 0.6940252 0.6347146 0.8547439 0.8547439 0.5313905
#> [43] 0.7208180 0.8547439 0.7926715 0.4037275 0.8547439 0.8734837 0.8547439
#> [50] 0.9569914 0.3992100 0.2390949 0.4909899 0.8547439 0.2390949 0.9569914
#> [57] 0.6940252 0.2296293 0.2082643 0.7876928
```
