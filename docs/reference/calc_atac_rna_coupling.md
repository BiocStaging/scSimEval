# Chromatin Accessibility to RNA Regulatory Effect Coupling

Evaluates whether simulated paired scATAC-seq and scRNA-seq profiles
capture the true regulatory coupling where regional chromatin opening
gates target gene transcription (scMultiSim; Li et al., Nat Methods
2023).

## Usage

``` r
calc_atac_rna_coupling(atac_data, rna_data, linked_pairs = NULL)
```

## Arguments

- atac_data:

  Matrix of chromatin peak accessibilities (peaks x cells).

- rna_data:

  Matrix of gene expression counts (genes x cells).

- linked_pairs:

  Optional 2-column data.frame of known linked peak-gene pairs. If NULL,
  assumes 1-to-1 matching by row index.

## Value

A list containing mean coupling correlation, positive coupling ratio,
and mean R-squared.

## Examples

``` r
data(example_multiomics, package = "scSimEval")
r_rna <- example_multiomics$ref_multi$rna
r_atac <- example_multiomics$ref_multi$atac
calc_atac_rna_coupling(r_rna, r_atac)
#> $mean_coupling_cor
#> [1] 0.05053426
#> 
#> $median_coupling_cor
#> [1] 0.05450756
#> 
#> $positive_coupling_ratio
#> [1] 0.7
#> 
#> $mean_r_squared
#> [1] 0.01392181
#> 
#> $coupling_correlations
#>  [1]  0.0312845269  0.1368683965  0.0406160087  0.1460827275  0.0301611844
#>  [6]  0.1108771568 -0.1339658616  0.0290878130  0.0897693201  0.1918758033
#> [11]  0.0986173172  0.1688360062  0.0234875888 -0.0166480296  0.2157394943
#> [16]  0.1913605228  0.1488709820  0.1206082132 -0.0934819349  0.0076918670
#> [21] -0.1168458056  0.1034349937  0.1685482917 -0.2398086524 -0.0335838915
#> [26] -0.0726869023 -0.0138701627  0.1120988700  0.1053941276  0.0159061794
#> [31] -0.0550171161  0.2320145448  0.0881769999 -0.1814790757  0.1619686757
#> [36]  0.0344227017 -0.2376761892  0.1180201262  0.0587262523  0.1385930997
#> [41]  0.1643753303  0.0502888627  0.0946995555 -0.0725053307 -0.0004954917
#> [46]  0.0440690631  0.1595968730  0.0301951830  0.1907721401 -0.0370862423
#> [51]  0.0268144452 -0.0127408267 -0.0935032115  0.1065337917 -0.0143846380
#> [56] -0.0015464598  0.1449270021  0.1334874707  0.0838144647  0.1106676300
#> 
```
