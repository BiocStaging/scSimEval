# Calculate Cross-Modality Correlation Fidelity

Evaluates whether cross-modality feature relationships (e.g.
peak-to-gene links, promoter accessibility vs. gene expression, or mRNA
vs. surface protein) are accurately captured by the simulation method
compared to the empirical reference.

## Usage

``` r
calc_cross_modality_correlation(
  ref_mod1,
  ref_mod2,
  sim_mod1,
  sim_mod2,
  feature_pairs = NULL,
  method = c("spearman", "pearson")
)
```

## Arguments

- ref_mod1:

  Reference matrix for Modality 1 (features x cells).

- ref_mod2:

  Reference matrix for Modality 2 (features x cells).

- sim_mod1:

  Simulated matrix for Modality 1.

- sim_mod2:

  Simulated matrix for Modality 2.

- feature_pairs:

  Optional 2-column data.frame/matrix of paired feature names/indices.

- method:

  Correlation method: "spearman" (default) or "pearson".

## Value

A named list of the 7 univariate accuracy metrics comparing the
reference vs. simulated cross-modality correlation distributions.

## Examples

``` r
data(example_multiomics, package = "scSimEval")
m <- example_multiomics
calc_cross_modality_correlation(m$ref_multi$rna, m$ref_multi$atac,
                                m$sim_multi$rna, m$sim_multi$atac)
#> $cross_modality_cor_MAD
#> [1] 0.06236202
#> 
#> $cross_modality_cor_KS
#> [1] 0.3166667
#> 
#> $cross_modality_cor_MAE
#> [1] 0.06352891
#> 
#> $cross_modality_cor_RMSE
#> [1] 0.06626403
#> 
#> $cross_modality_cor_OV
#> [1] 0.7340243
#> 
#> $cross_modality_cor_Bhattacharyya
#> [1] 0.001537345
#> 
#> $cross_modality_cor_Wasserstein
#> [1] 0.06352891
#> 
#> $cross_modality_cor_ECDF_DiffArea
#> [1] 0.1199979
#> 
#> $cross_modality_cor_Runs_Statistic
#> [1] -3.300231
#> 
#> $cross_modality_cor_Runs_PValue
#> [1] 0.0004830262
#> 
#> $cross_modality_cor_NN_Mismatch
#> [1] 0.1166667
#> 
#> $cross_modality_cor_Between_Dataset_Silh_Global
#> [1] 0.05927015
#> 
#> $cross_modality_cor_Between_Dataset_Silh_Local
#> [1] 0.1521508
#> 
```
