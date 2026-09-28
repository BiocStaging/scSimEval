# Evaluate Chromatin Accessibility Profile Concordance (DiTSim)

Computes global and cell-type-stratified Pearson and Spearman
correlations and Kullback-Leibler (KL) divergence between empirical
reference and simulated single-cell chromatin accessibility profiles.

## Usage

``` r
calc_accessibility_profile_concordance(
  ref_data,
  sim_data,
  cell_types = NULL,
  use_tfidf = TRUE
)
```

## Arguments

- ref_data:

  Reference count/accessibility matrix (peaks x cells).

- sim_data:

  Simulated count/accessibility matrix (peaks x cells).

- cell_types:

  Optional factor or vector of cell-type annotations for cells.

- use_tfidf:

  Logical, whether to apply TF-IDF transformation prior to mean profile
  calculation (default TRUE).

## Value

A list containing global PCC, global SCC, per-cell-type mean PCC/SCC,
and mean accessibility KL divergence.

## Examples

``` r
data(example_multiomics, package = "scSimEval")
r_rna <- example_multiomics$ref_multi$rna
r_atac <- example_multiomics$ref_multi$atac
calc_accessibility_profile_concordance(r_rna, r_atac)
#> $global_pcc
#> [1] 0.2982499
#> 
#> $global_scc
#> [1] 0.3483745
#> 
#> $kl_divergence
#> [1] 0.03359463
#> 
#> $mean_celltype_pcc
#> [1] NA
#> 
#> $mean_celltype_scc
#> [1] NA
#> 
#> $celltype_summary
#> NULL
#> 
```
