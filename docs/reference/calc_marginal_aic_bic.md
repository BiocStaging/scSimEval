# Marginal Model Goodness of Fit and Information Criteria for Single-Cell Simulators

Evaluates gene-wise and aggregate marginal goodness-of-fit (Poisson,
Negative Binomial, or Gaussian) for simulated or fitted single-cell
count matrices, computing total and mean AIC and BIC across genes
(scDesign3; Song et al., 2024).

## Usage

``` r
calc_marginal_aic_bic(
  counts,
  fitted_means,
  dispersions = NULL,
  distribution = c("poisson", "nb", "gaussian"),
  n_params_per_gene = NULL
)
```

## Arguments

- counts:

  Matrix or data.frame of observed or simulated counts (genes x cells).

- fitted_means:

  Matrix or data.frame of model fitted means or expectations (genes x
  cells).

- dispersions:

  Optional vector of gene-level dispersion parameters for Negative
  Binomial model. If NULL, estimated via method-of-moments.

- distribution:

  Parametric distribution: "poisson", "nb" (Negative Binomial), or
  "gaussian".

- n_params_per_gene:

  Number of estimated parameters per gene. If NULL, defaults to 1 for
  Poisson, 2 for NB, and 2 for Gaussian.

## Value

A list containing:

- total_loglik:

  Sum of marginal log-likelihoods across all genes and cells.

- total_aic:

  Aggregate AIC across all genes.

- total_bic:

  Aggregate BIC across all genes.

- mean_gene_aic:

  Mean AIC per gene.

- mean_gene_bic:

  Mean BIC per gene.

- gene_loglik:

  Vector of log-likelihoods per gene.

- gene_aic:

  Vector of AIC per gene.

- gene_bic:

  Vector of BIC per gene.

## Examples

``` r
calc_marginal_aic_bic(stats::rnorm(50), stats::rnorm(50))
#> $total_loglik
#> [1] 3.69921
#> 
#> $total_aic
#> [1] 92.60158
#> 
#> $total_bic
#> [1] -7.398419
#> 
#> $mean_gene_aic
#> [1] 1.852032
#> 
#> $mean_gene_bic
#> [1] -0.1479684
#> 
#> $gene_loglik
#>  [1] -26.1851034  -0.6910669  -1.7282698  -1.9540175  -4.4133334  11.6399656
#>  [7]   7.3789649   7.1207077  26.8902097  14.0419244  -0.1154853 -13.6697367
#> [13]  -3.5918324 -23.7036849  -3.5800488  -1.9615111  -1.3957709  -5.5900839
#> [19]  -3.0909186 -19.0917742  32.4137527   2.8331813  -1.8290634 -43.5217262
#> [25]  -3.4582560  -0.9881907  -2.3156039  -0.7528467   1.1362431   2.5525351
#> [31]  -0.2169932  12.2240701  14.8001893 -12.9214260  18.9000073  25.4761761
#> [37] -28.6642438 -41.8550556  14.4320399   0.4633477  -1.5488071  19.4345210
#> [43]  -3.3559335  -0.9259248  -1.4541922  17.8605151  14.9022871  17.4935149
#> [49]  -1.1615942  -2.5624482
#> 
#> $gene_aic
#>  [1]  54.3702069   3.3821338   5.4565396   5.9080350  10.8266668 -21.2799312
#>  [7] -12.7579297 -12.2414154 -51.7804194 -26.0838488   2.2309705  29.3394734
#> [13]   9.1836648  49.4073698   9.1600977   5.9230223   4.7915418  13.1801679
#> [19]   8.1818371  40.1835484 -62.8275054  -3.6663626   5.6581267  89.0434524
#> [25]   8.9165120   3.9763814   6.6312078   3.5056933  -0.2724862  -3.1050702
#> [31]   2.4339863 -22.4481402 -27.6003786  27.8428520 -35.8000146 -48.9523522
#> [37]  59.3284876  85.7101113 -26.8640797   1.0733045   5.0976142 -36.8690421
#> [43]   8.7118669   3.8518495   4.9083845 -33.7210301 -27.8045743 -32.9870298
#> [49]   4.3231885   7.1248964
#> 
#> $gene_bic
#>  [1]  52.3702069   1.3821338   3.4565396   3.9080350   8.8266668 -23.2799312
#>  [7] -14.7579297 -14.2414154 -53.7804194 -28.0838488   0.2309705  27.3394734
#> [13]   7.1836648  47.4073698   7.1600977   3.9230223   2.7915418  11.1801679
#> [19]   6.1818371  38.1835484 -64.8275054  -5.6663626   3.6581267  87.0434524
#> [25]   6.9165120   1.9763814   4.6312078   1.5056933  -2.2724862  -5.1050702
#> [31]   0.4339863 -24.4481402 -29.6003786  25.8428520 -37.8000146 -50.9523522
#> [37]  57.3284876  83.7101113 -28.8640797  -0.9266955   3.0976142 -38.8690421
#> [43]   6.7118669   1.8518495   2.9083845 -35.7210301 -29.8045743 -34.9870298
#> [49]   2.3231885   5.1248964
#> 
```
