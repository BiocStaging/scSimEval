# Perform Likelihood Ratio Test for Comparing Nested Single-Cell Simulation Models

Direct port and generalization of
[`scDesign3::perform_lrt`](https://rdrr.io/pkg/scDesign3/man/perform_lrt.html)
(Song et al., Nat Biotechnol 2024). Performs the likelihood ratio test
to compare two nested simulation models (e.g., cell-type/covariate model
vs intercept-only null model, or spline trajectory vs linear model).

## Usage

``` r
calc_likelihood_ratio_test(
  alter_model,
  null_model,
  df_alter = NULL,
  df_null = NULL
)
```

## Arguments

- alter_model:

  Alternative model (more complex) or list of alternative models per
  gene, or numeric log-likelihoods.

- null_model:

  Null model (simpler, strictly nested) or list of null models per gene,
  or numeric log-likelihoods.

- df_alter:

  Degrees of freedom for alternative model (used if models are numeric
  log-likelihoods).

- df_null:

  Degrees of freedom for null model (used if models are numeric
  log-likelihoods).

## Value

A data.frame containing:

- LogLik_alter:

  Log-likelihood under alternative model.

- LogLik_null:

  Log-likelihood under null model.

- df_alter:

  Degrees of freedom of alternative model.

- df_null:

  Degrees of freedom of null model.

- LR_statistic:

  Likelihood ratio statistic (-2 \* (LL_null - LL_alter)).

- delta_df:

  Difference in degrees of freedom (df_alter - df_null).

- p_value:

  P-value from chi-squared test with delta_df degrees of freedom.

## Examples

``` r
calc_likelihood_ratio_test(stats::rnorm(50), stats::rnorm(50))
#>    LogLik_alter  LogLik_null df_alter df_null LR_statistic delta_df p_value
#> 1  -0.417221145 -0.135185828       NA      NA -0.564070633       NA      NA
#> 2  -0.229395923 -0.640606348       NA      NA  0.822420849       NA      NA
#> 3   0.937472402 -2.606257861       NA      NA  7.087460527       NA      NA
#> 4   0.271706851 -1.292585490       NA      NA  3.128584681       NA      NA
#> 5  -0.595170773  0.764718765       NA      NA -2.719779076       NA      NA
#> 6   1.498313908  0.854546636       NA      NA  1.287534543       NA      NA
#> 7  -0.634377790 -0.639918641       NA      NA  0.011081701       NA      NA
#> 8  -0.310062141  0.242779652       NA      NA -1.105683585       NA      NA
#> 9   1.415969444 -0.184792816       NA      NA  3.201524520       NA      NA
#> 10  1.027384325  0.435881528       NA      NA  1.183005594       NA      NA
#> 11 -1.489916488  0.457453759       NA      NA -3.894740494       NA      NA
#> 12  0.097752902  0.096573751       NA      NA  0.002358303       NA      NA
#> 13 -0.274708562  0.658762688       NA      NA -1.866942500       NA      NA
#> 14 -0.821170042  0.638641991       NA      NA -2.919624065       NA      NA
#> 15 -0.988737886 -0.121423381       NA      NA -1.734629010       NA      NA
#> 16  1.135918365 -1.508831879       NA      NA  5.289500488       NA      NA
#> 17  1.089650728  1.837303392       NA      NA -1.495305329       NA      NA
#> 18 -0.124177823  0.863936442       NA      NA -1.976228530       NA      NA
#> 19 -0.444504065  0.652696617       NA      NA -2.194401364       NA      NA
#> 20  0.045719649 -0.008536259       NA      NA  0.108511816       NA      NA
#> 21 -1.651363708 -0.224889813       NA      NA -2.852947790       NA      NA
#> 22 -0.329350529  1.281120043       NA      NA -3.220941143       NA      NA
#> 23  0.168828394  0.546208001       NA      NA -0.754759214       NA      NA
#> 24 -0.434508733 -0.847064275       NA      NA  0.825111085       NA      NA
#> 25 -0.998421828  0.383729104       NA      NA -2.764301865       NA      NA
#> 26  0.194344736 -0.961418464       NA      NA  2.311526399       NA      NA
#> 27 -0.427944031  0.665976373       NA      NA -2.187840808       NA      NA
#> 28 -0.111651217  1.754012957       NA      NA -3.731328346       NA      NA
#> 29 -1.980097643 -0.138019777       NA      NA -3.684155731       NA      NA
#> 30  0.009597104 -1.927338579       NA      NA  3.873871365       NA      NA
#> 31 -0.062404642 -0.159374712       NA      NA  0.193940140       NA      NA
#> 32 -0.296573292 -0.709438106       NA      NA  0.825729629       NA      NA
#> 33 -1.393751041  0.618534679       NA      NA -4.024571439       NA      NA
#> 34 -0.468199511  1.490079899       NA      NA -3.916558821       NA      NA
#> 35  0.036386642 -0.278114794       NA      NA  0.629002872       NA      NA
#> 36 -0.298738168 -0.300270279       NA      NA  0.003064222       NA      NA
#> 37  0.295046284 -0.308454776       NA      NA  1.207002120       NA      NA
#> 38 -0.690588407  0.556190278       NA      NA -2.493557372       NA      NA
#> 39 -0.972875414  0.640017144       NA      NA -3.225785116       NA      NA
#> 40 -0.094946186 -1.801798590       NA      NA  3.413704808       NA      NA
#> 41 -1.882785636 -0.736105850       NA      NA -2.293359571       NA      NA
#> 42  0.329577842  1.571742307       NA      NA -2.484328929       NA      NA
#> 43  0.102973499  1.432726985       NA      NA -2.659506972       NA      NA
#> 44 -0.748540920 -1.622148705       NA      NA  1.747215571       NA      NA
#> 45 -0.035888525 -0.042089873       NA      NA  0.012402696       NA      NA
#> 46 -1.369178436  0.466922587       NA      NA -3.672202045       NA      NA
#> 47  0.212756151  0.761926923       NA      NA -1.098341546       NA      NA
#> 48 -0.887887107 -1.322229972       NA      NA  0.868685730       NA      NA
#> 49 -0.358430879 -0.323562015       NA      NA -0.069737728       NA      NA
#> 50 -1.025289292 -0.246355790       NA      NA -1.557867005       NA      NA
```
