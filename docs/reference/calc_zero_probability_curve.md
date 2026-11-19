# Fit Zero-Probability Dropout Curve (Splatter & ZINB-WaVE)

Fits an empirical logistic dropout curve: logit(P(Y=0)) = beta_0 +
beta_1 \* log(mu) to model the scRNA-seq dropout relationship with mean
expression.

## Usage

``` r
calc_zero_probability_curve(counts)
```

## Arguments

- counts:

  Count matrix (genes x cells) or SingleCellExperiment.

## Value

A named list with the fitted logistic dropout model parameters:

- `intercept`:

  Estimated logistic regression intercept (`beta_0`), or `NA_real_` if
  fewer than 10 valid genes.

- `slope`:

  Estimated logistic regression slope (`beta_1`) on `log(mean)`.

- `midpoint`:

  Inflection point (`-beta_0 / beta_1`) on the log-mean scale.

- `r_squared`:

  Coefficient of determination (\\R^2\\) of the linear fit.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_zero_probability_curve(example_scrna$ref)
#> $intercept
#> [1] 1.349168
#> 
#> $slope
#> [1] -0.7125713
#> 
#> $midpoint
#> [1] 1.893379
#> 
#> $r_squared
#> [1] 0.8012893
#> 
```
