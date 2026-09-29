# Predictive Cell Identity Classification Using DE / Top Variable Features

Trains a classifier on simulated features (80% train / 20% test) to
verify whether the simulated biological signals are sufficiently
predictive of cell identities. If `de_features` is NULL, the top
variable features are automatically selected.

## Usage

``` r
evaluate_predictive_de_model(
  data,
  group,
  de_features = NULL,
  n_top = 50,
  method = c("knn", "rf", "svm")
)
```

## Arguments

- data:

  Matrix (features x cells).

- group:

  Cell type labels.

- de_features:

  Optional character vector of selected DE feature names.

- n_top:

  Number of top features to auto-select if `de_features` is NULL
  (default 50).

- method:

  Classifier: "knn" (default), "rf", or "svm".

## Value

A list containing Accuracy, Macro Precision, Macro Recall, and Macro F1
score.

## Examples

``` r
data(example_scrna, package = "scSimEval")
evaluate_predictive_de_model(example_scrna$ref, group = example_scrna$cell_types)
#> $accuracy
#> [1] 0.8125
#> 
#> $precision
#> [1] 0.8846154
#> 
#> $recall
#> [1] 0.75
#> 
#> $F1
#> [1] 0.7681159
#> 
```
