# Evaluate Multiple Single-Cell Datasets Simultaneously

Orchestrates simultaneous benchmarking across a collection of
single-cell `SingleCellExperiment` datasets, including paired/unpaired
multiomics and unimodal scRNA-seq or scATAC-seq experiments.

## Usage

``` r
evaluate_multiple_datasets(
  datasets,
  pair_by_prefix = TRUE,
  compute_bivariate = FALSE,
  threads = 1L,
  verbose = TRUE
)
```

## Arguments

- datasets:

  A named list of dataset specifications (`Format A` grouped
  multiomics/unimodal lists, or `Format B` flat modality-suffixed
  lists).

- pair_by_prefix:

  Logical; if `TRUE` (default), datasets sharing a root prefix across
  RNA and ATAC entries are paired for multiomics evaluation.

- compute_bivariate:

  Logical; whether to compute 2D bivariate tests (default `FALSE`).

- threads:

  Integer number of CPU threads. Default is `1L`.

- verbose:

  Logical; whether to print progress messages. Default is `TRUE`.

## Value

An S3 object of class `c("scSimEval_consolidated", "list")` containing:

- consolidated_summary_table:

  Master tidy data.frame with columns `Dataset`, `Data_Name`,
  `Modality`, `Category`, `Property`, `Metric`, and `Value`.

- consolidated_score_matrix:

  Numeric matrix comparing datasets across all evaluated metric keys.

- dataset_overview:

  Data frame summarizing metric counts and mean KS, Wasserstein, and
  RMSE per dataset.

- dataset_results:

  Named list of detailed per-dataset evaluation outputs.

## Examples

``` r
data(example_scrna, package = "scSimEval")
ds <- list("DemoCohort" = list(ref = example_scrna$ref, sim = example_scrna$sim))
res <- evaluate_multiple_datasets(ds, verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> calcNormFactors has been renamed to normLibSizes
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
```
