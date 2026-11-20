# Plot Benchmark Metric Performance Barplots by Category

Generates publication-ready faceted barplots comparing simulator methods
across all individual evaluation metrics belonging to a specified
benchmark category. Each metric is presented in its own sub-panel with
exact numeric score labels, direction-awareness indicators (indicating
whether higher or lower values are optimal), and standardized or
original raw score formulation.

## Usage

``` r
plot_category_metric_bars(
  benchmark_data,
  category,
  score_type = c("normalized", "raw"),
  palette = NULL,
  ncol = NULL,
  base_size = 11
)
```

## Arguments

- benchmark_data:

  A benchmark data.frame (e.g. `demo$benchmark_summary_table` or output
  from
  [`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md))
  or a named list of benchmark result tables.

- category:

  Character. The benchmark category to visualize (e.g.,
  `"(I) Distributional Properties"` or `"Distributional Properties"`).

- score_type:

  Character. Either `"normalized"` (standardized \[0, 1\] fidelity score
  where 1.0 is optimal; default) or `"raw"` (unnormalized original
  metric values in native measurement units).

- palette:

  Optional named character vector of simulator colors.

- ncol:

  Integer. Optional number of columns in facet layout. Defaults to
  dynamic layout based on metric count.

- base_size:

  Numeric. Base font size. Default `11`.

## Value

A `ggplot` object.

## Examples

``` r
bm_df <- data.frame(
  Method   = rep(c("scDesign3", "scMultiSim", "Splatter", "SymSim"), each = 4),
  Category = "(I) Distributional Properties",
  Metric   = rep(c("KS", "Wasserstein", "MAD", "RMSE"), 4),
  Value    = c(
    0.10, 0.16, 0.12, 0.19,
    0.15, 0.22, 0.16, 0.25,
    0.23, 0.34, 0.26, 0.36,
    0.31, 0.44, 0.35, 0.47
  ),
  stringsAsFactors = FALSE
)
p <- plot_category_metric_bars(
  bm_df,
  category = "(I) Distributional Properties",
  base_size = 9.5
)
print(p)
```
