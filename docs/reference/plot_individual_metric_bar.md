# Plot Individual Metric Performance Barplot

Generates a publication-quality barplot comparing simulator methods on a
single evaluation metric, with exact score labels, direction awareness,
and performance ranking.

## Usage

``` r
plot_individual_metric_bar(
  benchmark_data,
  metric,
  score_type = c("normalized", "raw"),
  palette = NULL,
  base_size = 12
)
```

## Arguments

- benchmark_data:

  A benchmark data.frame or named list containing benchmark tables.

- metric:

  Character. The name of the metric to plot (e.g. `"KS"`,
  `"Wasserstein"`, `"ARI"`).

- score_type:

  Character. Either `"normalized"` (standardized \[0, 1\] fidelity
  score; default) or `"raw"` (original unnormalized metric value).

- palette:

  Optional named character vector of simulator colors.

- base_size:

  Numeric. Base font size. Default `12`.

## Value

A `ggplot` object.

## Examples

``` r
bm_df <- data.frame(
  Method   = c("scDesign3", "scMultiSim", "MOSim", "Splatter", "SymSim"),
  Category = "(I) Distributional Properties",
  Metric   = "KS",
  Value    = c(0.10, 0.14, 0.18, 0.23, 0.31),
  stringsAsFactors = FALSE
)
p <- plot_individual_metric_bar(bm_df, metric = "KS", base_size = 10)
print(p)
```
