# Plot Single-Cell Simulator Evaluation Summary

Produces a 600 DPI horizontal bar matrix ranking single-cell simulators
across the 8 canonical evaluation categories and overall composite
performance.

## Usage

``` r
plot_evaluation_summary(
  data,
  method_classes = NULL,
  category_colors = NULL,
  normalize_scores = TRUE,
  simulators_order = NULL,
  show_labels = TRUE,
  title = "Single-Cell Simulator Evaluation Summary",
  subtitle = "Average fidelity scores across 8 evaluation categories",
  base_size = 11
)

plot_benchmark_summary_bars(
  data,
  method_classes = NULL,
  category_colors = NULL,
  normalize_scores = TRUE,
  simulators_order = NULL,
  show_labels = TRUE,
  title = "Single-Cell Simulator Evaluation Summary",
  subtitle = "Average fidelity scores across 8 evaluation categories",
  base_size = 11
)

plot_summary_bars(
  data,
  method_classes = NULL,
  category_colors = NULL,
  normalize_scores = TRUE,
  simulators_order = NULL,
  show_labels = TRUE,
  title = "Single-Cell Simulator Evaluation Summary",
  subtitle = "Average fidelity scores across 8 evaluation categories",
  base_size = 11
)
```

## Arguments

- data:

  A benchmark summary data frame, a list containing
  `benchmark_summary_table`, or the output of
  [`evaluate_simulation_accuracy`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md).

- method_classes:

  Optional named vector mapping method names to class/category labels
  (e.g., `c("Splatter" = "Class 1", "scDesign3" = "Class 2")`). If
  provided, rows are grouped into class panels.

- category_colors:

  Optional named vector specifying custom colors for category panels and
  overall score.

- normalize_scores:

  Logical indicating whether to normalize scores to \[0, 1\]
  (direction-aware where larger is better). Default is `TRUE`.

- simulators_order:

  Optional vector specifying a custom ordering of simulators. Default is
  ranking by descending Overall Score.

- show_labels:

  Logical indicating whether to display numerical score labels on/beside
  each bar. Default is `TRUE`.

- title:

  Character plot title. Default is
  `"Single-Cell Simulator Evaluation Summary"`.

- subtitle:

  Character plot subtitle.

- base_size:

  Numeric base font size (default 11).

## Value

A `ggplot` object representing the multi-panel evaluation summary.

## Details

![Single-Cell Simulator Evaluation
Summary](figures/evaluation_summary_bars.png)

## Examples

``` r
bm_df <- data.frame(
  Method = rep(c("scDesign3", "scMultiSim", "Splatter", "SymSim"), each = 6),
  Category = rep(c(
    "(I) Distributional Properties", "(I) Distributional Properties",
    "(II) Sparsity & Zero-Inflation Dynamics",
    "(III) Cellular Manifold & Clustering Concordance",
    "(IV) Technical Variation & Batch Mixing Dynamics",
    "(V) Biological Signal & Differential Fidelity"
  ), 4),
  Metric = rep(c("KS", "Wasserstein", "Zero_Prob_RMSE",
                 "ARI", "CMS", "DE_F1"), 4),
  Value = c(
    0.10, 0.16, 0.07, 0.85, 0.88, 0.82,
    0.14, 0.21, 0.11, 0.79, 0.83, 0.77,
    0.23, 0.34, 0.17, 0.67, 0.72, 0.66,
    0.31, 0.44, 0.24, 0.56, 0.62, 0.55
  ),
  stringsAsFactors = FALSE
)
p <- plot_evaluation_summary(bm_df, base_size = 9.5)
print(p)
```
