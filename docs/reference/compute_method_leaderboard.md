# Compute Simulator Method Performance Leaderboard

Computes the overall performance leaderboard ranking simulation methods
across all benchmark evaluation metrics. Metrics are first
direction-inverted and min-max standardized into fidelity scores in \[0,
1\] (where 1.0 represents best observed performance). Simulators are
then rank-ordered by mean overall fidelity score.

## Usage

``` r
compute_method_leaderboard(benchmark_data)
```

## Arguments

- benchmark_data:

  A benchmark data.frame (such as `demo$benchmark_summary_table` or
  output from
  [`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md))
  or a named list of benchmark result tables.

## Value

A `data.frame` with columns:

- `Overall_Rank`: Integer rank (1 = top-performing simulator).

- `Method`: Simulator method name.

- `Average_Fidelity`: Formatted percentage string (e.g., `"58.5%"`).

- `Fidelity_Score`: Numeric composite fidelity score in \[0, 1\] rounded
  to 4 decimals.

## Examples

``` r
bm_df <- data.frame(
  Method = rep(c("scDesign3", "scMultiSim", "Splatter", "SymSim"), each = 4),
  Category = rep(c(
    "(I) Distributional Properties", "(I) Distributional Properties",
    "(III) Cellular Manifold & Clustering Concordance",
    "(IV) Technical Variation & Batch Mixing Dynamics"
  ), 4),
  Metric = rep(c("KS", "Wasserstein", "ARI", "CMS"), 4),
  Value = c(
    0.10, 0.16, 0.85, 0.88,
    0.14, 0.21, 0.79, 0.83,
    0.23, 0.34, 0.67, 0.72,
    0.31, 0.44, 0.56, 0.62
  ),
  stringsAsFactors = FALSE
)
lb <- compute_method_leaderboard(bm_df)
print(lb)
#>   Overall_Rank     Method Average_Fidelity Fidelity_Score
#> 3            1  scDesign3            75.0%         0.7500
#> 4            2 scMultiSim            65.4%         0.6541
#> 1            3   Splatter            43.3%         0.4332
#> 2            4     SymSim            25.0%         0.2500
```
