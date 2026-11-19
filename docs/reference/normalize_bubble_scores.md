# Normalize Bubble Matrix Scores Across Benchmark Datasets and Simulators

Rescales metric raw scores to a 0-1 range where 1 is always best, taking
into account whether the metric is higher-is-better or lower-is-better,
naturally bounded, or unbounded.

## Usage

``` r
normalize_bubble_scores(df)
```

## Arguments

- df:

  A data frame of benchmark results containing at least `Metric` and
  `Score_Raw` (or `Score`).

## Value

The input data frame with `Score_Norm` populated in `[0, 1]` (higher is
always better).
