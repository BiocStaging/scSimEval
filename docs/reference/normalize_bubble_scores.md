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

## Examples

``` r
df <- data.frame(
  Method = rep(c("scDesign3", "Splatter", "SymSim"), each = 2),
  Metric = rep(c("KS", "ARI"), 3),
  Value  = c(0.12, 0.85, 0.28, 0.68, 0.41, 0.52),
  stringsAsFactors = FALSE
)
norm_df <- normalize_bubble_scores(df)
print(norm_df)
#>      Method Metric Value Score_Raw Score_Norm
#> 1 scDesign3     KS  0.12      0.12  1.0000000
#> 2 scDesign3    ARI  0.85      0.85  1.0000000
#> 3  Splatter     KS  0.28      0.28  0.4482759
#> 4  Splatter    ARI  0.68      0.68  0.4848485
#> 5    SymSim     KS  0.41      0.41  0.0000000
#> 6    SymSim    ARI  0.52      0.52  0.0000000
```
