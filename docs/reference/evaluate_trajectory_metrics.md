# Full Trajectory Accuracy Evaluation

Evaluates trajectory preservation using pseudotime correlation, branch
height discrepancy, and univariate distribution distance metrics on
pseudotime. If raw count matrices are provided, pseudotime and lineage
trees are automatically inferred directly from the scRNA-seq data.

## Usage

``` r
evaluate_trajectory_metrics(
  ref_data,
  sim_data,
  cell_types_ref = NULL,
  cell_types_sim = NULL
)
```

## Arguments

- ref_data:

  Reference count matrix (genes x cells) OR numeric vector of
  pseudotimes.

- sim_data:

  Simulated count matrix (genes x cells) OR numeric vector of
  pseudotimes.

- cell_types_ref:

  Optional vector of reference cell types (for lineage tree inference).

- cell_types_sim:

  Optional vector of simulated cell types (for lineage tree inference).

## Value

A named list of trajectory accuracy metrics.

## Examples

``` r
data(example_scrna, package = "scSimEval")
evaluate_trajectory_metrics(example_scrna$ref, example_scrna$sim)
#> $pseudotime_correlation
#> [1] 1
#> 
#> $tree_height_rmse
#> [1] NA
#> 
#> $pseudotime_MAD
#> [1] 0.4272106
#> 
#> $pseudotime_KS
#> [1] 0.8
#> 
#> $pseudotime_MAE
#> [1] 0.3871033
#> 
#> $pseudotime_RMSE
#> [1] 0.3982072
#> 
#> $pseudotime_OV
#> [1] 0.2176649
#> 
#> $pseudotime_Bhattacharyya
#> [1] 0.006268418
#> 
#> $pseudotime_Wasserstein
#> [1] 0.3871033
#> 
#> $pseudotime_ECDF_DiffArea
#> [1] 0.3877229
#> 
#> $pseudotime_Runs_Statistic
#> [1] -7.772059
#> 
#> $pseudotime_Runs_PValue
#> [1] 3.861015e-15
#> 
#> $pseudotime_NN_Mismatch
#> [1] 0.80625
#> 
#> $pseudotime_Between_Dataset_Silh_Global
#> [1] 0.5144246
#> 
#> $pseudotime_Between_Dataset_Silh_Local
#> [1] 0.6821646
#> 
```
