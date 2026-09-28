# Evaluate Gene Co-Expression Module Fidelity (ESCO)

Automatically discovers or evaluates co-expression modules and measures
the preservation of intra-module vs. inter-module correlation structure.

## Usage

``` r
calc_coexpression_module_fidelity(
  ref_mat,
  sim_mat,
  modules = NULL,
  n_modules = 5,
  top_genes = 200
)
```

## Arguments

- ref_mat:

  Reference gene expression matrix (genes x cells).

- sim_mat:

  Simulated gene expression matrix (genes x cells).

- modules:

  Optional named list of gene modules. If NULL, auto-detected via
  hierarchical clustering.

- n_modules:

  Number of modules to detect if modules is NULL (default 5).

- top_genes:

  Number of top variable genes to consider for clustering (default 200).

## Value

A list summarizing module correlation preservation, RMSE, and modularity
ratio fidelity.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_coexpression_module_fidelity(example_scrna$ref, example_scrna$sim)
#> $module_correlation_r
#> [1] -0.01457616
#> 
#> $module_correlation_rmse
#> [1] 0.1690945
#> 
#> $module_correlation_mae
#> [1] 0.1360054
#> 
#> $ref_modularity_ratio
#> [1] 0.5568422
#> 
#> $sim_modularity_ratio
#> [1] -0.02798568
#> 
#> $modularity_fidelity
#> [1] 0
#> 
```
