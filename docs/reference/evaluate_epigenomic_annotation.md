# Evaluate Epigenomic Supervised Cell-Type Annotation (EpiAnno / SCAN-ATAC-Sim)

Trains a supervised classification model on reference single-cell
chromatin accessibility profiles (scATAC-seq / scCAS) and evaluates
cell-type prediction / projection accuracy on simulated cells against
ground-truth cell-type labels (Chen et al., Nat Mach Intell 2022).

## Usage

``` r
evaluate_epigenomic_annotation(
  ref_data,
  sim_data,
  ref_celltypes,
  sim_celltypes,
  method = c("knn", "centroid"),
  n_pcs = 20,
  k = 5
)
```

## Arguments

- ref_data:

  Reference count/accessibility matrix (features x cells).

- sim_data:

  Simulated count/accessibility matrix (features x cells).

- ref_celltypes:

  Factor or character vector of cell types for reference cells.

- sim_celltypes:

  Factor or character vector of true cell types for simulated cells.

- method:

  Classification method: "knn" (k-nearest neighbors on PCA, default) or
  "centroid" (nearest centroid classifier).

- n_pcs:

  Number of principal components for dimensionality reduction (default
  20).

- k:

  Number of nearest neighbors for k-NN (default 5).

## Value

A list containing overall accuracy, balanced accuracy, macro F1, macro
precision, macro recall, Cohen's kappa, and per-class performance
metrics.

## Examples

``` r
data(example_multiomics, package = "scSimEval")
evaluate_epigenomic_annotation(
  example_multiomics$ref_multi$atac,
  example_multiomics$sim_multi$atac,
  ref_celltypes = example_multiomics$cell_types,
  sim_celltypes = example_multiomics$cell_types
)
#> $accuracy
#> [1] 0.8625
#> 
#> $balanced_accuracy
#> [1] 0.8625
#> 
#> $macro_f1
#> [1] 0.8624785
#> 
#> $macro_precision
#> [1] 0.8627267
#> 
#> $macro_recall
#> [1] 0.8625
#> 
#> $cohen_kappa
#> [1] 0.725
#> 
#> $per_cell_type_f1
#>     TypeA     TypeB 
#> 0.8607595 0.8641975 
#> 
#> $confusion_matrix
#>        
#>         TypeA TypeB
#>   TypeA    34     5
#>   TypeB     6    35
#> 
#> $predicted_labels
#>  [1] "TypeA" "TypeB" "TypeA" "TypeA" "TypeA" "TypeA" "TypeA" "TypeA" "TypeA"
#> [10] "TypeB" "TypeA" "TypeA" "TypeA" "TypeA" "TypeA" "TypeA" "TypeB" "TypeA"
#> [19] "TypeB" "TypeA" "TypeA" "TypeB" "TypeA" "TypeA" "TypeA" "TypeA" "TypeA"
#> [28] "TypeA" "TypeA" "TypeA" "TypeA" "TypeA" "TypeA" "TypeA" "TypeA" "TypeB"
#> [37] "TypeA" "TypeA" "TypeA" "TypeA" "TypeB" "TypeB" "TypeB" "TypeB" "TypeB"
#> [46] "TypeA" "TypeB" "TypeB" "TypeA" "TypeB" "TypeA" "TypeB" "TypeB" "TypeB"
#> [55] "TypeB" "TypeB" "TypeA" "TypeB" "TypeB" "TypeB" "TypeB" "TypeB" "TypeB"
#> [64] "TypeB" "TypeB" "TypeB" "TypeB" "TypeB" "TypeB" "TypeB" "TypeB" "TypeA"
#> [73] "TypeB" "TypeB" "TypeB" "TypeB" "TypeB" "TypeB" "TypeB" "TypeB"
#> 
```
