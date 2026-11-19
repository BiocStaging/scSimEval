# Generative Precision and Recall for Single-Cell Manifolds

Evaluates the realism (Precision) and coverage/diversity (Recall) of
generated single-cell distributions relative to empirical reference data
using k-NN hyperspheres.

## Usage

``` r
calc_generative_precision_recall(
  ref_mat,
  sim_mat,
  cells_as_cols = TRUE,
  k = 5,
  max_cells = 500,
  max_features = 1000
)
```

## Arguments

- ref_mat:

  Matrix of reference cells (cells x features or features x cells).

- sim_mat:

  Matrix of simulated cells (cells x features or features x cells).

- cells_as_cols:

  Logical, whether cells are columns (default TRUE).

- k:

  Number of nearest neighbors to define the local manifold radius
  (default 5).

- max_cells:

  Maximum number of cells to subsample (default 500).

- max_features:

  Maximum number of features to retain for speed (default 1000).

## Value

A named list summarizing k-NN hypersphere manifold overlap:

- `generative_precision`:

  Fraction of simulated cells falling within the `k`-NN radius of at
  least one reference cell.

- `generative_recall`:

  Fraction of reference cells falling within the `k`-NN radius of at
  least one simulated cell.

- `generative_f1`:

  Harmonic mean of `generative_precision` and `generative_recall`.

## Examples

``` r
ref <- matrix(stats::rnorm(200), 20, 10)
sim <- matrix(stats::rnorm(200), 20, 10)
calc_generative_precision_recall(ref, sim)
#> $generative_precision
#> [1] 1
#> 
#> $generative_recall
#> [1] 1
#> 
#> $generative_f1
#> [1] 1
#> 
```
