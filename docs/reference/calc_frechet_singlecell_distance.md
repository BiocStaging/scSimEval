# Frechet Single-Cell Distance (FSD)

Computes the single-cell analogue of Frechet Inception Distance (FID) on
low-dimensional PCA embeddings between reference and simulated cell
populations.

## Usage

``` r
calc_frechet_singlecell_distance(
  ref_mat,
  sim_mat,
  cells_as_cols = TRUE,
  n_pcs = 15
)
```

## Arguments

- ref_mat:

  Matrix of reference cells (cells x features or features x cells).

- sim_mat:

  Matrix of simulated cells (cells x features or features x cells).

- cells_as_cols:

  Logical, whether cells are columns (default TRUE).

- n_pcs:

  Number of principal components to evaluate (default 15).

## Value

A list with Frechet distance (FSD), mean discrepancy, and covariance
trace discrepancy.

## Examples

``` r
ref <- matrix(stats::rnorm(100), 50, 2)
sim <- matrix(stats::rnorm(100), 50, 2)
calc_frechet_singlecell_distance(ref, sim)
#> $fsd
#> [1] 12.23358
#> 
#> $fsd_squared
#> [1] 149.6605
#> 
#> $mean_discrepancy
#> [1] 55.32323
#> 
#> $cov_discrepancy
#> [1] 94.33724
#> 
```
