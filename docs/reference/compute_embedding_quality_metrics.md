# Compute Quantitative Quality Metrics for Low-Dimensional Cell Embeddings

Computes quantitative cluster separability (Silhouette width),
unsupervised cluster recovery (Adjusted Rand Index), and distribution
preservation metrics (Mean UMI, Mean detected features, sparsity
percentage, and relative percentage error) between reference and
simulated datasets.

## Usage

``` r
compute_embedding_quality_metrics(
  embedding_data,
  metric_coords = c("Dim1", "Dim2")
)
```

## Arguments

- embedding_data:

  A `data.frame` produced by
  [`compute_dataset_embeddings`](https://kabilanbio.github.io/scSimEval/reference/compute_dataset_embeddings.md).

- metric_coords:

  Character vector of the 2 coordinate columns used to compute distances
  (default: `c("Dim1", "Dim2")`).

## Value

A `data.frame` summarizing quantitative embedding quality metrics for
each dataset.

## Examples

``` r
data(example_scrna, package = "scSimEval")
emb <- compute_dataset_embeddings(
  reference = example_scrna$ref,
  simulated = list("Splatter" = example_scrna$sim),
  reduction = "umap",
  cell_types = example_scrna$cell_types
)
#> Warning: You're computing too large a percentage of total singular values, use a standard svd instead.
#> Warning: You're computing too large a percentage of total singular values, use a standard svd instead.
compute_embedding_quality_metrics(emb)
#>     Dataset      Role Cells (N) Mean Silhouette ARI (Cluster Fidelity)
#> 1 Reference Reference        80          0.0009                 0.0277
#> 2  Splatter Simulated        80          0.0074                 0.0188
#>   Mean Library Size Mean Detected Features Library Size Diff (%)
#> 1             256.3                   49.2                  0.00
#> 2             242.3                   48.3                  5.46
```
