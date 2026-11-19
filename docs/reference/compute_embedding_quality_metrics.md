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

A `data.frame` (with `check.names = FALSE`) containing one row per
dataset with columns:

- `Dataset`:

  Dataset label (e.g., `"Reference"` or simulator name).

- `Role`:

  Dataset role (`"Reference"` or `"Simulated"`).

- `Cells (N)`:

  Number of cells in the embedding.

- `Mean Silhouette`:

  Mean silhouette width across cell-type labels rounded to 4 decimals
  (or `NA_real_`).

- `ARI (Cluster Fidelity)`:

  Adjusted Rand Index between clustering and cell-type labels rounded to
  4 decimals (or `NA_real_`).

- `Mean Library Size`:

  Mean total counts per cell rounded to 1 decimal.

- `Mean Detected Features`:

  Mean detected features per cell rounded to 1 decimal.

- `Library Size Diff (%)`:

  Relative percent difference in mean library size vs. Reference
  (present when a `"Reference"` row exists).

## Examples

``` r
data(example_scrna, package = "scSimEval")
emb <- compute_dataset_embeddings(
  reference = example_scrna$ref,
  simulated = list("Splatter" = example_scrna$sim),
  reduction = "umap",
  cell_types = example_scrna$cell_types
)
compute_embedding_quality_metrics(emb)
#>     Dataset      Role Cells (N) Mean Silhouette ARI (Cluster Fidelity)
#> 1 Reference Reference       100         -0.0391                 0.1203
#> 2  Splatter Simulated       100         -0.0191                 0.0725
#>   Mean Library Size Mean Detected Features Library Size Diff (%)
#> 1              37.6                    7.9                  0.00
#> 2              59.7                    7.1                 58.78
```
