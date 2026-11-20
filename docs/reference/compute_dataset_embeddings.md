# Compute Low-Dimensional Embeddings for Reference and Simulated Single-Cell Datasets

Projects empirical biological reference and simulated count matrices
into low-dimensional representations (UMAP, t-SNE, or PCA) for
side-by-side and multi-method topological benchmarking. Automatically
performs library-size normalization, variance-based feature selection,
principal component analysis, and optional unsupervised clustering.

## Usage

``` r
compute_dataset_embeddings(
  reference,
  simulated,
  reduction = c("umap", "tsne", "pca"),
  n_pcs = 30,
  perplexity = 30,
  n_neighbors = 15,
  min_dist = 0.3,
  seed = 42,
  cell_types = NULL,
  batch = NULL
)
```

## Arguments

- reference:

  Biological reference `SingleCellExperiment` object or count matrix
  (features x cells).

- simulated:

  A single simulated count matrix/`SingleCellExperiment` or a named list
  of them.

- reduction:

  Character string specifying the dimensionality reduction method:
  `"umap"` (default), `"tsne"`, or `"pca"`.

- n_pcs:

  Integer specifying the number of principal components to calculate
  (default: 30).

- perplexity:

  Numeric perplexity for t-SNE (default: 30; automatically adapted for
  small sample sizes).

- n_neighbors:

  Integer number of nearest neighbors for UMAP (default: 15).

- min_dist:

  Numeric minimum distance parameter for UMAP (default: 0.3).

- seed:

  Random seed for reproducibility (default: 42).

- cell_types:

  Optional factor or character vector of cell type annotations for
  reference cells (or named list if per-dataset).

- batch:

  Optional factor or character vector of batch annotations.

## Value

A tidy `data.frame` containing cell coordinates (`Dim1`, `Dim2`),
`Dataset` name, `Dataset_Type` ("Reference" vs "Simulated"),
`Cell_Type`, `Cluster`, `Library_Size`, `Detected_Features`, and
`Batch`.

## Examples

``` r
data(example_scrna, package = "scSimEval")
emb <- compute_dataset_embeddings(
  reference = example_scrna$ref,
  simulated = list("Splatter" = example_scrna$sim),
  reduction = "umap",
  cell_types = example_scrna$cell_types
)
head(emb)
#>              Cell_ID        Dim1       Dim2   Dataset      Role Dataset_Type
#> 1 AAACCGGCAAGTGTTT-1 -0.49065503  0.3989727 Reference Reference    Reference
#> 2 AAACGGATCGATCAGT-1  1.10728518  1.1790269 Reference Reference    Reference
#> 3 AACAGCAAGGTTAGAG-1 -1.10536298  2.2994250 Reference Reference    Reference
#> 4 AACCTAATCAAGTGAG-1 -0.04908856 -3.0515676 Reference Reference    Reference
#> 5 AAGGTGCAGGGTCCAC-1 -2.96960934 -0.3985537 Reference Reference    Reference
#> 6 AAGTTTGTCCACCTGT-1 -0.91385993  1.4372928 Reference Reference    Reference
#>           Cell_Type   Cluster Library_Size Detected_Features   Batch
#> 1   Parenchyma_cell Cluster_3            6                 5 Batch 1
#> 2   Parenchyma_cell Cluster_4            2                 2 Batch 1
#> 3         Epidermis Cluster_8           18                11 Batch 1
#> 4    Bud_primordium Cluster_4           17                11 Batch 1
#> 5 Vascular_cylinder Cluster_5          271                14 Batch 1
#> 6     Root_meristem Cluster_3           10                 8 Batch 1
```
