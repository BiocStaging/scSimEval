# Extract Dataset Dimensional Properties and Summary Statistics

Computes key dimensional and biological properties from single-cell or
multiomics count matrices, including cell counts, feature counts,
sparsity percentage, biological group/cell type counts, batch counts,
median library size, median detected features per cell, and mean
expression level.

## Usage

``` r
extract_dataset_summary(
  mat,
  role = "Reference",
  method_name = "Dataset",
  modality = "scRNA-seq",
  cell_types = NULL,
  batch_info = NULL
)
```

## Arguments

- mat:

  A count matrix (genes/features x cells), dgCMatrix,
  SingleCellExperiment, or Seurat object.

- role:

  Character string describing the role of the dataset (e.g.
  `"Biological Reference"` or `"Simulated"`). Default is `"Reference"`.

- method_name:

  Character string indicating the method or dataset name (e.g.
  `"Empirical Reference"`, `"Splatter"`). Default is `"Dataset"`.

- modality:

  Character string describing the molecular modality (e.g.
  `"scRNA-seq"`, `"scATAC-seq"`). Default is `"scRNA-seq"`.

- cell_types:

  Optional factor or character vector of cell type annotations.

- batch_info:

  Optional factor or character vector of batch annotations.

## Value

A 1-row `data.frame` (with `check.names = FALSE`) summarizing the
dataset with columns:

- `Dataset / Simulator`:

  Character label from `method_name`.

- `Role`:

  Dataset role label from `role`.

- `Modality`:

  Molecular modality label from `modality`.

- `Cells (N)`:

  Formatted total number of cells.

- `Features (P)`:

  Formatted total number of features.

- `Sparsity`:

  Percentage of zero entries formatted to two decimal places.

- `Cell Types (Groups)`:

  Formatted count of unique non-missing cell types, or
  `"Not specified"`.

- `Batches`:

  Formatted count of unique non-missing batches (defaults to `"1"`).

- `Median Lib Size`:

  Formatted median total counts per cell.

- `Median Detected Features`:

  Formatted median number of detected (non-zero) features per cell.

- `Mean Expression`:

  Formatted overall mean count across all matrix entries.

## Examples

``` r
data(example_scrna, package = "scSimEval")
extract_dataset_summary(example_scrna$ref, role = "Reference", method_name = "Empirical")
#>   Dataset / Simulator      Role  Modality Cells (N) Features (P) Sparsity
#> 1           Empirical Reference scRNA-seq       100          100   92.06%
#>   Cell Types (Groups) Batches Median Lib Size Median Detected Features
#> 1                   8       2            15.0                        7
#>   Mean Expression
#> 1           0.376
```
