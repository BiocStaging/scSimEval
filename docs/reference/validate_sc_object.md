# Validate and Standardize Single-Cell Input Objects

Strictly enforces `SingleCellExperiment` compliance in `scSimEval`.
Validates that input datasets are `SingleCellExperiment` objects,
extracting underlying count matrices, cell annotations, batch labels,
and dimensional reductions. Users starting from `Seurat` objects should
first convert via
[`Seurat::as.SingleCellExperiment()`](https://satijalab.org/seurat/reference/as.SingleCellExperiment.html).

## Usage

``` r
validate_sc_object(
  obj,
  arg_name = "data",
  reduction_name = NULL,
  modality = NULL,
  cell_type_col = "cell_type",
  batch_col = "batch"
)
```

## Arguments

- obj:

  A `SingleCellExperiment` object to validate.

- arg_name:

  Name of the argument being validated (for clear error messages).

- reduction_name:

  Optional name of precomputed dimensional reduction (e.g. `"PCA"`,
  `"integrated"`).

- modality:

  Optional character string specifying the modality (e.g. `"scRNA-seq"`
  or `"scATAC-seq"`) to guide assay/`altExp` selection.

- cell_type_col:

  Column name in `colData(obj)` containing cell-type labels (default
  `"cell_type"`).

- batch_col:

  Column name in `colData(obj)` containing batch/replicate labels
  (default `"batch"`).

## Value

A named list with 5 elements:

- object:

  The validated `SingleCellExperiment` object.

- counts:

  Extracted raw count matrix (`dgCMatrix` or numeric matrix, features x
  cells).

- cell_types:

  Character vector of extracted cell-type annotations, or `NULL` if
  absent.

- batch_info:

  Character vector of extracted batch annotations, or `NULL` if absent.

- embedding:

  Extracted dimensional reduction numeric matrix (cells x dims), or
  `NULL`.
