# Validate and Standardize Multiomics Input Objects

Validates that multiomics inputs conform to `SingleCellExperiment`
structures. Accepts either a single multi-modal `SingleCellExperiment`
object with an `altExp` for the secondary modality, a named list of two
`SingleCellExperiment` objects (e.g.
`list(RNA = rna_sce, ATAC = atac_sce)`), or a file path to an `.rds`
file containing either structure.

## Usage

``` r
validate_multiomics_input(
  obj,
  arg_name = "multi_data",
  reduction_name = NULL,
  cell_type_col = "cell_type",
  batch_col = "batch"
)
```

## Arguments

- obj:

  Multiomics `SingleCellExperiment` object, list of two
  `SingleCellExperiment` objects, or `.rds` path.

- arg_name:

  Name of argument for diagnostic error messages.

- reduction_name:

  Optional reduction name to extract.

- cell_type_col:

  Column name in `colData` containing cell-type labels (default
  `"cell_type"`).

- batch_col:

  Column name in `colData` containing batch labels (default `"batch"`).

## Value

A named list with 9 elements:

- m1:

  Count matrix for modality 1 (typically RNA).

- m2:

  Count matrix for modality 2 (typically ATAC).

- m1_obj:

  `SingleCellExperiment` object for modality 1.

- m2_obj:

  `SingleCellExperiment` object for modality 2.

- m1_name:

  Character name of modality 1.

- m2_name:

  Character name of modality 2.

- cell_types:

  Character vector (paired) or named list (unpaired) of cell-type
  labels, or `NULL`.

- batch_info:

  Character vector (paired) or named list (unpaired) of batch labels, or
  `NULL`.

- embedding:

  Extracted dimensional reduction matrix (cells x dims), or `NULL`.
