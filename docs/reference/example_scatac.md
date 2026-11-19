# Example scATAC-seq Benchmarking Dataset (Rice Bud)

A real, subsampled single-cell ATAC sequencing peak accessibility
dataset from rice bud tissue (*Oryza sativa*) designed for testing,
benchmarking, and demonstrating evaluation metrics in `scSimEval`.
Contains empirical reference peak counts and `scDesign3`-simulated
counts across the top 100 highly variable peaks (HVPs) and the same 100
balanced paired cells as `example_scrna`.

## Usage

``` r
data(example_scatac)
```

## Format

A list with 4 elements:

- ref:

  A `SingleCellExperiment` object containing reference peak
  accessibility counts (100 HVPs x 100 cells), cell types and batch
  annotations in `colData`, and precomputed `PCA` in `reducedDims`.

- sim:

  A `SingleCellExperiment` object containing simulated peak
  accessibility counts (100 HVPs x 100 cells), cell types and batch
  annotations in `colData`, and precomputed `PCA` in `reducedDims`.

- cell_types:

  A named factor of cell-type identities across 8 Rice Bud cell types
  for the 100 cells.

- batch_info:

  A named factor of batch annotations (`Bud_21`, `Bud_22`) for the 100
  cells.

## Source

Derived from rice bud paired multiomics data (*Oryza sativa*, GEO
accession GSE232863).

## Examples

``` r
data(example_scatac)
dim(example_scatac$ref)
#> [1] 100 100
dim(example_scatac$sim)
#> [1] 100 100
table(example_scatac$batch_info)
#> 
#> Bud_21 Bud_22 
#>     61     39 
```
