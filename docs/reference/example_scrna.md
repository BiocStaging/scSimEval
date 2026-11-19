# Example scRNA-seq Benchmarking Dataset (Rice Bud)

A real, subsampled single-cell RNA sequencing dataset from rice bud
tissue (*Oryza sativa*) designed for testing, benchmarking, and
demonstrating evaluation metrics in `scSimEval`. Contains empirical
reference counts and `scDesign3`-simulated counts across the top 100
highly variable genes (HVGs) and 100 balanced cells covering all 8
biological cell types and 2 batches.

## Usage

``` r
data(example_scrna)
```

## Format

A list with 4 elements:

- ref:

  A `SingleCellExperiment` object containing reference scRNA-seq counts
  (100 HVGs x 100 cells), cell types and batch annotations in `colData`,
  and precomputed `PCA` in `reducedDims`.

- sim:

  A `SingleCellExperiment` object containing simulated scRNA-seq counts
  (100 HVGs x 100 cells), cell types and batch annotations in `colData`,
  and precomputed `PCA` in `reducedDims`.

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
data(example_scrna)
dim(example_scrna$ref)
#> [1] 100 100
dim(example_scrna$sim)
#> [1] 100 100
table(example_scrna$cell_types)
#> 
#>              Bud_axis        Bud_primordium           Bud-unknown 
#>                    13                    13                    13 
#>             Epidermis       Parenchyma_cell         Root_meristem 
#>                    13                    12                    12 
#> Shoot_apical_meristem     Vascular_cylinder 
#>                    12                    12 
```
