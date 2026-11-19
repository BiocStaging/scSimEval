# Example Paired Multiomics Benchmarking Dataset (Rice Bud)

A unified paired single-cell multiomics dataset combining
`example_scrna` (100 HVGs) and `example_scatac` (100 HVPs) across the
exact same 100 balanced rice bud cells, accompanied by shared cell-type
labels, batch annotations, and empirical computational scalability
metrics.

## Usage

``` r
data(example_multiomics)
```

## Format

A list with 5 elements:

- ref_multi:

  Named list of `SingleCellExperiment` objects: `rna` (100 HVGs x 100
  cells) and `atac` (100 HVPs x 100 cells), complete with `colData` and
  precomputed `PCA` in `reducedDims`.

- sim_multi:

  Named list of simulated `SingleCellExperiment` objects: `rna` (100
  x 100) and `atac` (100 x 100).

- cell_types:

  A named factor of shared cell-type identities across 8 Rice Bud cell
  types for the 100 paired cells.

- batch_info:

  A named factor of batch annotations (`Bud_21`, `Bud_22`) for the 100
  paired cells.

- resource_stats:

  A list of empirical scalability metrics: `elapsed_time`, `memory_mb`,
  `sim_cores`.

## Source

Derived from rice bud paired multiomics data (*Oryza sativa*, GEO
accession GSE232863).

## Examples

``` r
data(example_multiomics)
names(example_multiomics$ref_multi)
#> [1] "rna"  "atac"
names(example_multiomics$sim_multi)
#> [1] "rna"  "atac"
```
