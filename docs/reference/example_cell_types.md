# Example Cell-Type Annotations Vector (Rice Bud)

A named factor vector containing discrete cell-type labels across all 8
Rice Bud cell types for 100 balanced single cells (12–13 cells per
type), matching the column names of
[`example_scrna`](https://kabilanbio.github.io/scSimEval/reference/example_scrna.md)
and
[`example_scatac`](https://kabilanbio.github.io/scSimEval/reference/example_scatac.md).

## Usage

``` r
data(example_cell_types)
```

## Format

A named factor with 100 elements across 8 levels (`Bud-unknown`,
`Bud_axis`, `Bud_primordium`, `Epidermis`, `Parenchyma_cell`,
`Root_meristem`, `Shoot_apical_meristem`, `Vascular_cylinder`).

## Source

Derived from `example_scrna$cell_types`.

## Examples

``` r
data(example_cell_types)
head(example_cell_types)
#> AAACCGGCAAGTGTTT-1 AAACGGATCGATCAGT-1 AACAGCAAGGTTAGAG-1 AACCTAATCAAGTGAG-1 
#>    Parenchyma_cell    Parenchyma_cell          Epidermis     Bud_primordium 
#> AAGGTGCAGGGTCCAC-1 AAGTTTGTCCACCTGT-1 
#>  Vascular_cylinder      Root_meristem 
#> 8 Levels: Bud_axis Bud_primordium Bud-unknown Epidermis ... Vascular_cylinder
table(example_cell_types)
#> example_cell_types
#>              Bud_axis        Bud_primordium           Bud-unknown 
#>                    13                    13                    13 
#>             Epidermis       Parenchyma_cell         Root_meristem 
#>                    13                    12                    12 
#> Shoot_apical_meristem     Vascular_cylinder 
#>                    12                    12 
```
