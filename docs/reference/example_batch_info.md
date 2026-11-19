# Example Batch Annotations Vector (Rice Bud)

A named factor vector containing biological batch annotations (`Bud_21`,
`Bud_22`) for 100 single cells, matching the column names of
[`example_scrna`](https://kabilanbio.github.io/scSimEval/reference/example_scrna.md)
and
[`example_scatac`](https://kabilanbio.github.io/scSimEval/reference/example_scatac.md).

## Usage

``` r
data(example_batch_info)
```

## Format

A named factor with 100 elements across 2 levels: `Bud_21` and `Bud_22`.

## Source

Derived from `example_scrna$batch_info`.

## Examples

``` r
data(example_batch_info)
head(example_batch_info)
#> AAACCGGCAAGTGTTT-1 AAACGGATCGATCAGT-1 AACAGCAAGGTTAGAG-1 AACCTAATCAAGTGAG-1 
#>             Bud_21             Bud_21             Bud_21             Bud_21 
#> AAGGTGCAGGGTCCAC-1 AAGTTTGTCCACCTGT-1 
#>             Bud_21             Bud_21 
#> Levels: Bud_21 Bud_22
table(example_batch_info)
#> example_batch_info
#> Bud_21 Bud_22 
#>     61     39 
```
