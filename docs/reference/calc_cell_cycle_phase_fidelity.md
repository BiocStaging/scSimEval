# Evaluate Cell Cycle Phase Distribution Fidelity

Compares cell cycle phase proportions (e.g. G0/G1, S, G2/M) between
reference and simulated single-cell populations.

## Usage

``` r
calc_cell_cycle_phase_fidelity(ref_phases, sim_phases)
```

## Arguments

- ref_phases:

  Factor or character vector of cell cycle phase assignments in
  reference.

- sim_phases:

  Factor or character vector of cell cycle phase assignments in
  simulation.

## Value

A list with phase proportions, Jensen-Shannon divergence, and
chi-squared test p-value.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_cell_cycle_phase_fidelity(example_scrna$cell_types, example_scrna$cell_types)
#> $ref_proportions
#> TypeA TypeB 
#>   0.5   0.5 
#> 
#> $sim_proportions
#> TypeA TypeB 
#>   0.5   0.5 
#> 
#> $phase_mae
#> [1] 0
#> 
#> $phase_jsd
#> [1] 0
#> 
#> $homogeneity_pvalue
#> [1] 1
#> 
```
