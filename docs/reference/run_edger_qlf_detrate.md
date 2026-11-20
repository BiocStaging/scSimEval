# Run edgeR QLF Differential Expression with Cellular Detection Rate (edgeRQLFDetRate)

Implements the benchmark-optimal `edgeRQLFDetRate` algorithm (Soneson &
Robinson, 2018; Duo et al., 2024 / Simpipe) for single-cell differential
expression analysis. Fits a quasi-likelihood negative binomial
generalized linear model including the cellular detection rate (DetRate:
fraction of expressed genes per cell) as a covariate.

## Usage

``` r
run_edger_qlf_detrate(
  counts,
  group,
  cdr = NULL,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0,
  min_cells = 3
)
```

## Arguments

- counts:

  Single-cell count matrix (genes x cells), or `SingleCellExperiment`
  object.

- group:

  Factor or vector specifying the comparison cell groups (must contain
  exactly 2 levels for pairwise DE).

- cdr:

  Optional numeric vector of cellular detection rates per cell. If NULL,
  computed automatically as fraction of non-zero genes.

- fdr_cutoff:

  FDR / adjusted p-value threshold for significance (default 0.05).

- logfc_cutoff:

  Absolute log2 fold-change threshold (default 0).

- min_cells:

  Minimum number of cells expressing a gene to test (default 3).

## Value

A data.frame with columns:

- gene:

  Gene identifier

- logFC:

  log2 fold change

- logCPM:

  log2 counts per million

- F_stat:

  Quasi-likelihood F-statistic (or t-statistic in fallback)

- PValue:

  Raw p-value

- FDR:

  Benjamini-Hochberg adjusted p-value

- is_deg:

  Logical flag indicating whether the gene meets the significance
  threshold

## Examples

``` r
data(example_scrna, package = "scSimEval")
top2 <- names(sort(table(example_scrna$cell_types), decreasing = TRUE))[1:2]
idx <- which(example_scrna$cell_types %in% top2)
res <- run_edger_qlf_detrate(example_scrna$ref[, idx], example_scrna$cell_types[idx])
head(res)
#>                      gene      logFC   logCPM      F_stat       PValue
#> Os04g0188433 Os04g0188433  0.3993855 15.55724  0.06888415 0.7930375570
#> Os04g0496300 Os04g0496300 -3.0963319 16.29357  5.06806129 0.0246475243
#> Os04g0659300 Os04g0659300  2.8811191 15.67757  2.85801360 0.0913182157
#> Os01g0756900 Os01g0756900 -4.5981140 16.23515 11.56257739 0.0007070579
#> Os01g0910900 Os01g0910900  2.7969064 15.44732  3.20103952 0.0739784000
#> Os10g0377150 Os10g0377150 -2.7209843 15.95230  4.30062802 0.0384250890
#>                      FDR is_deg
#> Os04g0188433 0.817068998  FALSE
#> Os04g0496300 0.055867722  FALSE
#> Os04g0659300 0.119096056  FALSE
#> Os01g0756900 0.003746631   TRUE
#> Os01g0910900 0.114330254  FALSE
#> Os10g0377150 0.076850178  FALSE
```
