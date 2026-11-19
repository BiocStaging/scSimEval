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
#> calcNormFactors has been renamed to normLibSizes
head(res)
#>                      gene      logFC   logCPM      F_stat       PValue
#> Os04g0188433 Os04g0188433 -1.6973682 14.99692  0.93341220 3.417600e-01
#> Os04g0496300 Os04g0496300 -5.4762732 17.57997 21.01082531 4.433661e-05
#> Os04g0659300 Os04g0659300 -0.6178729 15.19454  0.05069294 8.231181e-01
#> Os01g0756900 Os01g0756900 -7.0743346 17.54568 35.66166925 3.414173e-07
#> Os01g0910900 Os01g0910900  0.4332316 14.93205  0.21850884 6.440717e-01
#> Os10g0377150 Os10g0377150 -5.1144155 17.05487 17.09207411 2.055816e-04
#>                       FDR is_deg
#> Os04g0188433 4.841600e-01  FALSE
#> Os04g0496300 1.884306e-04   TRUE
#> Os04g0659300 9.328672e-01  FALSE
#> Os01g0756900 5.804095e-06   TRUE
#> Os01g0910900 8.422476e-01  FALSE
#> Os10g0377150 6.989775e-04   TRUE
```
