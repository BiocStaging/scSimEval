# Run edgeR-QLF Differential Accessibility (DA) Analysis on Pseudobulk Data (Libra Engine)

Performs differential accessibility (DA) analysis for scATAC-seq or
differential expression (DE) for scRNA-seq on replicate pseudobulk count
matrices using edgeR's Quasi-Likelihood F-test (QLF) model (Squair et
al., Nature Communications 2021; Teo et al., Nature Communications 2024
/ Libra).

## Usage

``` r
run_edger_qlf_da(
  pseudobulk_mat,
  groups,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0,
  norm_method = c("TMM", "RLE", "none")
)
```

## Arguments

- pseudobulk_mat:

  Pseudobulk count matrix (features/peaks x sample_replicates).

- groups:

  Factor or character vector of comparison group conditions for each
  column in `pseudobulk_mat`.

- fdr_cutoff:

  Adjusted p-value significance threshold (default 0.05).

- logfc_cutoff:

  Absolute log2 fold-change threshold (default 0).

- norm_method:

  Normalization method for edgeR calcNormFactors: "TMM" (default),
  "RLE", or "none".

## Value

A data.frame of differential accessibility results:

- feature:

  Peak/gene identifier

- logFC:

  log2 fold change between conditions

- logCPM:

  Average log2 counts per million

- F_stat:

  Quasi-likelihood F-statistic

- PValue:

  Raw p-value from F-test

- FDR:

  Benjamini-Hochberg adjusted p-value

- is_da:

  Logical flag indicating whether the region is significantly
  differentially accessible

## Details

Fits negative binomial generalized linear models using Empirical Bayes
quasi-likelihood dispersion shrinkage and F-tests (\`edgeR::glmQLFit\`
and \`edgeR::glmQLFTest\`), offering benchmark-proven control over false
discovery rates across biological replicates.

## Examples

``` r
mat <- matrix(rpois(100 * 4, lambda = 10), nrow = 100, ncol = 4)
rownames(mat) <- paste0("Peak_", 1:100)
colnames(mat) <- c("Rep1:A", "Rep2:A", "Rep1:B", "Rep2:B")
groups <- c("A", "A", "B", "B")
res <- run_edger_qlf_da(mat, groups)
#> calcNormFactors has been renamed to normLibSizes
head(res)
#>        feature      logFC   logCPM     F_stat    PValue       FDR is_da
#> Peak_1  Peak_1 -0.3019143 13.54416 0.46741686 0.4949693 0.8971892 FALSE
#> Peak_2  Peak_2  0.1974071 13.45300 0.18414511 0.6682968 0.8971892 FALSE
#> Peak_3  Peak_3 -0.2887477 13.32175 0.37540544 0.5407685 0.8971892 FALSE
#> Peak_4  Peak_4  0.1137972 13.60187 0.07144893 0.7895141 0.8971892 FALSE
#> Peak_5  Peak_5  0.4949255 13.65736 1.42380272 0.2341921 0.7888139 FALSE
#> Peak_6  Peak_6 -0.1030907 13.95182 0.08846914 0.7664411 0.8971892 FALSE
```
