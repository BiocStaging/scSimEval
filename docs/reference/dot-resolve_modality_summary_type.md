# Resolve Summary Feature Label Based on Assay Modality

Resolve Summary Feature Label Based on Assay Modality

## Usage

``` r
.resolve_modality_summary_type(
  df,
  assay_modality = c("auto", "scRNA-seq", "scATAC-seq", "multiomics")
)
```

## Arguments

- df:

  Benchmark data frame

- assay_modality:

  Optional modality string ("auto", "scRNA-seq", "scATAC-seq",
  "multiomics")

## Value

A character scalar: `"peak"`, `"gene"`, or `"features"`
