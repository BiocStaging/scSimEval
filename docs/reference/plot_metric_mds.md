# Multi-Dimensional Scaling (MDS) Ordination of Evaluation Metrics or Simulators

Projects evaluation metric profiles or simulator performances into a 2D
MDS space for benchmarking single-cell simulation methods.

## Usage

``` r
plot_metric_mds(
  benchmark_data,
  ordination_by = c("simulators", "summaries"),
  assay_modality = c("auto", "scRNA-seq", "scATAC-seq", "multiomics"),
  group_by = c("category", "summary_type"),
  by_category = FALSE,
  as_list = FALSE,
  category = NULL,
  exclude_categories = c("(VII) Trajectory & Lineage Topology",
    "(VIII) Computational Scalability & Efficiency"),
  metrics = NULL,
  normalize_scores = TRUE,
  palette = NULL,
  title = NULL,
  base_size = 11,
  ncol = 2
)
```

## Arguments

- benchmark_data:

  Benchmark summary data.frame, matrix, or list.

- ordination_by:

  Character. Either `"simulators"` (ordinated points represent candidate
  simulators) or `"summaries"` (default; ordinated points represent
  biological and computational summaries).

- assay_modality:

  Character string specifying the molecular modality: `"auto"` (default;
  automatically detects from data attributes or metrics), `"scRNA-seq"`,
  `"scATAC-seq"`, or `"multiomics"`. When `"scATAC-seq"`, feature-level
  metrics are labeled `"peak"`. When `"scRNA-seq"`, feature-level
  metrics are labeled `"gene"`. When `"multiomics"`, feature-level
  metrics are labeled `"features"`.

- group_by:

  Character. Grouping variable for metric points when
  `ordination_by = "summaries"`: `"category"` (default) groups points by
  evaluation category (I–VIII); `"summary_type"` groups by
  feature/cell/global level.

- by_category:

  Logical. If `TRUE` and `ordination_by = "summaries"`, generates a
  multi-panel grid displaying an MDS ordination for each evaluation
  category separately (skipping categories with fewer than 3 measures).
  Default is `FALSE`.

- as_list:

  Logical. If `TRUE` and `by_category = TRUE`, returns a named list of
  ggplot objects. Default `FALSE`.

- category:

  Optional character vector to filter by evaluation category (e.g.
  `"(I) Distributional Properties"`).

- exclude_categories:

  Character vector of categories to exclude from the ordination. Default
  is
  `c("(VII) Trajectory & Lineage Topology", "(VIII) Computational Scalability & Efficiency")`
  because these two categories contain only a few measures each, leaving
  the 6 canonical categories.

- metrics:

  Optional character vector of specific metrics to include.

- normalize_scores:

  Logical. If `TRUE` (default), scores are direction-aware normalized to
  \[0, 1\] so that metrics on disparate scales (e.g. MB vs distance) do
  not artificially distort Euclidean distance.

- palette:

  Optional named character vector of colors. For
  `ordination_by = "summaries"`, maps feature-level (`"peak"`, `"gene"`,
  or `"features"`), `"cell"`, and `"global"`.

- title:

  Optional plot title.

- base_size:

  Numeric. Base font size. Default `11`.

- ncol:

  Integer. Number of columns when `by_category = TRUE`. Default `2`.

## Value

A `ggplot` or `patchwork` object.

## Details

When `ordination_by = "summaries"`, each point represents an evaluated
property/metric, colored by its biological level: `"gene"` (red),
`"peak"` (red, for scATAC-seq), `"features"` (red, for multiomics),
`"cell"` (blue), or `"global"` (green). If `by_category = TRUE`, a
multi-panel dashboard is produced showing one MDS ordination per
category (by default across the 6 canonical categories with \>= 3
measures, excluding Trajectory and Scalability which each contain only 2
measures).

## Examples

``` r
bm_df <- data.frame(
  Method = rep(c("scDesign3", "scMultiSim", "Splatter", "SymSim", "MOSim"), each = 8),
  Category = rep(c(
    "(I) Distributional Properties", "(I) Distributional Properties",
    "(I) Distributional Properties",
    "(II) Sparsity & Zero-Inflation Dynamics",
    "(III) Cellular Manifold & Clustering Concordance",
    "(III) Cellular Manifold & Clustering Concordance",
    "(IV) Technical Variation & Batch Mixing Dynamics",
    "(V) Biological Signal & Differential Fidelity"
  ), 5),
  Metric = rep(c("KS", "Wasserstein", "RMSE", "Zero_Prob_RMSE",
                 "ARI", "Silhouette_ASW", "CMS", "DE_F1"), 5),
  Value = c(
    0.10, 0.16, 0.22, 0.07, 0.85, 0.74, 0.88, 0.84,
    0.18, 0.20, 0.24, 0.15, 0.82, 0.81, 0.76, 0.79,
    0.21, 0.35, 0.31, 0.12, 0.64, 0.70, 0.80, 0.61,
    0.34, 0.42, 0.48, 0.26, 0.59, 0.52, 0.61, 0.68,
    0.15, 0.29, 0.26, 0.19, 0.75, 0.66, 0.70, 0.77
  ),
  stringsAsFactors = FALSE
)
p <- plot_metric_mds(bm_df, base_size = 9.5)
print(p)
```
