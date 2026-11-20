# Full Category VII Trajectory & Lineage Topology Evaluation (dynverse Workflow)

Evaluates trajectory and lineage preservation using the standardized
*dynverse* benchmarking framework (*Saelens et al., Nature
Biotechnology, 2019*):

1.  Standardized cell outlier QC (3 MADs), low-expression gene filtering
    (\\\ge 5\\\\ cells, mean \\\> 0.02\\), and HVG selection (empirical
    mean-variance decomposition).

2.  Shared HVG PCA projection and Docker-free Minimum Spanning Tree
    (MST) milestone trajectory inference with continuous
    cell-to-milestone geodesic distances.

3.  Five complementary trajectory metrics: `him`
    (Hamming-Ipsen-Mikhailov similarity), `edge_flip` (minimum edge flip
    similarity), `geodesic_correlation`, `F1_branches`, and
    `wcor_features` (ranger Random Forest feature importance weighted
    correlation with `num_trees = 10000` and 1 percent `mtry`).

## Usage

``` r
evaluate_trajectory_metrics(
  ref_data,
  sim_data,
  cell_types_ref = NULL,
  cell_types_sim = NULL,
  num_trees = 10000L,
  nmads = 3,
  min_cell_frac = 0.05,
  min_mean_expr = 0.02,
  fdr_threshold = 0.05,
  bio_var_threshold = 0.5,
  threads = 1L,
  seed = 42L,
  verbose = FALSE
)
```

## Arguments

- ref_data:

  Reference count matrix or `SingleCellExperiment` object.

- sim_data:

  Simulated count matrix or `SingleCellExperiment` object.

- cell_types_ref:

  Optional vector of reference cell types.

- cell_types_sim:

  Optional vector of simulated cell types.

- num_trees:

  Number of ranger trees per milestone for `wcor_features` (default
  `10000L`).

- nmads:

  MAD cutoff for outlier cell QC (default `3`).

- min_cell_frac:

  Minimum cell fraction per gene (default `0.05`).

- min_mean_expr:

  Minimum average gene count (default `0.02`).

- fdr_threshold:

  FDR threshold for HVG selection (default `0.05`).

- bio_var_threshold:

  Biological variance threshold for HVG selection (default `0.5`).

- threads:

  Number of CPU threads for ranger (default `1L`).

- seed:

  Random seed (default `42L`).

- verbose:

  Logical; print progress messages.

## Value

A named list of Category VII *dynverse* trajectory accuracy metrics:

- him:

  Hamming-Ipsen-Mikhailov topology similarity (`[0, 1]`).

- edge_flip:

  Edge flip topology similarity (`[0, 1]`).

- geodesic_correlation:

  Spearman correlation of aligned geodesic distances (`[0, 1]`).

- F1_branches:

  Macro F1 overlap of MST branch assignments (`[0, 1]`).

- wcor_features:

  Reference-weighted Pearson correlation of Random Forest gene
  importances (`[0, 1]`).

- cor_features:

  Unweighted Pearson correlation of Random Forest gene importances.

- n_hvgs:

  Number of HVGs retained for trajectory evaluation.

- ref_trajectory:

  Inferred reference milestone trajectory object.

- sim_trajectory:

  Inferred simulated milestone trajectory object.

- preprocessing_summary:

  Cell and HVG filtering summary list.

## Examples

``` r
data(example_scrna, package = "scSimEval")
res <- evaluate_trajectory_metrics(
  example_scrna$ref, example_scrna$sim,
  cell_types_ref = example_scrna$cell_types,
  cell_types_sim = example_scrna$cell_types,
  num_trees = 200L
)
unlist(res[c("him", "edge_flip", "geodesic_correlation",
             "F1_branches", "wcor_features", "n_hvgs")])
#>                  him            edge_flip geodesic_correlation 
#>            0.5893282            0.4000000            0.5980102 
#>          F1_branches        wcor_features               n_hvgs 
#>            0.3948717            0.4157592           47.0000000 
```
