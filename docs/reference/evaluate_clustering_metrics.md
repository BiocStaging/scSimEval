# Full Clustering Performance Evaluation

Evaluates both unsupervised cluster separation (ASW, DB, Neighborhood
Purity, CDI) and supervised concordance against ground truth labels
(ARI, NMI).

## Usage

``` r
evaluate_clustering_metrics(
  data,
  cluster_info = NULL,
  pred_clusters = NULL,
  dist_mat = NULL,
  cell_types = NULL,
  ref_data = NULL,
  ref_cell_types = NULL,
  batch_info = NULL,
  threads = 1L
)
```

## Arguments

- data:

  Count or normalized matrix (features x cells).

- cluster_info:

  Cluster or cell type labels.

- pred_clusters:

  Optional predicted cluster labels (if different from ground truth). If
  NULL, k-means is automatically run.

- dist_mat:

  Optional precomputed distance matrix.

- cell_types:

  Optional alias for cluster_info.

- ref_data:

  Optional reference data to compute reference clustering quality.

- ref_cell_types:

  Optional reference cell type labels (if ref_data has different cells
  than data).

- batch_info:

  Optional batch labels per cell passed to
  [`calc_cdi`](https://kabilanbio.github.io/scSimEval/reference/calc_cdi.md).

- threads:

  Integer number of parallel CPU threads passed to
  [`calc_cdi`](https://kabilanbio.github.io/scSimEval/reference/calc_cdi.md)
  (default `1L`).

## Value

A named list of clustering quality and cluster concordance metrics (keys
present depend on inputs):

- `silhouette_sim`:

  Mean silhouette width across ground-truth labels in simulated data.

- `silhouette`:

  Alias of `silhouette_sim`.

- `davies_bouldin`:

  Davies-Bouldin index computed on simulated data.

- `davies_bouldin_sim`:

  Alias of `davies_bouldin`.

- `neighborhood_purity`:

  Mean k-nearest-neighbor label purity in simulated data.

- `silhouette_ref`:

  Mean silhouette width in reference data (present when `ref_data` is
  supplied).

- `davies_bouldin_ref`:

  Davies-Bouldin index in reference data (present when `ref_data` is
  supplied).

- `cdi`:

  Normalized CDI-AIC score from
  [`calc_cdi`](https://kabilanbio.github.io/scSimEval/reference/calc_cdi.md)
  (when finite).

- `CDI_AIC`:

  Alias of `cdi`.

- `CDI_BIC`:

  Normalized CDI-BIC score from
  [`calc_cdi`](https://kabilanbio.github.io/scSimEval/reference/calc_cdi.md)
  (when finite).

- `ari`:

  Adjusted Rand Index between `pred_clusters` (or k-means clusters) and
  `cluster_info`.

- `ARI`:

  Alias of `ari`.

- `nmi`:

  Normalized Mutual Information between `pred_clusters` and
  `cluster_info`.

- `NMI`:

  Alias of `nmi`.

## Examples

``` r
data <- matrix(stats::rpois(200, lambda = 3), 20, 10)
cl <- factor(rep(c("A", "B"), each = 5))
dist_mat <- stats::dist(t(data))
evaluate_clustering_metrics(data, cl)
#> $silhouette_sim
#> [1] -0.04865282
#> 
#> $silhouette
#> [1] -0.04865282
#> 
#> $davies_bouldin
#> [1] 3.904203
#> 
#> $davies_bouldin_sim
#> [1] 3.904203
#> 
#> $neighborhood_purity
#> [1] 0.38
#> 
#> $cdi
#> [1] 4.413008
#> 
#> $CDI_AIC
#> [1] 4.413008
#> 
#> $CDI_BIC
#> [1] 4.534042
#> 
#> $ari
#> [1] -0.05882353
#> 
#> $ARI
#> [1] -0.05882353
#> 
#> $nmi
#>          1 
#> 0.03705068 
#> 
#> $NMI
#>          1 
#> 0.03705068 
#> 
```
