# Full Clustering Performance Evaluation

Evaluates both unsupervised cluster separation (ASW, Dunn, Connectivity,
DB, CH, Neighborhood Purity, CDI) and supervised concordance against
ground truth labels (Clustering Accuracy ACC, Hungarian
F1/Precision/Recall, ARI, NMI, AMI, FMI, Homogeneity, Completeness,
V-measure).

## Usage

``` r
evaluate_clustering_metrics(
  data,
  cluster_info = NULL,
  pred_clusters = NULL,
  dist_mat = NULL,
  cell_types = NULL,
  ref_data = NULL
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

## Value

A named list of all clustering metrics.

## Examples

``` r
data <- matrix(stats::rnorm(200), 20, 10)
cl <- factor(rep(c("A", "B"), each = 5))
dist_mat <- stats::dist(data)
evaluate_clustering_metrics(data, cl)
#> Warning: NaNs produced
#> $silhouette_sim
#> [1] 0.0251254
#> 
#> $silhouette
#> [1] 0.0251254
#> 
#> $dunn_sim
#> [1] 0.5719799
#> 
#> $dunn
#> [1] 0.5719799
#> 
#> $connectivity_sim
#> [1] NA
#> 
#> $connectivity
#> [1] NA
#> 
#> $davies_bouldin_sim
#> [1] 2.529681
#> 
#> $davies_bouldin
#> [1] 2.529681
#> 
#> $calinski_harabasz_sim
#> [1] 1.249101
#> 
#> $calinski_harabasz
#> [1] 1.249101
#> 
#> $neighborhood_purity
#> [1] 0.48
#> 
#> $cdi
#> [1] -5.460536
#> 
#> $CDI_AIC
#> [1] -5.460536
#> 
#> $CDI_BIC
#> [1] -5.400018
#> 
#> $clustering_accuracy
#> [1] 0.7
#> 
#> $accuracy
#> [1] 0.7
#> 
#> $hungarian_F1
#> [1] 0.6969697
#> 
#> $hungarian_precision
#> [1] 0.7083333
#> 
#> $hungarian_recall
#> [1] 0.7
#> 
#> $ari
#> [1] 0.05970149
#> 
#> $ARI
#> [1] 0.05970149
#> 
#> $nmi
#>         1 
#> 0.1263464 
#> 
#> $NMI
#>         1 
#> 0.1263464 
#> 
#> $ami
#> [1] 0.0403207
#> 
#> $AMI
#> [1] 0.0403207
#> 
#> $fmi
#> [1] 0.48795
#> 
#> $FMI
#> [1] 0.48795
#> 
#> $homogeneity
#> [1] 0.1245112
#> 
#> $completeness
#> [1] 0.1282364
#> 
#> $v_measure
#> [1] 0.1263464
#> 
#> $matched_pairs
#>   pred truth precision recall        F1
#> 1    1     A 0.7500000    0.6 0.6666667
#> 2    2     B 0.6666667    0.8 0.7272727
#> 
```
