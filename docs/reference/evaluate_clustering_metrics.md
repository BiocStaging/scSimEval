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
#> [1] -0.01727743
#> 
#> $silhouette
#> [1] -0.01727743
#> 
#> $dunn_sim
#> [1] 0.5508905
#> 
#> $dunn
#> [1] 0.5508905
#> 
#> $connectivity_sim
#> [1] NA
#> 
#> $connectivity
#> [1] NA
#> 
#> $davies_bouldin_sim
#> [1] 3.22902
#> 
#> $davies_bouldin
#> [1] 3.22902
#> 
#> $calinski_harabasz_sim
#> [1] 0.7656334
#> 
#> $calinski_harabasz
#> [1] 0.7656334
#> 
#> $neighborhood_purity
#> [1] 0.44
#> 
#> $cdi
#> [1] -6.167219
#> 
#> $CDI_AIC
#> [1] -6.167219
#> 
#> $CDI_BIC
#> [1] -6.106702
#> 
#> $clustering_accuracy
#> [1] 0.6
#> 
#> $accuracy
#> [1] 0.6
#> 
#> $hungarian_F1
#> [1] 0.5833333
#> 
#> $hungarian_precision
#> [1] 0.6190476
#> 
#> $hungarian_recall
#> [1] 0.6
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
#> $ami
#> [1] 0
#> 
#> $AMI
#> [1] 0
#> 
#> $fmi
#> [1] 0.4564355
#> 
#> $FMI
#> [1] 0.4564355
#> 
#> $homogeneity
#> [1] 0.03485155
#> 
#> $completeness
#> [1] 0.03954603
#> 
#> $v_measure
#> [1] 0.03705068
#> 
#> $matched_pairs
#>   pred truth precision recall        F1
#> 1    1     B 0.5714286    0.8 0.6666667
#> 2    2     A 0.6666667    0.4 0.5000000
#> 
```
