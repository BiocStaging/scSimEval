# Evaluate Differentially Expressed Gene (DEG) Fidelity

Evaluates the preservation of differentially expressed genes (DEGs) and
biological signal between empirical reference and simulated datasets
across three landmark single-cell benchmarking frameworks:

- **Simpipe (Duo et al., 2024):** True DEG ratio, Distribution score
  (p-value uniformity via Pearson Chi-Square goodness-of-fit test on
  remaining genes after DEG removal), and supervised machine learning
  classification (Accuracy, Precision, Recall, F1) using simulated DEGs
  to predict cell identity.

- **SimBench (Cao et al., 2021):** Symmetric Mean Absolute Percentage
  Error (SMAPE) on DEG proportions, SimBench DE fidelity score, log2
  fold-change effect size Pearson and Spearman correlation, and top-N
  DEG Jaccard overlap.

- **Shaky Foundations (Crowell et al., 2023):** Group silhouette
  separation width, silhouette discrepancy, and Percent Variance
  Explained (PVE) by group assignment.

## Usage

``` r
evaluate_deg_fidelity(
  ref_data,
  sim_data,
  ref_celltypes,
  sim_celltypes,
  group1 = NULL,
  group2 = NULL,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0.5,
  top_n_de = 100,
  classifier = c("knn", "svm", "rf"),
  run_pvalue_uniformity = TRUE
)
```

## Arguments

- ref_data:

  Reference single-cell expression or count matrix (features x cells).

- sim_data:

  Simulated single-cell expression or count matrix (features x cells).

- ref_celltypes:

  Factor or character vector of cell types for reference cells.

- sim_celltypes:

  Factor or character vector of cell types for simulated cells.

- group1:

  Optional character string specifying group 1 name. Default NULL
  (auto-selects top abundant type).

- group2:

  Optional character string specifying group 2 name. Default NULL
  (auto-selects second most abundant type).

- fdr_cutoff:

  Adjusted p-value significance threshold. Default is 0.05.

- logfc_cutoff:

  Absolute log2 fold-change cutoff. Default is 0.5.

- top_n_de:

  Number of top ranked DEGs to evaluate for Jaccard overlap and
  classification. Default is 100.

- classifier:

  Supervised classifier: "knn" (default), "svm", or "rf".

- run_pvalue_uniformity:

  Logical, whether to run Simpipe's Chi-square p-value uniformity test.
  Default is TRUE.

## Value

A list containing:

- `deg_summary_table`: Comprehensive data.frame uniting all 15 metrics
  across Simpipe, SimBench, and Shaky Foundations.

- `deg_genes_ref`: Character vector of significant DEGs in reference.

- `deg_genes_sim`: Character vector of significant DEGs in simulation.

- `contrast`: String describing the evaluated group contrast.

- `logfc_ref`: Named numeric vector of reference log2 fold-changes.

- `logfc_sim`: Named numeric vector of simulated log2 fold-changes.

- `ml_classification`: Detailed classifier performance metrics.

- `pvalue_uniformity`: Chi-square goodness-of-fit test results.

## Examples

``` r
data(example_scrna, package = "scSimEval")
evaluate_deg_fidelity(example_scrna$ref, example_scrna$sim,
                      ref_celltypes = example_scrna$cell_types,
                      sim_celltypes = example_scrna$cell_types)
#> $deg_summary_table
#>                                   Framework                     Metric
#> 1                Simpipe (Duo et al., 2024)                  DEG_Ratio
#> 2                Simpipe (Duo et al., 2024)    PValue_Uniformity_Chisq
#> 3                Simpipe (Duo et al., 2024)         Distribution_Score
#> 4                Simpipe (Duo et al., 2024)        Classifier_Accuracy
#> 5                Simpipe (Duo et al., 2024)        Classifier_Macro_F1
#> 6                Simpipe (Duo et al., 2024)    Classifier_Macro_Recall
#> 7               SimBench (Cao et al., 2021)             SimBench_SMAPE
#> 8               SimBench (Cao et al., 2021) SimBench_DE_Fidelity_Score
#> 9               SimBench (Cao et al., 2021)        Log2FC_Pearson_Corr
#> 10              SimBench (Cao et al., 2021)       Log2FC_Spearman_Corr
#> 11              SimBench (Cao et al., 2021)    Top_DEG_Jaccard_Overlap
#> 12 Shaky Foundations (Crowell et al., 2023)             Silhouette_Sim
#> 13 Shaky Foundations (Crowell et al., 2023)     Silhouette_Discrepancy
#> 14 Shaky Foundations (Crowell et al., 2023)              PVE_Group_Sim
#> 15 Shaky Foundations (Crowell et al., 2023)            PVE_Discrepancy
#>         Value                                      Direction
#> 1  0.00000000                 Target = 1.0 (Optimal balance)
#> 2  8.00000000            Lower is better (closer to Uniform)
#> 3  1.00000000          Higher is better (1.0 = Uniform null)
#> 4  0.56250000        Higher is better (Group predictability)
#> 5  0.54655870                 Higher is better (Balanced F1)
#> 6  0.56250000                 Higher is better (Sensitivity)
#> 7  0.00000000             Lower is better (0.0 = zero error)
#> 8  1.00000000         Higher is better (1.0 = perfect match)
#> 9  0.22342423           Higher is better (Effect size match)
#> 10 0.21994999            Higher is better (Rank order match)
#> 11 1.00000000            Higher is better (Set intersection)
#> 12 0.01251165            Higher is better (Group separation)
#> 13 0.01672510       Lower is better (Preserves ref geometry)
#> 14 0.02717485                         Target ~ Reference PVE
#> 15 0.02090127 Lower is better (Preserves variance explained)
#> 
#> $deg_genes_ref
#> character(0)
#> 
#> $deg_genes_sim
#> character(0)
#> 
#> $contrast
#> [1] "TypeA vs TypeB"
#> 
#> $logfc_ref
#>     Gene_01     Gene_02     Gene_03     Gene_04     Gene_05     Gene_06 
#> -0.60637028 -0.29796278 -0.50452220 -2.31778333 -1.55309970  0.31769193 
#>     Gene_07     Gene_08     Gene_09     Gene_10     Gene_11     Gene_12 
#> -1.99123826 -1.50367977 -0.34437550 -0.55530844 -1.79516352  0.72415535 
#>     Gene_13     Gene_14     Gene_15     Gene_16     Gene_17     Gene_18 
#> -2.10870020 -2.26939347 -0.97881556  0.36870164 -0.36217332  1.49073810 
#>     Gene_19     Gene_20     Gene_21     Gene_22     Gene_23     Gene_24 
#>  0.58308124  0.41152149  2.37362792  1.85517047 -0.72881718  1.68441295 
#>     Gene_25     Gene_26     Gene_27     Gene_28     Gene_29     Gene_30 
#>  1.43737516 -0.01609177  0.58626822 -1.31435042  0.81561993  2.58656260 
#>     Gene_31     Gene_32     Gene_33     Gene_34     Gene_35     Gene_36 
#> -0.32363607  0.56013593  0.76828399  0.99379745 -0.19140583 -0.15726218 
#>     Gene_37     Gene_38     Gene_39     Gene_40     Gene_41     Gene_42 
#>  0.38752957 -0.16621122 -0.56841979  0.19971306  0.56360424 -0.72796963 
#>     Gene_43     Gene_44     Gene_45     Gene_46     Gene_47     Gene_48 
#> -0.61769661 -0.02490899 -1.46089961  0.82490276 -0.59907429  1.21762555 
#>     Gene_49     Gene_50     Gene_51     Gene_52     Gene_53     Gene_54 
#> -0.15212162  1.50959792 -1.83635434  2.76079252  1.08418629  0.69835960 
#>     Gene_55     Gene_56     Gene_57     Gene_58     Gene_59     Gene_60 
#> -1.19379914 -0.18456258 -0.85965997 -2.42344584  1.18381091 -0.07550262 
#> 
#> $logfc_sim
#>     Gene_01     Gene_02     Gene_03     Gene_04     Gene_05     Gene_06 
#>  0.20098688 -1.54962939 -2.20141611 -2.51167683 -1.39182712 -3.16418085 
#>     Gene_07     Gene_08     Gene_09     Gene_10     Gene_11     Gene_12 
#>  1.00706753 -0.45864794 -1.47302161 -0.33737362 -1.88340364 -1.67832193 
#>     Gene_13     Gene_14     Gene_15     Gene_16     Gene_17     Gene_18 
#> -2.27706281  2.36600321 -0.47130545  1.76857445 -0.78360988  0.43650258 
#>     Gene_19     Gene_20     Gene_21     Gene_22     Gene_23     Gene_24 
#>  1.66269550 -1.79751652 -0.92217045  0.92912292  0.11788119  2.34336152 
#>     Gene_25     Gene_26     Gene_27     Gene_28     Gene_29     Gene_30 
#>  1.17430924 -0.72467111  2.14475554 -0.20251115  0.31885515 -0.36959198 
#>     Gene_31     Gene_32     Gene_33     Gene_34     Gene_35     Gene_36 
#>  1.01761989  0.15765163 -0.82733615  1.41685213 -1.47575878  0.92888116 
#>     Gene_37     Gene_38     Gene_39     Gene_40     Gene_41     Gene_42 
#>  0.03081620 -1.19994730  0.07958345 -0.61287496  0.45192353  1.02685338 
#>     Gene_43     Gene_44     Gene_45     Gene_46     Gene_47     Gene_48 
#>  1.49253585  3.63529833 -0.09558368 -0.71368698  1.15139639 -0.04050831 
#>     Gene_49     Gene_50     Gene_51     Gene_52     Gene_53     Gene_54 
#> -1.61361316 -0.63260072 -1.31050468  1.62231995  0.35737120 -0.16798107 
#>     Gene_55     Gene_56     Gene_57     Gene_58     Gene_59     Gene_60 
#>  0.75307391  1.42370854 -0.12720692 -0.39370357  0.17658538 -1.34196978 
#> 
#> $ml_classification
#> $ml_classification$accuracy
#> [1] 0.5625
#> 
#> $ml_classification$precision
#> [1] 0.5727273
#> 
#> $ml_classification$recall
#> [1] 0.5625
#> 
#> $ml_classification$F1
#> [1] 0.5465587
#> 
#> 
#> $pvalue_uniformity
#> $pvalue_uniformity$statistic
#> [1] 8
#> 
#> $pvalue_uniformity$pvalue
#> [1] 0.5341462
#> 
#> $pvalue_uniformity$distribution_score
#> [1] 1
#> 
#> 
```
