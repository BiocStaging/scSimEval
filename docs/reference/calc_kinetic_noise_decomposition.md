# Kinetic Noise Decomposition (SymSim & Elowitz et al.)

Decomposes gene expression variance / squared coefficient of variation
(CV^2) into intrinsic transcriptional bursting noise and extrinsic
cell-state noise.

## Usage

``` r
calc_kinetic_noise_decomposition(counts, cell_states = NULL)
```

## Arguments

- counts:

  Matrix of expression counts (genes x cells) or SingleCellExperiment.

- cell_states:

  Optional factor or vector of cell states / subpopulation clusters. If
  provided, noise is partitioned into within-state (intrinsic) and
  between-state (extrinsic) components. If NULL, intrinsic noise is
  estimated via Poisson shot-noise expectation (1 / mean).

## Value

A list containing mean intrinsic noise, mean extrinsic noise, noise
ratio, and gene-level vectors.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_kinetic_noise_decomposition(example_scrna$ref)
#> $mean_intrinsic_noise
#> [1] 0.2613753
#> 
#> $mean_extrinsic_noise
#> [1] 0.8192499
#> 
#> $mean_total_cv2
#> [1] 1.080625
#> 
#> $mean_intrinsic_fraction
#> [1] 0.2483567
#> 
#> $gene_intrinsic_noise
#>   Gene_01   Gene_02   Gene_03   Gene_04   Gene_05   Gene_06   Gene_07   Gene_08 
#> 0.1462523 0.2234637 0.1793722 0.2040816 0.1396161 0.2150538 0.1789709 0.2099738 
#>   Gene_09   Gene_10   Gene_11   Gene_12   Gene_13   Gene_14   Gene_15   Gene_16 
#> 0.1326700 0.1773836 0.1596806 0.1523810 0.3065134 0.1805869 0.2025316 0.2402402 
#>   Gene_17   Gene_18   Gene_19   Gene_20   Gene_21   Gene_22   Gene_23   Gene_24 
#> 0.2721088 0.4705882 0.4301075 0.1895735 0.3571429 0.3846154 0.2657807 0.1946472 
#>   Gene_25   Gene_26   Gene_27   Gene_28   Gene_29   Gene_30   Gene_31   Gene_32 
#> 0.1789709 0.2409639 0.2515723 0.4166667 0.2346041 0.2797203 0.2388060 0.2197802 
#>   Gene_33   Gene_34   Gene_35   Gene_36   Gene_37   Gene_38   Gene_39   Gene_40 
#> 0.2339181 0.3448276 0.2298851 0.4232804 0.3088803 0.2622951 0.5031447 0.2185792 
#>   Gene_41   Gene_42   Gene_43   Gene_44   Gene_45   Gene_46   Gene_47   Gene_48 
#> 0.2325581 0.3252033 0.2332362 0.3686636 0.3669725 0.2614379 0.3463203 0.2067183 
#>   Gene_49   Gene_50   Gene_51   Gene_52   Gene_53   Gene_54   Gene_55   Gene_56 
#> 0.2424242 0.1826484 0.2105263 0.2919708 0.1731602 0.1970443 0.2531646 0.2298851 
#>   Gene_57   Gene_58   Gene_59   Gene_60 
#> 0.1975309 0.1985112 0.6722689 0.2930403 
#> 
#> $gene_extrinsic_noise
#>  [1] 0.9944127 0.6690320 1.6270973 0.9395347 0.9097513 0.5671834 0.6033002
#>  [8] 0.5247425 0.9844321 0.7038289 0.6902993 1.1296214 2.2545140 0.7547838
#> [15] 0.8207059 0.4369200 0.6009920 0.9331173 1.2220627 0.6993030 1.1392405
#> [22] 0.8460040 0.5447728 0.8611620 1.0176706 0.6661510 1.0317850 0.9177215
#> [29] 0.5091115 0.7280849 0.9524266 0.8011935 0.6401450 0.9784915 1.0171736
#> [36] 0.6755024 0.5636564 1.0006293 0.8010419 0.5806554 0.6689213 0.5657673
#> [43] 0.4583813 0.9704865 0.8773502 0.9909637 0.4643794 0.6148743 0.9022794
#> [50] 0.8933102 0.6697290 0.8470471 0.6633489 0.6117012 0.9863317 0.7402279
#> [57] 0.7017146 0.7227510 0.6434502 0.8237249
#> 
```
