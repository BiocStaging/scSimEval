# Automatically Infer Pseudotime Trajectory from scRNA-seq Counts

Estimates cell differentiation pseudotime directly from single-cell
expression counts using diffusion/principal curve projection along the
first principal component of top highly variable genes.

## Usage

``` r
infer_scrna_pseudotime(data, n_top = 500)
```

## Arguments

- data:

  Expression or count matrix (genes x cells).

- n_top:

  Number of top variable genes to use (default 500).

## Value

Numeric vector of pseudotime values in \[0, 1\] for each cell.

## Examples

``` r
data(example_scrna, package = "scSimEval")
infer_scrna_pseudotime(example_scrna$ref)
#>   Cell_01   Cell_02   Cell_03   Cell_04   Cell_05   Cell_06   Cell_07   Cell_08 
#> 0.4077206 0.3824215 0.2083128 0.3588559 0.3305395 0.2891701 0.1925207 0.4602022 
#>   Cell_09   Cell_10   Cell_11   Cell_12   Cell_13   Cell_14   Cell_15   Cell_16 
#> 0.2697856 0.2915467 0.2126668 0.5736764 0.5735789 0.2865348 0.2232895 0.2670922 
#>   Cell_17   Cell_18   Cell_19   Cell_20   Cell_21   Cell_22   Cell_23   Cell_24 
#> 0.3264004 0.3412905 0.3647760 0.5085966 0.2759201 0.4556877 0.2732168 0.4202623 
#>   Cell_25   Cell_26   Cell_27   Cell_28   Cell_29   Cell_30   Cell_31   Cell_32 
#> 0.8008962 0.4600573 0.3223783 0.2694681 0.5034364 0.2783259 0.2966012 0.2199424 
#>   Cell_33   Cell_34   Cell_35   Cell_36   Cell_37   Cell_38   Cell_39   Cell_40 
#> 0.7308869 0.3190394 0.3035495 0.2611547 0.3065771 0.6495041 0.2249674 0.3793338 
#>   Cell_41   Cell_42   Cell_43   Cell_44   Cell_45   Cell_46   Cell_47   Cell_48 
#> 0.5402216 0.2549464 0.4271279 0.4228344 0.2297063 0.7789738 0.3981666 0.3743102 
#>   Cell_49   Cell_50   Cell_51   Cell_52   Cell_53   Cell_54   Cell_55   Cell_56 
#> 0.3095495 0.3068042 0.3956266 0.5233562 0.2976209 0.2968928 0.0000000 0.4711851 
#>   Cell_57   Cell_58   Cell_59   Cell_60   Cell_61   Cell_62   Cell_63   Cell_64 
#> 0.4100254 0.2786418 0.3588908 0.4419328 0.4710782 0.4995378 0.1989761 0.5595960 
#>   Cell_65   Cell_66   Cell_67   Cell_68   Cell_69   Cell_70   Cell_71   Cell_72 
#> 0.3185634 0.1639920 1.0000000 0.5901425 0.5127609 0.4472541 0.6446821 0.3932446 
#>   Cell_73   Cell_74   Cell_75   Cell_76   Cell_77   Cell_78   Cell_79   Cell_80 
#> 0.4083712 0.4369437 0.2657939 0.2581369 0.7184793 0.5168538 0.3310625 0.5513640 
```
