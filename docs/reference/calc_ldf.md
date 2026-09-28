# Local Density Factor (LDF)

Computes the Local Density Estimate (LDE) and Local Density Factor (LDF)
using a Gaussian kernel over reachability distances in the k-NN
neighborhood, adapted from Latecki et al. and CellMixS.

## Usage

``` r
calc_ldf(coords, k = 15, h = 1, c = 1)
```

## Arguments

- coords:

  Matrix of coordinates / embeddings (cells x dimensions).

- k:

  Number of nearest neighbors. Default is 15.

- h:

  Bandwidth parameter for Gaussian kernel. Default is 1.

- c:

  Scaling constant for comparison of LDE to neighboring observations.
  Default is 1.

## Value

A named list:

- lde:

  Local density estimate for each cell

- ldf:

  Local density factor for each cell

## Examples

``` r
coords <- matrix(stats::rnorm(100), 50, 2)
batch <- factor(rep(c("B1", "B2"), length.out = 50))
calc_ldf(coords)
#> $lde
#>  [1] 0.04392640 0.06863963 0.07059773 0.05947652 0.05519008 0.05542242
#>  [7] 0.07058913 0.02104998 0.04236533 0.05086161 0.06004533 0.06706947
#> [13] 0.02617925 0.05902243 0.06475610 0.07136396 0.06796081 0.07091042
#> [19] 0.06591993 0.06829885 0.07105228 0.07011438 0.05890794 0.07273537
#> [25] 0.06861190 0.06566178 0.06952166 0.06971279 0.07029205 0.06471481
#> [31] 0.05678199 0.06018493 0.05093880 0.04956676 0.06798634 0.07399568
#> [37] 0.06650839 0.07151670 0.06759566 0.07084081 0.06184846 0.06095425
#> [43] 0.06768948 0.06682576 0.06786471 0.07292855 0.07268008 0.05611276
#> [49] 0.07272839 0.06724810
#> 
#> $ldf
#>  [1] 0.5808307 0.4986130 0.4914289 0.5229573 0.5421240 0.5304372 0.4888909
#>  [8] 0.7417946 0.5880810 0.5634012 0.5209854 0.4969409 0.7121219 0.5113306
#> [15] 0.5106579 0.4809139 0.4954857 0.4923190 0.4975911 0.4926206 0.4899978
#> [22] 0.4888430 0.5154732 0.4889349 0.4955440 0.4953893 0.4864983 0.4937084
#> [29] 0.4943349 0.5108274 0.5422449 0.5226592 0.5398283 0.5736045 0.4978009
#> [36] 0.4814279 0.5028557 0.4920200 0.4956763 0.4905105 0.5228449 0.5063961
#> [43] 0.5023297 0.4979123 0.5009823 0.4867952 0.4877077 0.5413959 0.4904279
#> [50] 0.5034150
#> 
```
