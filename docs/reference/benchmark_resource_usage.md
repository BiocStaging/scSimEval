# Benchmark Execution Runtime and Memory Usage

Times an expression and measures change in memory allocation.

## Usage

``` r
benchmark_resource_usage(expr)
```

## Arguments

- expr:

  Expression or function to benchmark.

## Value

A list with elapsed seconds and peak memory allocated (MB).

## Examples

``` r
benchmark_resource_usage({ Sys.sleep(0.01); 1 + 1 })
#> $result
#> [1] 2
#> 
#> $elapsed_seconds
#> [1] 0.02
#> 
#> $memory_mb
#> [1] 0.003444672
#> 
```
