# Benchmark Execution Runtime and Memory Usage

Times an expression and measures change in memory allocation, computing
both raw runtime/memory metrics and core-normalized Category VIII
scalability indicators (`core_seconds` and `memory_per_core_mib`).

## Usage

``` r
benchmark_resource_usage(expr, sim_cores = NULL)
```

## Arguments

- expr:

  Expression or function to benchmark.

- sim_cores:

  Optional positive integer specifying CPU cores used during execution
  (default `1L`).

## Value

A named list containing the evaluated expression result and resource
metrics:

- `result`:

  Return value of `expr` (or `NULL` if evaluation errored).

- `core_seconds`:

  Total CPU core-seconds (`elapsed_seconds * cores`).

- `core_hours`:

  Total CPU core-hours (`core_seconds / 3600`).

- `memory_per_core_mib`:

  Peak R heap memory per allocated core in MiB.

- `memory_per_core_gib`:

  Peak R heap memory per allocated core in GiB.

- `memory_per_core_mb`:

  Peak R heap memory per allocated core in decimal MB.

- `memory_per_core_gb`:

  Peak R heap memory per allocated core in decimal GB.

- `elapsed_seconds`:

  Wall-clock elapsed runtime in seconds.

- `elapsed_minutes`:

  Wall-clock elapsed runtime in minutes.

- `memory_mb`:

  Peak R heap memory allocated in MiB.

- `memory_mib`:

  Alias of `memory_mb` in MiB.

- `memory_mb_si`:

  Peak R heap memory allocated in decimal MB.

- `sim_cores`:

  Normalized integer core count (included only when `sim_cores` is
  supplied).

## References

Lannelongue, L., et al. (2021). Green Algorithms. *Advanced Science*,
8(12), 2100707. Holtgrewe, M., et al. (2015). *Bioinformatics*, 31(12),
1904–1912.

## Examples

``` r
benchmark_resource_usage({ Sys.sleep(0.01); 1 + 1 })
#> $result
#> [1] 2
#> 
#> $core_seconds
#> [1] 0.011
#> 
#> $core_hours
#> [1] 3.055556e-06
#> 
#> $memory_per_core_mib
#> [1] 0.004438782
#> 
#> $memory_per_core_gib
#> [1] 4.334748e-06
#> 
#> $memory_per_core_mb
#> [1] 0.0046544
#> 
#> $memory_per_core_gb
#> [1] 4.6544e-06
#> 
#> $elapsed_seconds
#> [1] 0.011
#> 
#> $elapsed_minutes
#> [1] 0.0001833333
#> 
#> $memory_mb
#> [1] 0.004438782
#> 
#> $memory_mib
#> [1] 0.004438782
#> 
#> $memory_mb_si
#> [1] 0.0046544
#> 
```
