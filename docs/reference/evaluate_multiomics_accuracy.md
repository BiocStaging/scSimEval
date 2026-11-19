# Master Multiomics Simulation Benchmarking Pipeline

Orchestrates comprehensive benchmarking across all 8 evaluation
categories using `SingleCellExperiment` inputs (either a multi-modal SCE
containing an `altExp` for chromatin accessibility, or a named list of
two SCE objects such as `list(RNA = rna_sce, ATAC = atac_sce)`).
Cell-type and batch labels are read directly from `colData` via
`cell_type_col` and `batch_col`. Users starting from `Seurat` objects
should convert via
[`Seurat::as.SingleCellExperiment()`](https://satijalab.org/seurat/reference/as.SingleCellExperiment.html).

## Usage

``` r
evaluate_multiomics_accuracy(
  ref_multi,
  sim_multi,
  cell_type_col = "cell_type",
  batch_col = "batch",
  pairing = c("paired", "unpaired"),
  memory_mb = NULL,
  elapsed_time = NULL,
  time_unit = c("seconds", "minutes"),
  memory_unit = c("MiB", "MB"),
  sim_cores = NULL,
  separate_resources = NULL,
  feature_pairs = NULL,
  compute_bivariate = FALSE,
  threads = 1L,
  trajectory_num_trees = 10000L,
  verbose = TRUE
)
```

## Arguments

- ref_multi:

  Reference multiomics `SingleCellExperiment` (with `altExp`) or named
  list of two SCE objects.

- sim_multi:

  Simulated multiomics `SingleCellExperiment` (with `altExp`) or named
  list of two SCE objects.

- cell_type_col:

  Column name in `colData` containing cell-type labels (default
  `"cell_type"`).

- batch_col:

  Column name in `colData` containing batch/replicate labels (default
  `"batch"`).

- pairing:

  Character string specifying `"paired"` (default) or `"unpaired"` assay
  structure.

- memory_mb:

  Optional numeric peak memory allocation (in `memory_unit`).

- elapsed_time:

  Optional numeric wall-clock runtime (in `time_unit`).

- time_unit:

  Unit of `elapsed_time`: `"seconds"` (default) or `"minutes"`.

- memory_unit:

  Unit of `memory_mb`: `"MiB"` (default) or `"MB"`.

- sim_cores:

  Optional positive integer specifying CPU cores used by the simulator.

- separate_resources:

  Optional named list (e.g., `list(rna = list(...), atac = list(...))`)
  giving per-modality `elapsed_time`, `memory_mb`, `time_unit`,
  `memory_unit`, and `sim_cores` when RNA and ATAC layers are simulated
  separately.

- feature_pairs:

  Optional 2-column data.frame of linked feature pairs.

- compute_bivariate:

  Logical, whether to compute 2D bivariate tests. Default `FALSE`.

- threads:

  Integer number of CPU threads. Default is `1L`.

- trajectory_num_trees:

  Number of Random Forest trees for Category VII `wcor_features`
  (default `10000L`).

- verbose:

  Logical, whether to print progress messages. Default is `TRUE`.

## Value

A named list containing:

- benchmark_summary_table:

  Master tidy data frame with columns `Category`, `Property`, `Metric`,
  `Value`, `Modality` (and `Sim_Cores` if provided).

- metrics_summary_table:

  Alias of `benchmark_summary_table`.

- summary_table:

  Alias of `benchmark_summary_table`.

- pairing:

  Resolved pairing mode (`"paired"` or `"unpaired"`).

- mod1_unimodal:

  Full unimodal evaluation output for modality 1.

- mod2_unimodal:

  Full unimodal evaluation output for modality 2.

- clustering:

  List with elements `m1` and `m2` containing Category III results.

- batch:

  Named list of Category IV batch mixing results per modality.

- signal:

  Named list of Category V differential fidelity results per modality.

- trajectory:

  Category VII *dynverse* trajectory results inferred from the scRNA-seq
  layer.

- cross_modality:

  List storing Category VI integration output under `correlation`.

- resource_usage:

  Normalized Category VIII scalability list supporting both `"combined"`
  and `"separate"` modes.

## Details

![scSimEval Unified Benchmarking Workflow
Diagram](figures/workflow_diagram.png)

**Multiomics Metric Measurement Capacity (56 Evaluation Measures
Suite):**

- **Paired Multiomics (56 Metrics):** Dual unimodal properties (Cat
  I–II), dual clustering (Cat III), dual batch mixing (Cat IV), dual
  biological signal (Cat V), cross-modal integration & pairing (Cat VI),
  scRNA trajectory (Cat VII), and scalability (Cat VIII).

- **Unpaired Multiomics (51 Metrics):** Evaluates both modalities
  independently with population-level integration & alignment, omitting
  cell-matched pairing metrics.

## Examples

``` r
data(example_multiomics, package = "scSimEval")
res <- evaluate_multiomics_accuracy(
  example_multiomics$ref_multi,
  example_multiomics$sim_multi,
  verbose = FALSE
)
#> calcNormFactors has been renamed to normLibSizes
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> calcNormFactors has been renamed to normLibSizes
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> calcNormFactors has been renamed to normLibSizes
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> calcNormFactors has been renamed to normLibSizes
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> Warning: Zero sample variances detected, have been offset away from zero
#> calcNormFactors has been renamed to normLibSizes
#> Warning: Zero sample variances detected, have been offset away from zero
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> Warning: procv: no points with non-zero weight
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> Warning: Zero sample variances detected, have been offset away from zero
#> calcNormFactors has been renamed to normLibSizes
#> Warning: Zero sample variances detected, have been offset away from zero
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> Warning: Zero sample variances detected, have been offset away from zero
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> Warning: Zero sample variances detected, have been offset away from zero
#> Warning: shoving 'interior' knots matching boundary knots to inside
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> Warning: Zero sample variances detected, have been offset away from zero
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
```
