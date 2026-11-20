# Real Unpaired Single-Cell Multiomics Benchmark Dataset (Soybean Root, SingleCellExperiment)

A real, pre-processed unpaired single-cell RNA-seq and ATAC-seq
benchmark dataset of soybean root tissue (*Glycine max*), standardized
using Bioconductor `SingleCellExperiment` objects. Contains independent
empirical reference profiles (top 100 HVGs x 100 balanced RNA cells and
top 100 HVPs x 100 independently sampled balanced ATAC cells) covering
all 8 biological cell types (`Cortex`, `Endodermis`, `Epidermis`,
`Pericycle`, `Phloem`, `Procambium`, `Root_hair`, `Xylem`; 12–13 cells
per type), alongside simulated outputs from 3 unpaired multiomics
simulators (`scMultiSim`, `scDesign3`, and `scCross`) and observed
scalability benchmarks.

## Usage

``` r
data(example_soybean_root)
```

## Format

A list with 4 elements:

- metadata:

  List containing dataset name, organism, tissue, pairing mode
  (`"unpaired"`), dimensions, and simulators.

- reference:

  List containing standardized `SingleCellExperiment` objects:
  `combined` (multi-modal `SingleCellExperiment` with `colData$modality`
  distinguishing independent RNA and ATAC cells across
  `altExp(combined, "ATAC")`), `rna` (unimodal `SingleCellExperiment`
  for scRNA-seq), and `atac` (unimodal `SingleCellExperiment` for
  scATAC-seq).

- simulated:

  Named list of simulated outputs for 3 simulators: `scMultiSim`,
  `scDesign3`, and `scCross`. Each simulator entry contains `combined`,
  `rna`, and `atac` `SingleCellExperiment` objects matching the
  reference structure.

- scalability_benchmarks:

  A data.frame with benchmarked runtime, peak RAM, and CPU cores for
  RNA, ATAC, and combined execution.

## Source

Preprocessed from soybean root single-cell atlas datasets.

## Examples

``` r
data(example_soybean_root)
example_soybean_root$reference$combined
#> class: SingleCellExperiment 
#> dim: 100 200 
#> metadata(0):
#> assays(1): counts
#> rownames(100): Glyma.06G242900 Glyma.12G154300 ... Glyma.13G272000
#>   Glyma.04G048100
#> rowData names(0):
#> colnames(200): RNA_1 RNA_2 ... ATAC_99 ATAC_100
#> colData names(8): orig.ident nCount_RNA ... nCount_ATAC nFeature_ATAC
#> reducedDimNames(1): PCA
#> mainExpName: NULL
#> altExpNames(1): ATAC
dim(example_soybean_root$reference$rna)
#> [1] 100 100
dim(example_soybean_root$reference$atac)
#> [1] 100 100
names(example_soybean_root$simulated)
#> [1] "scMultiSim" "scDesign3"  "scCross"   
example_soybean_root$scalability_benchmarks
#>       Method Runtime_Seconds PeakRAM_MiB Sim_Cores RNA_Runtime_Seconds
#> 1 scMultiSim          107.55    1354.500        16              107.55
#> 2  scDesign3          504.30   14073.900         4               69.20
#> 3    scCross         1791.52    2053.928        16             1791.52
#>   RNA_PeakRAM_MiB RNA_Sim_Cores ATAC_Runtime_Seconds ATAC_PeakRAM_MiB
#> 1        1354.500            16               107.55         1354.500
#> 2        2736.100            16               435.10        14073.900
#> 3        2053.928            16              1791.52         2053.928
#>   ATAC_Sim_Cores
#> 1             16
#> 2              4
#> 3             16
# \donttest{
res <- evaluate_multiomics_accuracy(
  ref_multi = list(RNA = example_soybean_root$reference$rna,
                   ATAC = example_soybean_root$reference$atac),
  sim_multi = list(RNA = example_soybean_root$simulated$scMultiSim$rna,
                   ATAC = example_soybean_root$simulated$scMultiSim$atac),
  pairing = "unpaired",
  verbose = FALSE
)
head(res$benchmark_summary_table, 8)
#>                    Category     Property                      Metric
#> 1 Distributional Properties library_size                 Wasserstein
#> 2 Distributional Properties library_size                          KS
#> 3 Distributional Properties library_size                         MAD
#> 4 Distributional Properties library_size               Bhattacharyya
#> 5 Distributional Properties library_size              Runs_Statistic
#> 6 Distributional Properties library_size                 NN_Mismatch
#> 7 Distributional Properties library_size Between_Dataset_Silh_Global
#> 8 Distributional Properties library_size                         MAE
#>         Value Modality
#> 1 37.61000000      RNA
#> 2  0.65000000      RNA
#> 3 38.00000000      RNA
#> 4  0.07370043      RNA
#> 5 -7.37245762      RNA
#> 6  0.50500000      RNA
#> 7  0.27304482      RNA
#> 8 37.61000000      RNA
# }
```
