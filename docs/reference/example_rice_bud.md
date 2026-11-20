# Real Paired Single-Cell Multiomics Benchmark Dataset (Rice Bud, SingleCellExperiment)

A real, pre-processed paired single-cell RNA-seq and ATAC-seq benchmark
dataset of rice bud tissue (*Oryza sativa*, GEO accession GSE232863),
standardized using Bioconductor `SingleCellExperiment` objects. Contains
paired empirical reference profiles (top 100 HVGs x 100 balanced cells
for RNA, top 100 HVPs x 100 identical paired cells for ATAC) covering
all 8 biological cell types equally (12–13 cells per type) and 2
biological batches (`Bud_21`, `Bud_22`), alongside simulated outputs
from 6 single-cell multiomics simulators (`scMultiSim`, `scDesign3`,
`MOSim`, `scCross`, `scMoMtF`, and `Matilda`) and observed scalability
benchmarks (runtime in seconds, peak RAM in MiB, and CPU cores).

## Usage

``` r
data(example_rice_bud)
```

## Format

A list with 4 elements:

- metadata:

  List containing dataset name, organism, tissue, pairing mode
  (`"paired"`), dimensions, and simulators.

- reference:

  List containing standardized `SingleCellExperiment` objects:
  `combined` (multi-modal `SingleCellExperiment` with RNA counts and
  `altExp(combined, "ATAC")` sharing identical paired cells and
  `colData`), `rna` (unimodal `SingleCellExperiment` for scRNA-seq), and
  `atac` (unimodal `SingleCellExperiment` for scATAC-seq).

- simulated:

  Named list of simulated outputs for 6 simulators: `scMultiSim`,
  `scDesign3`, `MOSim`, `scCross`, `scMoMtF`, and `Matilda`. Each
  simulator entry contains `combined`, `rna`, and `atac`
  `SingleCellExperiment` objects matching the reference structure.

- scalability_benchmarks:

  A data.frame with benchmarked `Method`, `Runtime_Seconds`,
  `PeakRAM_MiB`, and `Sim_Cores`.

## Source

Gene Expression Omnibus (GEO) GSE232863. Preprocessed via standard rice
tissue pipeline.

## Examples

``` r
data(example_rice_bud)
example_rice_bud$reference$combined
#> class: SingleCellExperiment 
#> dim: 100 100 
#> metadata(0):
#> assays(1): counts
#> rownames(100): Os04g0188433 Os05g0202800 ... Os06g0553200 Os04g0578600
#> rowData names(0):
#> colnames(100): AAACCGGCAAGTGTTT-1 AAACGGATCGATCAGT-1 ...
#>   TTTGGTGCATGTCAGC-1 TTTGTGTTCATAACGC-1
#> colData names(7): orig.ident nCount_RNA ... nCount_ATAC nFeature_ATAC
#> reducedDimNames(1): PCA
#> mainExpName: NULL
#> altExpNames(1): ATAC
dim(example_rice_bud$reference$rna)
#> [1] 100 100
dim(example_rice_bud$reference$atac)
#> [1] 100 100
names(example_rice_bud$simulated)
#> [1] "scMultiSim" "scDesign3"  "MOSim"      "scCross"    "scMoMtF"   
#> [6] "Matilda"   
example_rice_bud$scalability_benchmarks
#>       Method Runtime_Seconds PeakRAM_MiB Sim_Cores
#> 1 scMultiSim          191.89   1139.0000        16
#> 2  scDesign3           87.87   5713.9000        16
#> 3      MOSim          294.45   1399.5000        16
#> 4    scCross          867.68   1261.5585        16
#> 5    scMoMtF           13.84    323.9536        16
#> 6    Matilda           57.32    335.9413        16
# \donttest{
res <- evaluate_multiomics_accuracy(
  ref_multi = list(RNA = example_rice_bud$reference$rna,
                   ATAC = example_rice_bud$reference$atac),
  sim_multi = list(RNA = example_rice_bud$simulated$scDesign3$rna,
                   ATAC = example_rice_bud$simulated$scDesign3$atac),
  pairing = "paired",
  verbose = FALSE
)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Warning: Zero sample variances detected, have been offset away from zero
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
#> 1 22.87000000      RNA
#> 2  0.20000000      RNA
#> 3 10.00000000      RNA
#> 4  0.00564741      RNA
#> 5 -2.83556062      RNA
#> 6  0.06500000      RNA
#> 7  0.02451093      RNA
#> 8 22.87000000      RNA
# }
```
