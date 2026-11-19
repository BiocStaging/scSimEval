#' Example scRNA-seq Benchmarking Dataset (Rice Bud)
#'
#' A real, subsampled single-cell RNA sequencing dataset from rice bud tissue
#' (\emph{Oryza sativa}) designed for testing, benchmarking, and demonstrating
#' evaluation metrics in \code{scSimEval}. Contains empirical reference counts
#' and \code{scDesign3}-simulated counts across the top 100 highly variable genes
#' (HVGs) and 100 balanced cells covering all 8 biological cell types and 2 batches.
#'
#' @docType data
#' @name example_scrna
#' @usage data(example_scrna)
#' @format A list with 4 elements:
#' \describe{
#'   \item{ref}{A \code{SingleCellExperiment} object containing reference scRNA-seq counts (100 HVGs x 100 cells), cell types and batch annotations in \code{colData}, and precomputed \code{PCA} in \code{reducedDims}.}
#'   \item{sim}{A \code{SingleCellExperiment} object containing simulated scRNA-seq counts (100 HVGs x 100 cells), cell types and batch annotations in \code{colData}, and precomputed \code{PCA} in \code{reducedDims}.}
#'   \item{cell_types}{A named factor of cell-type identities across 8 Rice Bud cell types for the 100 cells.}
#'   \item{batch_info}{A named factor of batch annotations (\code{Bud_21}, \code{Bud_22}) for the 100 cells.}
#' }
#' @source Derived from rice bud paired multiomics data (\emph{Oryza sativa}, GEO accession GSE232863).
#' @examples
#' data(example_scrna)
#' dim(example_scrna$ref)
#' dim(example_scrna$sim)
#' table(example_scrna$cell_types)
"example_scrna"

#' Example scATAC-seq Benchmarking Dataset (Rice Bud)
#'
#' A real, subsampled single-cell ATAC sequencing peak accessibility dataset from
#' rice bud tissue (\emph{Oryza sativa}) designed for testing, benchmarking, and
#' demonstrating evaluation metrics in \code{scSimEval}. Contains empirical reference
#' peak counts and \code{scDesign3}-simulated counts across the top 100 highly variable
#' peaks (HVPs) and the same 100 balanced paired cells as \code{example_scrna}.
#'
#' @docType data
#' @name example_scatac
#' @usage data(example_scatac)
#' @format A list with 4 elements:
#' \describe{
#'   \item{ref}{A \code{SingleCellExperiment} object containing reference peak accessibility counts (100 HVPs x 100 cells), cell types and batch annotations in \code{colData}, and precomputed \code{PCA} in \code{reducedDims}.}
#'   \item{sim}{A \code{SingleCellExperiment} object containing simulated peak accessibility counts (100 HVPs x 100 cells), cell types and batch annotations in \code{colData}, and precomputed \code{PCA} in \code{reducedDims}.}
#'   \item{cell_types}{A named factor of cell-type identities across 8 Rice Bud cell types for the 100 cells.}
#'   \item{batch_info}{A named factor of batch annotations (\code{Bud_21}, \code{Bud_22}) for the 100 cells.}
#' }
#' @source Derived from rice bud paired multiomics data (\emph{Oryza sativa}, GEO accession GSE232863).
#' @examples
#' data(example_scatac)
#' dim(example_scatac$ref)
#' dim(example_scatac$sim)
#' table(example_scatac$batch_info)
"example_scatac"

#' Example Paired Multiomics Benchmarking Dataset (Rice Bud)
#'
#' A unified paired single-cell multiomics dataset combining \code{example_scrna}
#' (100 HVGs) and \code{example_scatac} (100 HVPs) across the exact same 100 balanced
#' rice bud cells, accompanied by shared cell-type labels, batch annotations, and
#' empirical computational scalability metrics.
#'
#' @docType data
#' @name example_multiomics
#' @usage data(example_multiomics)
#' @format A list with 5 elements:
#' \describe{
#'   \item{ref_multi}{Named list of \code{SingleCellExperiment} objects: \code{rna} (100 HVGs x 100 cells) and \code{atac} (100 HVPs x 100 cells), complete with \code{colData} and precomputed \code{PCA} in \code{reducedDims}.}
#'   \item{sim_multi}{Named list of simulated \code{SingleCellExperiment} objects: \code{rna} (100 x 100) and \code{atac} (100 x 100).}
#'   \item{cell_types}{A named factor of shared cell-type identities across 8 Rice Bud cell types for the 100 paired cells.}
#'   \item{batch_info}{A named factor of batch annotations (\code{Bud_21}, \code{Bud_22}) for the 100 paired cells.}
#'   \item{resource_stats}{A list of empirical scalability metrics: \code{elapsed_time}, \code{memory_mb}, \code{sim_cores}.}
#' }
#' @source Derived from rice bud paired multiomics data (\emph{Oryza sativa}, GEO accession GSE232863).
#' @examples
#' data(example_multiomics)
#' names(example_multiomics$ref_multi)
#' names(example_multiomics$sim_multi)
#' @export
"example_multiomics"

#' Example Cell-Type Annotations Vector (Rice Bud)
#'
#' A named factor vector containing discrete cell-type labels across all 8 Rice Bud
#' cell types for 100 balanced single cells (12--13 cells per type), matching the
#' column names of \code{\link{example_scrna}} and \code{\link{example_scatac}}.
#'
#' @docType data
#' @name example_cell_types
#' @usage data(example_cell_types)
#' @format A named factor with 100 elements across 8 levels (\code{Bud-unknown},
#'   \code{Bud_axis}, \code{Bud_primordium}, \code{Epidermis}, \code{Parenchyma_cell},
#'   \code{Root_meristem}, \code{Shoot_apical_meristem}, \code{Vascular_cylinder}).
#' @source Derived from \code{example_scrna$cell_types}.
#' @examples
#' data(example_cell_types)
#' head(example_cell_types)
#' table(example_cell_types)
#' @export
"example_cell_types"

#' Example Batch Annotations Vector (Rice Bud)
#'
#' A named factor vector containing biological batch annotations (\code{Bud_21},
#' \code{Bud_22}) for 100 single cells, matching the column names of
#' \code{\link{example_scrna}} and \code{\link{example_scatac}}.
#'
#' @docType data
#' @name example_batch_info
#' @usage data(example_batch_info)
#' @format A named factor with 100 elements across 2 levels: \code{Bud_21} and \code{Bud_22}.
#' @source Derived from \code{example_scrna$batch_info}.
#' @examples
#' data(example_batch_info)
#' head(example_batch_info)
#' table(example_batch_info)
#' @export
"example_batch_info"

#' Real Paired Single-Cell Multiomics Benchmark Dataset (Rice Bud, Seurat v5)
#'
#' A real, pre-processed paired single-cell RNA-seq and ATAC-seq benchmark dataset of
#' rice bud tissue (\emph{Oryza sativa}, GEO accession GSE232863), standardized
#' using \code{Seurat} v5 objects.
#' Contains paired empirical reference profiles (top 100 HVGs x 100 balanced cells for RNA,
#' top 100 HVPs x 100 identical paired cells for ATAC) covering all 8 biological cell types
#' equally (12--13 cells per type) and 2 biological batches (\code{Bud_21}, \code{Bud_22}),
#' alongside simulated outputs from 6 single-cell multiomics simulators (\code{scMultiSim},
#' \code{scDesign3}, \code{MOSim}, \code{scCross}, \code{scMoMtF}, and \code{Matilda})
#' and observed scalability benchmarks (runtime in seconds, peak RAM in MiB, and CPU cores).
#'
#' @docType data
#' @name example_rice_bud
#' @usage data(example_rice_bud)
#' @format A list with 4 elements:
#' \describe{
#'   \item{metadata}{List containing dataset name, organism, tissue, pairing mode (\code{"paired"}), dimensions, and simulators.}
#'   \item{reference}{List containing standardized \code{Seurat} v5 objects:
#'     \code{combined} (multi-modal \code{Seurat} with \code{RNA} 100 x 100 and \code{ATAC} 100 x 100 sharing identical paired cells and cell types),
#'     \code{rna} (unimodal \code{Seurat} object for scRNA-seq, Data Mode 1), and
#'     \code{atac} (unimodal \code{Seurat} object for scATAC-seq, Data Mode 2).}
#'   \item{simulated}{Named list of simulated outputs for 6 simulators: \code{scMultiSim}, \code{scDesign3}, \code{MOSim}, \code{scCross}, \code{scMoMtF}, and \code{Matilda}. Each simulator entry contains \code{combined}, \code{rna}, and \code{atac} \code{Seurat} objects matching the reference structure.}
#'   \item{scalability_benchmarks}{A data.frame with benchmarked \code{Method}, \code{Runtime_Seconds}, \code{PeakRAM_MiB}, and \code{Sim_Cores}.}
#' }
#' @source Gene Expression Omnibus (GEO) GSE232863. Preprocessed via standard rice tissue pipeline.
#' @examples
#' data(example_rice_bud)
#' example_rice_bud$reference$combined
#' dim(example_rice_bud$reference$rna)
#' dim(example_rice_bud$reference$atac)
#' names(example_rice_bud$simulated)
#' \donttest{
#' evaluate_multiomics_accuracy(
#'   ref_multi = list(RNA = example_rice_bud$reference$rna, ATAC = example_rice_bud$reference$atac),
#'   sim_multi = list(RNA = example_rice_bud$simulated$scDesign3$rna, ATAC = example_rice_bud$simulated$scDesign3$atac),
#'   is_paired = TRUE
#' )
#' }
"example_rice_bud"

#' Real Unpaired Single-Cell Multiomics Benchmark Dataset (Soybean Root, Seurat v5)
#'
#' A real, pre-processed unpaired single-cell RNA-seq and ATAC-seq benchmark dataset of
#' soybean root tissue (\emph{Glycine max}), standardized using \code{Seurat} v5 objects.
#' Contains independent empirical reference profiles (top 100 HVGs x 100 balanced RNA cells
#' and top 100 HVPs x 100 independently sampled balanced ATAC cells) covering all 8 biological
#' cell types (\code{Cortex}, \code{Endodermis}, \code{Epidermis}, \code{Pericycle},
#' \code{Phloem}, \code{Procambium}, \code{Root_hair}, \code{Xylem}; 12--13 cells per type),
#' alongside simulated outputs from 3 unpaired multiomics simulators (\code{scMultiSim},
#' \code{scDesign3}, and \code{scCross}) and observed scalability benchmarks.
#'
#' @docType data
#' @name example_soybean_root
#' @usage data(example_soybean_root)
#' @format A list with 4 elements:
#' \describe{
#'   \item{metadata}{List containing dataset name, organism, tissue, pairing mode (\code{"unpaired"}), dimensions, and simulators.}
#'   \item{reference}{List containing standardized \code{Seurat} v5 objects:
#'     \code{combined} (multi-modal \code{Seurat} spanning 200 independent cells with \code{RNA} 100 x 100 and \code{ATAC} 100 x 100),
#'     \code{rna} (unimodal \code{Seurat} object for scRNA-seq), and
#'     \code{atac} (unimodal \code{Seurat} object for scATAC-seq).}
#'   \item{simulated}{Named list of simulated outputs for 3 simulators: \code{scMultiSim}, \code{scDesign3}, and \code{scCross}. Each simulator entry contains \code{combined}, \code{rna}, and \code{atac} \code{Seurat} objects matching the reference structure.}
#'   \item{scalability_benchmarks}{A data.frame with benchmarked runtime, peak RAM, and CPU cores for RNA, ATAC, and combined execution.}
#' }
#' @source Preprocessed from soybean root single-cell atlas datasets.
#' @examples
#' data(example_soybean_root)
#' example_soybean_root$reference$combined
#' dim(example_soybean_root$reference$rna)
#' dim(example_soybean_root$reference$atac)
#' names(example_soybean_root$simulated)
#' \donttest{
#' evaluate_multiomics_accuracy(
#'   ref_multi = list(RNA = example_soybean_root$reference$rna, ATAC = example_soybean_root$reference$atac),
#'   sim_multi = list(RNA = example_soybean_root$simulated$scMultiSim$rna, ATAC = example_soybean_root$simulated$scMultiSim$atac),
#'   is_paired = FALSE
#' )
#' }
"example_soybean_root"
