# Interactive Shiny Studio — Complete User Guide

## Overview

The **scSimEval Shiny Studio** is a fully self-contained, browser-based
graphical interface that gives researchers access to all 62 evaluation
measures and every visualization in the package — without writing a
single line of R code.

Launch it with:

``` r
library(scSimEval)
launch_scSimEval_app()          # opens at http://127.0.0.1:<port>
```

> **File size limit:** The app accepts uploads up to **500 MB** per
> file, which handles most real-world single-cell datasets.

The studio is organized into **six primary navigation tabs** in the top
bar, plus an **eight-panel Visualizations module** inside Tab 4. The
sections below document every tab, every control, and all export options
in the order they appear in the application.

------------------------------------------------------------------------

## Tab 1 — Home

The landing page gives a quick orientation and provides **one-click jump
buttons** to any part of the workflow:

| Button                           | Destination                |
|----------------------------------|----------------------------|
| **1. Data Upload & Evaluation**  | Data Hub tab               |
| **2. Comparative Bubble Matrix** | Bubble Matrix tab          |
| **3. Diagnostic Visualizations** | Visualizations tab         |
| **4. Download Results**          | Download Results tab       |
| **5. Help & Manual**             | Help & Getting Started tab |
| **6. Team & Contact**            | Team & Contact tab         |

The four **KPI cards** below the hero section display: 62 curated
evaluation measures, 8 canonical categories, 3 supported modalities, and
the \[0, 1\] standardized fidelity scale.

------------------------------------------------------------------------

## Tab 2 — Data Hub

Data ingestion and evaluation control centre. Choose one of **four
evaluation modes** from the radio button panel on the left sidebar.

### Mode 1 — Explore Demo Benchmark (6 Simulators)

Click **Load Demo Benchmark (6 Simulators)** to instantly populate the
app with pre-computed results for **Splatter**, **scDesign3**,
**SCRIP**, **SymSim**, **dyngen**, and **simATAC** across all 62
measures. No file uploads required.

### Mode 2 — Single-Cell Evaluation (scRNA-seq / scATAC-seq)

Evaluate one or more unimodal count-matrix simulators.

| Step | Control | Notes |
|----|----|----|
| **1** | **Reference Count Matrix** | Real empirical cells. Formats: `.rds`, `.csv`, `.tsv`, `.txt`, `SingleCellExperiment`, `Seurat`. |
| **2** | **Simulated Count Matrices** (multi-select) | One file per simulator. Same formats. |
| **3** | **Simulator Name / Runtime (s) / Memory (MiB)** | Enter for each uploaded simulator. |
| **4** | **Cell Type Labels** (optional) | `.rds` (factor / vector / data.frame) or `.csv`/`.tsv`. |
| **4** | **Batch Annotations** (optional) | Same formats as cell types. |
| — | **Append to current benchmark** | Adds new results alongside existing ones. |
| — | **Evaluate Single-Cell Simulators** | Triggers the full 50+ metric pipeline. |

### Mode 3 — Multiomics Evaluation (scRNA-seq + scATAC-seq)

| Subtype | Description |
|----|----|
| **Paired Co-assay** | Same individual cells in both modalities (10x Multiome, SHARE-seq). All 62 metrics including FOSCTTM and <Match@1>. |
| **Unpaired Profiling** | Separate cells from same tissue. Population-level cross-modal metrics only; cell-pairing metrics safely skipped. |

| Step | Control | Notes |
|----|----|----|
| **1** | **Reference RNA Count Matrix** | Genes × cells. |
| **1** | **Reference ATAC Count Matrix** | Peaks × cells. |
| **2** | **Number of Multiomics Simulators** | 1 – 5. Generates matching upload fields. |
| **2** | **Per-Simulator Fields** | Simulated RNA matrix, simulated ATAC matrix, name, runtime (s), memory (MiB). |
| **3** | **Cell Type Labels / Batch** (optional) | Applied to RNA modality. |
| — | **Evaluate Multiomics Simulators** | Runs the full 62-metric multiomics pipeline. |

### Mode 4 — Upload Saved Results (.rds / .csv)

Load output previously saved from
[`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md),
[`evaluate_multiomics_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_accuracy.md),
or
[`evaluate_multiple_datasets()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiple_datasets.md):

``` r
bench <- evaluate_simulation_accuracy(ref, sim, method_name = "MyTool")
saveRDS(bench, "my_benchmark.rds")
```

Click **Load Saved File** to restore the entire benchmark session.

### Main Panel — Active Benchmark & Dataset Summary

| Section | Content |
|----|----|
| **Status Banner** | Number of simulators loaded and total metrics computed |
| **Dataset Properties Table** | Per-matrix: Cells (N), Features (P), Sparsity (%), Cell Types, Batches, Median Library Size, Median Detected Features, Mean Expression |
| **Download Dataset Summary** | Export properties table as CSV |
| **Benchmark Preview Table** | Searchable / sortable DT table of all raw and normalized metric scores |

------------------------------------------------------------------------

## Tab 3 — Comparative Bubble Matrix

The flagship visualization: every simulator × all 62 metrics in one
colour-coded matrix. Bubble **diameter** is proportional to the
normalized fidelity score. Scores ≥ 0.96 are rendered as bold square
glyphs. Simulators are rank-ordered top-to-bottom by composite average
score.

### Sidebar Controls

| Control | Options | Default |
|----|----|----|
| **Category Filter** | All 62 · (I) Distributional Properties · (II) Correlations & Zero-Inflation · (III) Cellular Structure & Concordance · (IV) Batch Effects & Confounder Mixing · (V) Biological Signal & Downstream Fidelity · (VI) Trajectory & Lineage Dynamics · (VII) Cross-Modal Coupling & Modularity · (VIII) Computational Scalability | All |
| **Method Selector** | Dynamic checkboxes per simulator | All selected |
| **Matrix Display Width (px)** | 1 200 – 3 200 (step 50) | 2 200 |
| **Matrix Display Height (px)** | 450 – 1 100 (step 20) | 680 |

> **Tip:** Default 2 200 × 680 px matches the publication format (23 : 7
> aspect ratio). Use the horizontal scrollbar to inspect all 62 measures
> without label squishing.

### Exports

| Button                      | Output                    |
|-----------------------------|---------------------------|
| **Download JPEG (600 DPI)** | Publication-grade `.jpeg` |
| **Download Vector PDF**     | Scalable `.pdf`           |

------------------------------------------------------------------------

## Tab 4 — Visualizations

Eight diagnostic sub-panels accessed via the **pill navigation bar**.
Every sub-panel exposes individual 600 DPI JPEG and vector PDF export
buttons.

### Sub-panel 1 — Evaluation Summary

[`plot_evaluation_summary()`](https://kabilanbio.github.io/scSimEval/reference/plot_evaluation_summary.md)
— rank-ordered horizontal bar chart.

| Control                       | Description                         |
|-------------------------------|-------------------------------------|
| **Show Score Labels**         | Toggle numerical scores on each bar |
| **Normalize Scores \[0, 1\]** | Standardized vs. raw display        |

### Sub-panel 2 — Distribution QC

[`plot_distribution_qc()`](https://kabilanbio.github.io/scSimEval/reference/plot_distribution_qc.md)
— 14-panel grid comparing expression densities, library sizes, and
zero-inflation curves between reference and simulated cells.

| Control       | Options                                         |
|---------------|-------------------------------------------------|
| **QC Layout** | Comprehensive (14 panels) · Density Curves Only |

### Sub-panel 3 — Scalability Benchmark

[`plot_scalability_benchmark()`](https://kabilanbio.github.io/scSimEval/reference/plot_scalability_benchmark.md)
— runtime and memory dashboards with Pareto efficiency frontiers.

| Control | Options |
|----|----|
| **Scalability View** | 4-Panel Comprehensive · Runtime Only · Peak RAM Only · Runtime vs Memory Trade-Off · Resource Cost Footprint |

### Sub-panel 4 — Metric Boxplots

[`plot_metric_boxplots()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_boxplots.md)
/
[`plot_individual_metric_bar()`](https://kabilanbio.github.io/scSimEval/reference/plot_individual_metric_bar.md)
— two view modes.

**Individual Metric mode:**

| Control                | Description                             |
|------------------------|-----------------------------------------|
| **1. Choose Category** | Filters the metric dropdown             |
| **2. Choose Metric**   | Dynamic; populated by selected category |
| **3. Score Type**      | Normalized \[0, 1\] · Raw Metric Value  |

**Category Group mode:**

| Control             | Description                              |
|---------------------|------------------------------------------|
| **Choose Category** | All 8 combined, or one specific category |
| **Score Type**      | Normalized \[0, 1\] · Raw Value          |

### Sub-panel 5 — Metric Heatmap

[`plot_metric_heatmap()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_heatmap.md)
— method × metric grid; exact raw scores in bold text with
direction-aware fill colours.

| Control                 | Options                              | Default |
|-------------------------|--------------------------------------|---------|
| **Category Filter**     | All 8 combined · individual category | All     |
| **Heatmap Height (px)** | 400 – 1 500 (step 50)                | 1 100   |
| **Heatmap Width (px)**  | 600 – 1 300 (step 20)                | 880     |

### Sub-panel 6 — PCA Ordination

[`plot_metric_pca()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_pca.md)
— simulators projected into performance space with discriminating metric
vector loadings.

| Control | Options | Default |
|----|----|----|
| **Category Filter** | All Combined · (I) Dist. Props · (II) Corr. & Zero-Inflation · (III) Cellular Structure · (IV) Batch Effects · (V) Bio. Signal · (VII) Cross-Modal | All |
| **Panel View** | Both (Biplot + Loadings) · Simulators Only · Loadings Only | Both |
| **Top Metrics (Loadings)** | 4 – 50 (step 1) | 14 |
| **PCA Height (px)** | 500 – 1 600 (step 50) | 1 100 |
| **PCA Width (px)** | 650 – 1 400 (step 25) | 950 |

> When **All Categories Combined** is selected, only top-ranked metrics
> by loading magnitude \[√(PC1² + PC2²)\] are shown to prevent visual
> clutter.

### Sub-panel 7 — MDS Metric Space

[`plot_metric_mds()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_mds.md)
— Multi-Dimensional Scaling ordination.

| Control             | Options                             | Default       |
|---------------------|-------------------------------------|---------------|
| **Category Filter** | Same six options as PCA             | All           |
| **MDS Target**      | By Simulators · By Metric Summaries | By Simulators |
| **MDS Height (px)** | 400 – 1 300 (step 20)               | 620           |
| **MDS Width (px)**  | 600 – 1 400 (step 20)               | 900           |

### Sub-panel 8 — Cell Embeddings (t-SNE & UMAP)

[`compute_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/compute_dataset_embeddings.md) +
[`plot_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/plot_dataset_embeddings.md).

**Reduction & Layout:**

| Control | Options | Default |
|----|----|----|
| **Reduction Technique** | UMAP · t-SNE · PCA | UMAP |
| **Comparison Layout** | Faceted Grid (All Simulators) · Side-by-Side (vs Single Simulator) | Faceted Grid |
| **Select Simulator** (side-by-side only) | Dynamic dropdown | — |

**Colour & Simulator:**

| Control | Options | Default |
|----|----|----|
| **Simulator Selector** | Dynamic checkboxes | All |
| **Color Cells By** | Biological Cell Type / Group · Unsupervised Cluster · Sequencing Depth (Library Size) · Technical Batch · Dataset Source | Cell Type |

**Fine-Tuning:**

| Control                           | Range                 | Default |
|-----------------------------------|-----------------------|---------|
| **Number of PCs**                 | 5 – 50 (step 5)       | 20      |
| **t-SNE Perplexity** (t-SNE only) | 5 – 50 (step 5)       | 15      |
| **UMAP Neighbors** (UMAP only)    | 5 – 50 (step 5)       | 15      |
| **Point Size**                    | 0.2 – 3.0 (step 0.1)  | 1.0     |
| **Alpha (transparency)**          | 0.2 – 1.0 (step 0.05) | 0.8     |

**Quantitative Embedding Quality Table** (below the plot):

| Column                 | Description                                    |
|------------------------|------------------------------------------------|
| Dataset                | Reference or simulator name                    |
| Cells                  | Number of cells                                |
| Mean Silhouette        | Average silhouette width of cell-type clusters |
| ARI (Cluster Fidelity) | ARI between k-means clusters and known labels  |
| Mean Library Size      | Mean total UMI / peak count per cell           |
| Mean Detected Features | Mean genes / peaks detected per cell           |
| Library Size Diff (%)  | Absolute % difference vs. reference mean       |

------------------------------------------------------------------------

## Tab 5 — Download Results

### 1. All-in-One Benchmark Archive (.zip)

Contents of the single zip:

- **Excel Workbook (.xlsx)** — 4 sheets: All Benchmark Metrics, Dataset
  Properties, Method Rankings, Scalability Metrics
- **Master Table (.csv)** — all 62 metrics, tidy long format
- **R Object (.rds)** — for downstream R analysis
- **Multi-Page PDF Report** — all figures in one file
- **Individual JPEG Figures** — 600 DPI per visualization

Click **Download Complete Results (.zip)**.

### 2. Spreadsheets & Data Files

| Button                      | Format                     |
|-----------------------------|----------------------------|
| Download Excel File (.xlsx) | Multi-sheet workbook       |
| Download CSV Table (.csv)   | Flat tidy table            |
| Download R Data File (.rds) | Native R serialized object |

### 3. Complete Multi-Page PDF Report

Pages compiled into one PDF:

1.  Comparative Bubble Matrix
2.  Evaluation Summary
3.  Scalability Benchmark
4.  Metric Boxplots
5.  Metric Heatmap
6.  PCA Ordination
7.  MDS Metric Space
8.  Distribution QC Curves

Click **Download All Figures (.pdf)**.

### Interactive Benchmark Data Table

A full `DT` table with **per-column search boxes** lets you filter by
Method Name, Category, or Metric before exporting.

------------------------------------------------------------------------

## Tab 6 — Help & Getting Started

An in-app reference guide structured into seven sections:

1.  **Workflow Architecture** — 4-step pipeline diagram
2.  **Eight Evaluation Categories** — tabular reference for all 62
    metrics
3.  **Score Normalization Pipeline** — direction inversion → min-max
    scaling diagrams
4.  **Studio Navigation Guide** — accordion-style guide to each tab
5.  **Accepted Input Formats** — matrix types, annotation vectors, saved
    files
6.  **Frequently Asked Questions**
7.  **Team, Contact & Citation**

The tab also prominently links to the full online documentation at
<https://kabilanbio.github.io/scSimEval>.

------------------------------------------------------------------------

## R API Cross-Reference

Every Shiny control maps directly to a documented R function:

| Shiny Panel | R Function |
|----|----|
| Data Hub → Single-Cell | [`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md) |
| Data Hub → Multiomics | [`evaluate_multiomics_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_accuracy.md) |
| Data Hub → Multiple Datasets | [`evaluate_multiple_datasets()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiple_datasets.md) |
| Bubble Matrix | [`plot_benchmark_bubble_matrix()`](https://kabilanbio.github.io/scSimEval/reference/plot_benchmark_bubble_matrix.md) |
| Evaluation Summary | [`plot_evaluation_summary()`](https://kabilanbio.github.io/scSimEval/reference/plot_evaluation_summary.md) |
| Distribution QC | [`plot_distribution_qc()`](https://kabilanbio.github.io/scSimEval/reference/plot_distribution_qc.md) |
| Scalability Benchmark | [`plot_scalability_benchmark()`](https://kabilanbio.github.io/scSimEval/reference/plot_scalability_benchmark.md) |
| Metric Boxplots (category) | [`plot_metric_boxplots()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_boxplots.md) |
| Metric Boxplots (individual) | [`plot_individual_metric_bar()`](https://kabilanbio.github.io/scSimEval/reference/plot_individual_metric_bar.md) |
| Metric Heatmap | [`plot_metric_heatmap()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_heatmap.md) |
| PCA Ordination | [`plot_metric_pca()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_pca.md) |
| MDS Metric Space | [`plot_metric_mds()`](https://kabilanbio.github.io/scSimEval/reference/plot_metric_mds.md) |
| Cell Embeddings | [`compute_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/compute_dataset_embeddings.md) + [`plot_dataset_embeddings()`](https://kabilanbio.github.io/scSimEval/reference/plot_dataset_embeddings.md) |
| Embedding Quality Metrics | [`compute_embedding_quality_metrics()`](https://kabilanbio.github.io/scSimEval/reference/compute_embedding_quality_metrics.md) |

For complete parameter descriptions see the
[Reference](https://kabilanbio.github.io/scSimEval/reference/index.md)
index.

------------------------------------------------------------------------

    #> R version 4.6.1 (2026-06-24 ucrt)
    #> Platform: x86_64-w64-mingw32/x64
    #> Running under: Windows 11 x64 (build 26200)
    #> 
    #> Matrix products: default
    #>   LAPACK version 3.12.1
    #> 
    #> locale:
    #> [1] LC_COLLATE=English_India.utf8  LC_CTYPE=English_India.utf8   
    #> [3] LC_MONETARY=English_India.utf8 LC_NUMERIC=C                  
    #> [5] LC_TIME=English_India.utf8    
    #> 
    #> time zone: Asia/Calcutta
    #> tzcode source: internal
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
    #>  [5] xfun_0.60         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
    #>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
    #> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
    #> [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.12.0     
    #> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
    #> [25] rlang_1.3.0       fs_2.1.0          htmlwidgets_1.6.4
