library(testthat)
library(Matrix)

test_that("example_scrna dataset has valid structure and dimensions", {
  data("example_scrna", package = "scSimEval", envir = environment())
  expect_true(exists("example_scrna"))
  expect_type(example_scrna, "list")
  expect_named(example_scrna, c("ref", "sim", "cell_types", "batch_info"), ignore.order = TRUE)
  
  # Matrix dimensions (100 HVGs x 100 balanced Rice Bud cells)
  expect_equal(dim(example_scrna$ref), c(100, 100))
  expect_equal(dim(example_scrna$sim), c(100, 100))
  
  # Cell annotations
  expect_length(example_scrna$cell_types, 100)
  expect_length(example_scrna$batch_info, 100)
  expect_true(is.factor(example_scrna$cell_types))
  expect_true(is.factor(example_scrna$batch_info))
  expect_equal(length(levels(example_scrna$cell_types)), 8)
  expect_equal(length(levels(example_scrna$batch_info)), 2)
})

test_that("example_scatac dataset has valid structure and dimensions", {
  data("example_scatac", package = "scSimEval", envir = environment())
  expect_true(exists("example_scatac"))
  expect_type(example_scatac, "list")
  expect_named(example_scatac, c("ref", "sim", "cell_types", "batch_info"), ignore.order = TRUE)
  
  # Matrix dimensions (100 HVPs x 100 balanced Rice Bud cells)
  expect_equal(dim(example_scatac$ref), c(100, 100))
  expect_equal(dim(example_scatac$sim), c(100, 100))
  
  # Cell annotations
  expect_length(example_scatac$cell_types, 100)
  expect_length(example_scatac$batch_info, 100)
  expect_true(is.factor(example_scatac$cell_types))
  expect_true(is.factor(example_scatac$batch_info))
})

test_that("example_multiomics dataset has paired modalities and resource stats", {
  data("example_multiomics", package = "scSimEval", envir = environment())
  expect_true(exists("example_multiomics"))
  expect_type(example_multiomics, "list")
  expect_true(all(c("ref_multi", "sim_multi", "cell_types", "batch_info", "resource_stats") %in% names(example_multiomics)))
  
  expect_named(example_multiomics$ref_multi, c("rna", "atac"), ignore.order = TRUE)
  expect_named(example_multiomics$sim_multi, c("rna", "atac"), ignore.order = TRUE)
  
  # Cross-modality cell alignment
  expect_equal(ncol(example_multiomics$ref_multi$rna), ncol(example_multiomics$ref_multi$atac))
  expect_equal(colnames(example_multiomics$ref_multi$rna), colnames(example_multiomics$ref_multi$atac))
  
  # Resource stats
  expect_true(is.numeric(example_multiomics$resource_stats$elapsed_time))
  expect_true(is.numeric(example_multiomics$resource_stats$memory_mb))
})

test_that("example_rice_bud dataset has standardized SingleCellExperiment paired multiomics structure and 6 simulators", {
  data("example_rice_bud", package = "scSimEval", envir = environment())
  expect_true(exists("example_rice_bud"))
  expect_type(example_rice_bud, "list")
  expect_named(example_rice_bud, c("metadata", "reference", "simulated", "scalability_benchmarks"), ignore.order = TRUE)
  
  # Reference dimensions: standardized SingleCellExperiment objects (100 HVGs x 100 cells, 100 HVPs x 100 cells)
  expect_s4_class(example_rice_bud$reference$combined, "SingleCellExperiment")
  expect_s4_class(example_rice_bud$reference$rna, "SingleCellExperiment")
  expect_s4_class(example_rice_bud$reference$atac, "SingleCellExperiment")
  expect_true("ATAC" %in% SingleCellExperiment::altExpNames(example_rice_bud$reference$combined))
  expect_equal(dim(example_rice_bud$reference$rna), c(100, 100))
  expect_equal(dim(example_rice_bud$reference$atac), c(100, 100))
  expect_equal(ncol(example_rice_bud$reference$rna), ncol(example_rice_bud$reference$atac))
  
  # Embedded cell_type and batch (all 8 Rice Bud cell types covered equally)
  cd_rb <- SummarizedExperiment::colData(example_rice_bud$reference$combined)
  expect_true("cell_type" %in% colnames(cd_rb))
  expect_true("batch" %in% colnames(cd_rb))
  expect_equal(length(unique(cd_rb$cell_type)), 8)
  expect_equal(length(unique(cd_rb$batch)), 2)
  
  # 6 Simulators
  sim_methods <- c("scMultiSim", "scDesign3", "MOSim", "scCross", "scMoMtF", "Matilda")
  expect_true(all(sim_methods %in% names(example_rice_bud$simulated)))
  for (m in sim_methods) {
    expect_s4_class(example_rice_bud$simulated[[m]]$combined, "SingleCellExperiment")
    expect_equal(dim(example_rice_bud$simulated[[m]]$rna), c(100, 100))
    expect_equal(dim(example_rice_bud$simulated[[m]]$atac), c(100, 100))
  }
})

test_that("example_soybean_root dataset has standardized SingleCellExperiment unpaired multiomics structure and 3 simulators", {
  data("example_soybean_root", package = "scSimEval", envir = environment())
  expect_true(exists("example_soybean_root"))
  expect_type(example_soybean_root, "list")
  expect_named(example_soybean_root, c("metadata", "reference", "simulated", "scalability_benchmarks"), ignore.order = TRUE)
  
  # Reference dimensions: standardized SingleCellExperiment objects (100 HVGs x 100 RNA cells, 100 HVPs x 100 ATAC cells)
  expect_s4_class(example_soybean_root$reference$combined, "SingleCellExperiment")
  expect_s4_class(example_soybean_root$reference$rna, "SingleCellExperiment")
  expect_s4_class(example_soybean_root$reference$atac, "SingleCellExperiment")
  expect_true("ATAC" %in% SingleCellExperiment::altExpNames(example_soybean_root$reference$combined))
  expect_equal(dim(example_soybean_root$reference$rna), c(100, 100))
  expect_equal(dim(example_soybean_root$reference$atac), c(100, 100))
  
  # 3 Simulators (MOSim excluded because it only supports paired co-assay simulation)
  sim_methods <- c("scMultiSim", "scDesign3", "scCross")
  expect_equal(length(example_soybean_root$simulated), 3)
  expect_true(all(sim_methods %in% names(example_soybean_root$simulated)))
  expect_false("MOSim" %in% names(example_soybean_root$simulated))
  for (m in sim_methods) {
    expect_s4_class(example_soybean_root$simulated[[m]]$combined, "SingleCellExperiment")
    expect_equal(dim(example_soybean_root$simulated[[m]]$rna), c(100, 100))
    expect_equal(dim(example_soybean_root$simulated[[m]]$atac), c(100, 100))
  }
})

test_that("tfidf_transform and scATAC-seq peak-peak correlations use TF-IDF (simPIC standard)", {
  set.seed(123)
  mat <- matrix(rpois(30 * 40, lambda = 0.6), nrow = 30, ncol = 40)
  mat[1, ] <- 0 # zero-sum peak row guard
  rownames(mat) <- paste0("chr1:", seq(100, by = 500, length.out = 30), "-", seq(400, by = 500, length.out = 30))
  sp_mat <- methods::as(mat, "dgCMatrix")

  tf_dense <- tfidf_transform(mat)
  tf_sparse <- tfidf_transform(sp_mat)
  expect_equal(tf_dense, as.matrix(tf_sparse), tolerance = 1e-12)
  expect_false(any(is.na(tf_dense) | is.infinite(tf_dense)))
  expect_equal(tf_dense[1, ], rep(0, 40))

  # Verify exact Signac method=1 formula: log1p((x / colSum) * (ncol / rowSum) * 1e4)
  cs <- colSums(mat)
  rs <- rowSums(mat)
  expected_2_3 <- log1p((mat[2, 3] / cs[3]) * (ncol(mat) / rs[2]) * 1e4)
  expect_equal(unname(tf_dense[2, 3]), unname(expected_2_3), tolerance = 1e-12)

  # Verify extract_feature_properties automatically uses TF-IDF on peaks vs log2-CPM on RNA
  fp_atac_auto <- extract_feature_properties(sp_mat, n_top_cor = 20)
  fp_atac_explicit <- extract_feature_properties(sp_mat, n_top_cor = 20, modality = "scATAC-seq")
  fp_rna_forced <- extract_feature_properties(sp_mat, n_top_cor = 20, modality = "scRNA-seq")

  expect_equal(fp_atac_auto$gene_cor, fp_atac_explicit$gene_cor, tolerance = 1e-12)
  expect_equal(unname(fp_atac_auto$mean_expression), unname(rowMeans(tf_dense)), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(fp_atac_auto$mean_expression, fp_rna_forced$mean_expression)))
  expect_false(isTRUE(all.equal(fp_atac_auto$gene_cor, fp_rna_forced$gene_cor)))
})
