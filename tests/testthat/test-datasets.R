library(testthat)
library(Matrix)

test_that("example_scrna dataset has valid structure and dimensions", {
  data("example_scrna", package = "scSimEval", envir = environment())
  expect_true(exists("example_scrna"))
  expect_type(example_scrna, "list")
  expect_named(example_scrna, c("ref", "sim", "cell_types", "batch_info"), ignore.order = TRUE)
  
  # Matrix dimensions
  expect_equal(dim(example_scrna$ref), c(60, 80))
  expect_equal(dim(example_scrna$sim), c(60, 80))
  
  # Cell annotations
  expect_length(example_scrna$cell_types, 80)
  expect_length(example_scrna$batch_info, 80)
  expect_true(is.factor(example_scrna$cell_types))
  expect_true(is.factor(example_scrna$batch_info))
  expect_equal(levels(example_scrna$cell_types), c("TypeA", "TypeB"))
  expect_equal(levels(example_scrna$batch_info), c("Batch1", "Batch2"))
})

test_that("example_scatac dataset has valid structure and dimensions", {
  data("example_scatac", package = "scSimEval", envir = environment())
  expect_true(exists("example_scatac"))
  expect_type(example_scatac, "list")
  expect_named(example_scatac, c("ref", "sim", "cell_types", "batch_info"), ignore.order = TRUE)
  
  # Matrix dimensions
  expect_equal(dim(example_scatac$ref), c(60, 80))
  expect_equal(dim(example_scatac$sim), c(60, 80))
  
  # Cell annotations
  expect_length(example_scatac$cell_types, 80)
  expect_length(example_scatac$batch_info, 80)
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
  expect_true(is.numeric(example_multiomics$resource_stats$cpu_time))
  expect_true(is.numeric(example_multiomics$resource_stats$memory_mb))
})

test_that("example_rice_bud dataset has real paired multiomics structure and 6 simulators", {
  data("example_rice_bud", package = "scSimEval", envir = environment())
  expect_true(exists("example_rice_bud"))
  expect_type(example_rice_bud, "list")
  expect_named(example_rice_bud, c("metadata", "reference", "simulated", "scalability_benchmarks"), ignore.order = TRUE)
  
  # Reference dimensions: 200 genes x 100 cells, 200 peaks x 100 cells
  expect_equal(dim(example_rice_bud$reference$rna), c(200, 100))
  expect_equal(dim(example_rice_bud$reference$atac), c(200, 100))
  expect_equal(ncol(example_rice_bud$reference$rna), ncol(example_rice_bud$reference$atac))
  
  # Metadata
  expect_equal(nrow(example_rice_bud$reference$cell_types), 100)
  expect_equal(nrow(example_rice_bud$reference$batches), 100)
  
  # 6 Simulators
  sim_methods <- c("scMultiSim", "scDesign3", "MOSim", "scCross", "scMoMtF", "Matilda")
  expect_true(all(sim_methods %in% names(example_rice_bud$simulated)))
  for (m in sim_methods) {
    expect_equal(dim(example_rice_bud$simulated[[m]]$rna), c(200, 100))
    expect_equal(dim(example_rice_bud$simulated[[m]]$atac), c(200, 100))
  }
})

test_that("example_soybean_root dataset has real unpaired multiomics structure and 3 simulators", {
  data("example_soybean_root", package = "scSimEval", envir = environment())
  expect_true(exists("example_soybean_root"))
  expect_type(example_soybean_root, "list")
  expect_named(example_soybean_root, c("metadata", "reference", "simulated", "scalability_benchmarks"), ignore.order = TRUE)
  
  # Reference dimensions: 200 genes x 100 cells, 200 peaks x 100 cells
  expect_equal(dim(example_soybean_root$reference$rna), c(200, 100))
  expect_equal(dim(example_soybean_root$reference$atac), c(200, 100))
  
  # 3 Simulators (MOSim excluded because it only supports paired co-assay simulation)
  sim_methods <- c("scMultiSim", "scDesign3", "scCross")
  expect_equal(length(example_soybean_root$simulated), 3)
  expect_true(all(sim_methods %in% names(example_soybean_root$simulated)))
  expect_false("MOSim" %in% names(example_soybean_root$simulated))
  for (m in sim_methods) {
    expect_equal(dim(example_soybean_root$simulated[[m]]$rna), c(200, 100))
    expect_equal(dim(example_soybean_root$simulated[[m]]$atac), c(200, 100))
  }
})
