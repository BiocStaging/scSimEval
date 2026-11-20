#' Append Rows from a Summary Data Frame with Category/Property Tagging
#' @keywords internal
#' @noRd
.append_rows_from_df <- function(df, category, property, modality = NULL) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(list())
  out <- data.frame(
    Category = category,
    Property = property,
    Metric = as.character(df$Metric),
    Value = as.numeric(df$Value),
    stringsAsFactors = FALSE
  )
  if (!is.null(modality)) out$Modality <- modality
  list(out)
}

#' Append Scalar Metrics from a Named List or Vector
#' @keywords internal
#' @noRd
.append_rows_from_fields <- function(obj, fields, category, property, modality = NULL) {
  if (is.null(obj)) return(list())
  rows <- list()
  for (m in fields) {
    val <- obj[[m]]
    if (!is.null(val) && length(val) == 1L && is.finite(val)) {
      df <- data.frame(
        Category = category,
        Property = property,
        Metric = m,
        Value = as.numeric(val),
        stringsAsFactors = FALSE
      )
      if (!is.null(modality)) df$Modality <- modality
      rows[[length(rows) + 1L]] <- df
    }
  }
  rows
}

#' Append Category VII dynverse Trajectory Metrics to Summary Rows
#' @keywords internal
#' @noRd
.append_trajectory_rows <- function(traj_res, modality = NULL) {
  if (is.null(traj_res)) return(list())
  spec <- list(
    list(field = "him", prop = "trajectory_topology"),
    list(field = "edge_flip", prop = "trajectory_topology"),
    list(field = "geodesic_correlation", prop = "geodesic_ordering"),
    list(field = "F1_branches", prop = "branch_concordance"),
    list(field = "wcor_features", prop = "dynamic_features")
  )
  rows <- list()
  for (s in spec) {
    val <- traj_res[[s$field]]
    if (!is.null(val) && length(val) == 1L && is.finite(val)) {
      df <- data.frame(
        Category = "(VII) Trajectory & Lineage Topology",
        Property = s$prop,
        Metric = s$field,
        Value = as.numeric(val),
        stringsAsFactors = FALSE
      )
      if (!is.null(modality)) df$Modality <- modality
      rows[[length(rows) + 1L]] <- df
    }
  }
  rows
}

#' Filter Rare Cell Types (< 10 Cells) for SimBench Compliance
#' @keywords internal
#' @noRd
.filter_rare_cell_types <- function(counts, cell_types, batch_info, label = "reference", verbose = TRUE) {
  if (is.null(cell_types) || length(cell_types) != ncol(counts)) {
    return(list(counts = counts, cell_types = cell_types, batch_info = batch_info))
  }
  tab <- table(cell_types)
  rare <- names(tab)[tab < 10L]
  if (length(rare) > 0L && length(tab) > length(rare)) {
    keep <- which(!cell_types %in% rare)
    if (length(keep) >= 10L) {
      if (verbose) {
        message(sprintf(
          "SimBench standard: skipping cell type(s) with < 10 cells in %s: %s",
          label, paste(rare, collapse = ", ")
        ))
      }
      counts <- counts[, keep, drop = FALSE]
      cell_types <- cell_types[keep]
      if (!is.null(batch_info) && length(batch_info) >= max(keep)) {
        batch_info <- batch_info[keep]
      }
    }
  }
  list(counts = counts, cell_types = cell_types, batch_info = batch_info)
}

#' Extract Cell- and Feature-Level Properties for Reference and Simulation
#' @keywords internal
#' @noRd
.compute_unimodal_properties <- function(ref_counts, sim_counts, modality_resolved, threads) {
  if (threads >= 2L && .Platform$OS.type != "windows") {
    inner_th <- max(1L, as.integer(floor(threads / 4L)))
    prop_ext <- tryCatch(
      .sc_parallel_lapply(seq_len(4L), function(task_i) {
        if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
          try(RhpcBLASctl::blas_set_num_threads(inner_th), silent = TRUE)
          try(RhpcBLASctl::omp_set_num_threads(inner_th), silent = TRUE)
        }
        switch(
          task_i,
          extract_cell_properties(ref_counts, modality = modality_resolved, verbose = FALSE),
          extract_cell_properties(sim_counts, modality = modality_resolved, verbose = FALSE),
          extract_feature_properties(ref_counts, modality = modality_resolved, verbose = FALSE),
          extract_feature_properties(sim_counts, modality = modality_resolved, verbose = FALSE)
        )
      }, threads = min(threads, 4L)),
      error = function(e) NULL
    )
    .set_parallel_threads(threads)
    if (!is.null(prop_ext) && length(prop_ext) == 4L &&
        !any(vapply(prop_ext, is.null, logical(1)))) {
      return(list(
        ref_cell = prop_ext[[1L]],
        sim_cell = prop_ext[[2L]],
        ref_feat = prop_ext[[3L]],
        sim_feat = prop_ext[[4L]]
      ))
    }
  }
  list(
    ref_cell = extract_cell_properties(ref_counts, modality = modality_resolved, verbose = FALSE),
    sim_cell = extract_cell_properties(sim_counts, modality = modality_resolved, verbose = FALSE),
    ref_feat = extract_feature_properties(ref_counts, modality = modality_resolved, verbose = FALSE),
    sim_feat = extract_feature_properties(sim_counts, modality = modality_resolved, verbose = FALSE)
  )
}

#' Run Univariate Cell and Feature Metric Evaluations
#' @keywords internal
#' @noRd
.run_univariate_eval <- function(ref_cell, sim_cell, ref_feat, sim_feat, is_atac, threads) {
  has_tmm <- !is_atac && !is.null(ref_cell$tmm_factor) &&
             length(ref_cell$tmm_factor) > 0L && !all(is.na(ref_cell$tmm_factor))
  cell_props <- if (has_tmm) {
    c("library_size", "log_library_size", "tmm_factor",
      "zero_fraction_cell", "detection_freq_cell", "cell_cor", "pca_dist")
  } else {
    c("library_size", "log_library_size",
      "zero_fraction_cell", "detection_freq_cell", "cell_cor", "pca_dist")
  }
  feat_props <- c("mean_expression", "variance", "scaled_variance", "cv",
                  "zero_fraction_feature", "dispersion", "gene_cor")

  all_props <- c(
    lapply(cell_props, function(p) list(type = "cell", name = p, r = ref_cell[[p]], s = sim_cell[[p]])),
    lapply(feat_props, function(p) list(type = "feat", name = p, r = ref_feat[[p]], s = sim_feat[[p]]))
  )
  all_props <- all_props[vapply(all_props, function(x) !is.null(x$r) && !is.null(x$s) && length(x$r) > 0L, logical(1))]

  prop_res <- tryCatch(
    .sc_parallel_lapply(all_props, function(item) calc_all_univariate_metrics(item$r, item$s, item$name), threads = threads),
    error = function(e) NULL
  )
  if (is.null(prop_res) || length(prop_res) != length(all_props)) {
    prop_res <- lapply(all_props, function(item) calc_all_univariate_metrics(item$r, item$s, item$name))
  }

  cell_results <- list()
  feat_results <- list()
  for (idx in seq_along(all_props)) {
    item <- all_props[[idx]]
    if (item$type == "cell") cell_results[[item$name]] <- prop_res[[idx]]
    else feat_results[[item$name]] <- prop_res[[idx]]
  }
  list(cell_results = cell_results, feat_results = feat_results)
}

#' Run Bivariate Metric Evaluations
#' @keywords internal
#' @noRd
.run_bivariate_eval <- function(ref_cell, sim_cell, ref_feat, sim_feat, threads) {
  biv_configs <- list(
    lib_vs_cellzero = list(
      ref = cbind(log10(ref_cell$library_size + 1), ref_cell$zero_fraction_cell),
      sim = cbind(log10(sim_cell$library_size + 1), sim_cell$zero_fraction_cell),
      prefix = "lib_vs_cellzero"
    ),
    mean_vs_sd = list(
      ref = cbind(ref_feat$mean_expression, ref_feat$sd),
      sim = cbind(sim_feat$mean_expression, sim_feat$sd),
      prefix = "mean_vs_sd"
    ),
    mean_vs_featzero = list(
      ref = cbind(ref_feat$mean_expression, ref_feat$zero_fraction_feature),
      sim = cbind(sim_feat$mean_expression, sim_feat$zero_fraction_feature),
      prefix = "mean_vs_featzero"
    ),
    mean_vs_dispersion = list(
      ref = cbind(ref_feat$mean_expression, ref_feat$dispersion),
      sim = cbind(sim_feat$mean_expression, sim_feat$dispersion),
      prefix = "mean_vs_dispersion"
    )
  )

  if (threads > 1L) {
    n_workers <- min(threads, length(biv_configs))
    inner_th <- max(1L, as.integer(floor(threads / n_workers)))
    biv_list <- tryCatch(
      .sc_parallel_lapply(biv_configs, function(cfg) {
        .set_parallel_threads(inner_th)
        calc_all_bivariate_metrics(cfg$ref, cfg$sim, cfg$prefix, threads = inner_th)
      }, threads = n_workers),
      error = function(e) NULL
    )
    .set_parallel_threads(threads)
    if (!is.null(biv_list) && length(biv_list) == length(biv_configs)) {
      return(stats::setNames(biv_list, names(biv_configs)))
    }
  }

  stats::setNames(
    lapply(biv_configs, function(cfg) calc_all_bivariate_metrics(cfg$ref, cfg$sim, cfg$prefix, threads = threads)),
    names(biv_configs)
  )
}

#' Evaluate Category II Sparsity Curves and Category III Manifold Distances
#' @keywords internal
#' @noRd
.run_sparsity_manifold_eval <- function(ref_counts, sim_counts, is_atac, threads) {
  run_curve <- function() {
    if (!is_atac) evaluate_zero_probability_curve(ref_counts, sim_counts)
    else evaluate_accessibility_sparsity_curve(ref_counts, sim_counts)
  }

  sp_res <- NULL
  if (threads >= 2L && .Platform$OS.type != "windows") {
    sp_res <- tryCatch(
      .sc_parallel_lapply(seq_len(5L), function(task_i) {
        switch(
          task_i,
          run_curve(),
          calc_excess_zero_weights(ref_counts),
          calc_excess_zero_weights(sim_counts),
          calc_generative_precision_recall(ref_counts, sim_counts),
          calc_mmd(ref_counts, sim_counts)
        )
      }, threads = min(threads, 5L)),
      error = function(e) NULL
    )
  }

  if (!is.null(sp_res) && length(sp_res) == 5L) {
    curve_obj <- sp_res[[1L]]
    ref_ex <- sp_res[[2L]]
    sim_ex <- sp_res[[3L]]
    gen_pr <- sp_res[[4L]]
    mmd_val <- sp_res[[5L]]
  } else {
    curve_obj <- run_curve()
    ref_ex <- calc_excess_zero_weights(ref_counts)
    sim_ex <- calc_excess_zero_weights(sim_counts)
    gen_pr <- calc_generative_precision_recall(ref_counts, sim_counts)
    mmd_val <- calc_mmd(ref_counts, sim_counts)
  }

  list(
    zero_prob_eval = if (!is_atac) curve_obj else NULL,
    acc_sparsity_eval = if (is_atac) curve_obj else NULL,
    ref_excess_w = ref_ex,
    sim_excess_w = sim_ex,
    gen_pr = gen_pr,
    mmd_eval = mmd_val
  )
}

#' Convert Core Distributional, Sparsity, and Manifold Results into Tidy Summary Rows
#' @keywords internal
#' @noRd
.build_unimodal_summary_rows <- function(cell_results, feat_results, bivariate_results, sp_m) {
  univ_metrics <- c("Wasserstein", "KS", "MAD", "Bhattacharyya", "Runs_Statistic",
                    "NN_Mismatch", "Between_Dataset_Silh_Global", "MAE", "RMSE", "KDE_ISE")

  append_prop_block <- function(res_list, metric_names) {
    out <- list()
    for (prop in names(res_list)) {
      for (m in metric_names) {
        val <- res_list[[prop]][[paste0(prop, "_", m)]]
        if (!is.null(val)) {
          out[[length(out) + 1L]] <- data.frame(
            Category = "Distributional Properties",
            Property = prop,
            Metric = m,
            Value = as.numeric(val),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    out
  }

  rows <- c(
    append_prop_block(cell_results, univ_metrics),
    append_prop_block(feat_results, univ_metrics)
  )
  if (length(bivariate_results) > 0L) {
    rows <- c(rows, append_prop_block(bivariate_results, c(
      "Fasano_Franceschini_2D_KS", "EMD_2D", "NN_Mismatch_2D", "Between_Dataset_Silh_2D_Global"
    )))
  }

  if (!is.null(sp_m$zero_prob_eval)) {
    rows <- c(rows, .append_rows_from_fields(
      sp_m$zero_prob_eval, names(sp_m$zero_prob_eval),
      "Sparsity & Zero-Inflation Dynamics", "dropout_curve"
    ))
  }
  if (!is.null(sp_m$acc_sparsity_eval)) {
    rows <- c(rows, .append_rows_from_fields(
      sp_m$acc_sparsity_eval, c("curve_rmse", "r2_discrepancy"),
      "Sparsity & Zero-Inflation Dynamics", "accessibility_sparsity_curve"
    ))
  }
  rows[[length(rows) + 1L]] <- data.frame(
    Category = "Sparsity & Zero-Inflation Dynamics",
    Property = "excess_zero_weights",
    Metric = "excess_zero_weight_discrepancy",
    Value = abs(sp_m$ref_excess_w$mean_excess_zero_weight - sp_m$sim_excess_w$mean_excess_zero_weight),
    stringsAsFactors = FALSE
  )
  rows <- c(rows, .append_rows_from_fields(
    list(MMD = sp_m$mmd_eval,
         Generative_Precision = sp_m$gen_pr$generative_precision,
         Generative_Recall = sp_m$gen_pr$generative_recall),
    c("MMD", "Generative_Precision", "Generative_Recall"),
    "Cellular Manifold & Clustering Concordance", "manifold"
  ))
  rows
}

#' Evaluate Category III Clustering, Category IV Batch Mixing, and Category V Biological Signal
#' @keywords internal
#' @noRd
.evaluate_unimodal_downstream <- function(ref_counts, sim_counts,
                                          cell_types_ref, cell_types_sim,
                                          batch_ref, batch_sim,
                                          is_atac, threads, verbose) {
  rows <- list()
  clu_eval <- tryCatch(
    evaluate_clustering_metrics(
      data = sim_counts, cell_types = cell_types_sim,
      ref_data = ref_counts, ref_cell_types = cell_types_ref,
      batch_info = batch_sim, threads = threads
    ),
    error = function(e) NULL
  )
  rows <- c(rows, .append_rows_from_fields(
    clu_eval, c("silhouette_sim", "cdi", "neighborhood_purity", "ari", "nmi"),
    "Cellular Manifold & Clustering Concordance", "clustering"
  ))

  b_eval <- NULL
  if (!is.null(batch_sim)) {
    if (verbose) message("Evaluating batch integration & confounder mixing...")
    b_eval <- tryCatch(
      evaluate_batch_metrics(data = sim_counts, batch_info = batch_sim,
                             cell_types = cell_types_sim, threads = threads, verbose = FALSE),
      error = function(e) NULL
    )
    rows <- c(rows, .append_rows_from_fields(
      b_eval, c("batch_silhouette", "pcr_r2", "cms", "isi", "cross_batch_F1"),
      "Technical Variation & Batch Mixing Dynamics", "batch_mixing"
    ))
  }

  signal_eval <- NULL
  if (!is.null(cell_types_ref) && !is.null(cell_types_sim) &&
      !is.null(batch_ref) && !is.null(batch_sim)) {
    if (is_atac) {
      if (verbose) message("Evaluating scATAC-seq biological signal fidelity (Libra Pseudobulk edgeR-QLF)...")
      signal_eval <- tryCatch(
        evaluate_atac_da_fidelity(
          ref_data = ref_counts, sim_data = sim_counts,
          ref_celltypes = cell_types_ref, sim_celltypes = cell_types_sim,
          batch_ref = batch_ref, batch_sim = batch_sim, threads = threads
        ),
        error = function(e) {
          if (verbose) message("scATAC-seq DA evaluation note: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(signal_eval)) {
        rows <- c(rows, .append_rows_from_df(
          signal_eval$da_summary_table,
          "Biological Signal & Differential Fidelity", "scATAC_differential_accessibility"
        ))
      }
    } else {
      if (verbose) message("Evaluating scRNA-seq biological signal fidelity (Libra Pseudobulk edgeR-QLF)...")
      signal_eval <- tryCatch(
        evaluate_deg_fidelity(
          ref_data = ref_counts, sim_data = sim_counts,
          ref_celltypes = cell_types_ref, sim_celltypes = cell_types_sim,
          batch_ref = batch_ref, batch_sim = batch_sim, threads = threads
        ),
        error = function(e) {
          if (verbose) message("scRNA-seq pseudobulk DE evaluation note: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(signal_eval)) {
        rows <- c(rows, .append_rows_from_df(
          signal_eval$deg_summary_table,
          "Biological Signal & Differential Fidelity", "scRNA_differential_expression"
        ))
      }
    }
  } else if (!is.null(cell_types_ref) && !is.null(cell_types_sim) && verbose) {
    mod_str <- if (is_atac) "scATAC-seq" else "scRNA-seq"
    message(sprintf("Biological replicate annotations ('%s') not provided: skipping Category V pseudobulk analysis for %s.", "batch", mod_str))
  }

  list(rows = rows, clustering = clu_eval, batch_effects = b_eval, biological_signal = signal_eval)
}

#' Compute Unimodal Computational Scalability Metrics (Category VIII)
#' @keywords internal
#' @noRd
.compute_unimodal_scalability <- function(elapsed_time, memory_mb, time_unit, memory_unit, sim_cores, verbose) {
  time_unit_norm <- .normalize_time_unit(time_unit)
  memory_unit_norm <- .normalize_memory_unit(memory_unit)
  elapsed_seconds <- .convert_time_to_seconds(elapsed_time, time_unit_norm)
  memory_mib <- .convert_memory_to_mib(memory_mb, memory_unit_norm)
  sim_cores_val <- .normalize_sim_cores(sim_cores)
  cores_mult <- if (!is.null(sim_cores_val)) sim_cores_val else 1L
  core_sec_val <- if (!is.null(elapsed_seconds) && !is.na(elapsed_seconds)) as.numeric(elapsed_seconds * cores_mult) else NULL
  mem_core_val <- if (!is.null(memory_mib) && !is.na(memory_mib)) as.numeric(memory_mib / cores_mult) else NULL

  rows <- list()
  if (!is.null(core_sec_val) || !is.null(mem_core_val)) {
    if (verbose) message("[7/7] Logging computational scalability & efficiency metrics...")
    rows <- .append_rows_from_fields(
      list(core_seconds = core_sec_val, memory_per_core_mib = mem_core_val),
      c("core_seconds", "memory_per_core_mib"),
      "(VIII) Computational Scalability & Efficiency", "resource_usage"
    )
  }

  list(
    rows = rows,
    sim_cores_val = sim_cores_val,
    resource_usage = list(
      core_seconds = core_sec_val,
      memory_per_core_mib = mem_core_val,
      memory_per_core_gib = if (!is.null(mem_core_val)) mem_core_val / 1024 else NULL,
      memory_mb = memory_mib,
      memory_input = memory_mb,
      memory_unit = memory_unit_norm,
      elapsed_time = elapsed_seconds,
      elapsed_time_input = elapsed_time,
      time_unit = time_unit_norm,
      sim_cores = sim_cores_val
    )
  )
}

#' Evaluate Comprehensive Simulation Accuracy (Unimodal SingleCellExperiment)
#'
#' Computes distribution distances, cellular and feature-level properties,
#' zero-inflation curves, clustering concordance, batch mixing, differential
#' expression/accessibility fidelity, trajectory topology, and computational
#' scalability between reference and simulated \code{SingleCellExperiment} objects.
#' Users starting from \code{Seurat} objects should convert via
#' \code{Seurat::as.SingleCellExperiment()}.
#'
#' \if{html}{\figure{workflow_diagram.png}{options: width="800" alt="scSimEval Unified Benchmarking Workflow Diagram"}}
#'
#' @details
#' \strong{Metric Measurement Capacity across 4 Analysis Options (56 Evaluation Measures Suite):}
#' Cell-type and batch annotations are retrieved directly from \code{colData(ref_sce)} and
#' \code{colData(sim_sce)} using \code{cell_type_col} and \code{batch_col}.
#' \itemize{
#'   \item \strong{1. scRNA-seq (45 Metrics):} Categories I, II, III, IV, V, VII, and VIII.
#'   \item \strong{2. scATAC-seq (32 Metrics):} Categories I, II, III, IV, V, and VIII (trajectory skipped for ATAC).
#'   \item \strong{3. Paired Multiomics (56 Metrics):} All 8 categories (I--VIII) via \code{\link{evaluate_multiomics_accuracy}}.
#'   \item \strong{4. Unpaired Multiomics (51 Metrics):} All 8 categories omitting 5 cell-matched pairing metrics.
#' }
#'
#' @param ref_sce Reference \code{SingleCellExperiment} object.
#' @param sim_sce Simulated \code{SingleCellExperiment} object.
#' @param cell_type_col Column name in \code{colData} containing cell-type labels (default \code{"cell_type"}).
#' @param batch_col Column name in \code{colData} containing batch/replicate labels (default \code{"batch"}).
#' @param modality Optional modality hint (\code{"scRNA-seq"} or \code{"scATAC-seq"}); auto-detected if \code{NULL}.
#' @param compute_bivariate Logical, whether to compute 2D bivariate tests. Default is \code{TRUE}.
#' @param threads Integer number of parallel worker threads. Default is \code{1L}.
#' @param memory_mb Optional numeric peak memory allocation (in \code{memory_unit}).
#' @param elapsed_time Optional numeric wall-clock elapsed runtime (in \code{time_unit}).
#' @param time_unit Unit of \code{elapsed_time}: \code{"seconds"} (default) or \code{"minutes"}.
#' @param memory_unit Unit of \code{memory_mb}: \code{"MiB"} (default) or \code{"MB"}.
#' @param sim_cores Optional positive integer specifying CPU cores used by the simulator.
#' @param trajectory_num_trees Number of Random Forest trees for Category VII \code{wcor_features} (default \code{10000L}).
#' @param include_trajectory Logical; whether to evaluate Category VII trajectory metrics for scRNA-seq. Default \code{TRUE}.
#' @param include_downstream Logical; whether to evaluate Categories III, IV, and V downstream metrics. Default \code{TRUE}.
#' @param verbose Logical, whether to print progress messages. Default is \code{TRUE}.
#'
#' @return A named list containing tidy summaries and detailed category results:
#' \describe{
#'   \item{metrics_summary_table}{Data frame with columns \code{Category}, \code{Property}, \code{Metric}, \code{Value} (and \code{Sim_Cores} if provided).}
#'   \item{cell_metrics}{Named list of per-property univariate metric lists at cell resolution.}
#'   \item{feature_metrics}{Named list of per-property univariate metric lists at feature resolution.}
#'   \item{bivariate_metrics}{Named list of 2D bivariate metric comparisons.}
#'   \item{zero_probability_curve}{Category II scRNA-seq dropout curve results, or \code{NULL}.}
#'   \item{accessibility_sparsity_curve}{Category II scATAC-seq sparsity curve results, or \code{NULL}.}
#'   \item{excess_zero_weights}{List with \code{ref} and \code{sim} excess-zero probability weights.}
#'   \item{generative_manifold}{List with \code{mmd} and \code{precision_recall} generative manifold scores.}
#'   \item{ref_properties}{List containing \code{cell} and \code{feature} property vectors extracted from \code{ref_sce}.}
#'   \item{sim_properties}{List containing \code{cell} and \code{feature} property vectors extracted from \code{sim_sce}.}
#'   \item{batch_effects}{Category IV batch mixing metrics list, or \code{NULL}.}
#'   \item{biological_signal}{Category V pseudobulk DE/DA fidelity results, or \code{NULL}.}
#'   \item{trajectory}{Category VII \emph{dynverse} trajectory evaluation results, or \code{NULL}.}
#'   \item{resource_usage}{Normalized Category VIII scalability list (\code{core_seconds}, \code{memory_per_core_mib}, \code{memory_mb}, \code{elapsed_time}, \code{sim_cores}).}
#' }
#' @examples
#' data(example_scrna, package = "scSimEval")
#' res <- evaluate_simulation_accuracy(
#'   example_scrna$ref, example_scrna$sim,
#'   include_trajectory = FALSE, verbose = FALSE
#' )
#' head(res$metrics_summary_table, 10)
#' @export
evaluate_simulation_accuracy <- function(
  ref_sce,
  sim_sce,
  cell_type_col = "cell_type",
  batch_col = "batch",
  modality = NULL,
  compute_bivariate = TRUE,
  threads = 1L,
  memory_mb = NULL,
  elapsed_time = NULL,
  time_unit = c("seconds", "minutes"),
  memory_unit = c("MiB", "MB"),
  sim_cores = NULL,
  trajectory_num_trees = 10000L,
  include_trajectory = TRUE,
  include_downstream = TRUE,
  verbose = TRUE
) {
  val_ref <- validate_sc_object(ref_sce, "ref_sce", modality = modality,
                                cell_type_col = cell_type_col, batch_col = batch_col)
  val_sim <- validate_sc_object(sim_sce, "sim_sce", modality = modality,
                                cell_type_col = cell_type_col, batch_col = batch_col)

  f_ref <- .filter_rare_cell_types(val_ref$counts, val_ref$cell_types, val_ref$batch_info, "reference", verbose)
  f_sim <- .filter_rare_cell_types(val_sim$counts, val_sim$cell_types, val_sim$batch_info, "simulation", verbose)
  ref_counts <- f_ref$counts
  sim_counts <- f_sim$counts
  cell_types_ref <- f_ref$cell_types
  cell_types_sim <- f_sim$cell_types
  batch_ref <- f_ref$batch_info
  batch_sim <- f_sim$batch_info

  threads <- .set_parallel_threads(threads)
  is_atac <- .is_atac_modality(val_ref$object, modality) || .is_atac_modality(val_sim$object, modality)
  modality_resolved <- if (is_atac) "scATAC-seq" else "scRNA-seq"

  if (verbose) message("[1/5] Extracting cell-level and feature-level properties...")
  props <- .compute_unimodal_properties(ref_counts, sim_counts, modality_resolved, threads)

  if (verbose) message("[3/5] Computing univariate accuracy metrics...")
  univ <- .run_univariate_eval(props$ref_cell, props$sim_cell, props$ref_feat, props$sim_feat, is_atac, threads)

  bivariate_results <- list()
  if (compute_bivariate) {
    if (verbose) message("[4/5] Computing bivariate metrics (Fasano-Franceschini, KDE zstat, 2D EMD)...")
    bivariate_results <- .run_bivariate_eval(props$ref_cell, props$sim_cell, props$ref_feat, props$sim_feat, threads)
  }

  if (verbose) message("[5/5] Computing zero-probability and manifold distances...")
  sp_m <- .run_sparsity_manifold_eval(ref_counts, sim_counts, is_atac, threads)
  summary_rows <- .build_unimodal_summary_rows(univ$cell_results, univ$feat_results, bivariate_results, sp_m)

  downstream <- list(rows = list(), batch_effects = NULL, biological_signal = NULL)
  if (isTRUE(include_downstream)) {
    downstream <- .evaluate_unimodal_downstream(
      ref_counts, sim_counts, cell_types_ref, cell_types_sim, batch_ref, batch_sim,
      is_atac, threads, verbose
    )
    summary_rows <- c(summary_rows, downstream$rows)
  }

  traj_res <- NULL
  if (!is_atac && isTRUE(include_trajectory)) {
    if (verbose) message("[6/7] Evaluating trajectory and lineage topology (HVG + dynverse, scRNA-seq only)...")
    traj_res <- tryCatch(
      evaluate_trajectory_metrics(
        ref_data = ref_counts, sim_data = sim_counts,
        cell_types_ref = cell_types_ref, cell_types_sim = cell_types_sim,
        num_trees = trajectory_num_trees, threads = threads, verbose = FALSE
      ),
      error = function(e) NULL
    )
    summary_rows <- c(summary_rows, .append_trajectory_rows(traj_res, modality = NULL))
  } else if (is_atac && verbose) {
    message("[6/7] Category VII (Trajectory & Lineage Topology) requires scRNA-seq data: skipping for scATAC-seq.")
  }

  scal <- .compute_unimodal_scalability(elapsed_time, memory_mb, time_unit, memory_unit, sim_cores, verbose)
  summary_rows <- c(summary_rows, scal$rows)

  summary_table <- do.call(rbind, summary_rows)
  if (!is.null(scal$sim_cores_val) && !is.null(summary_table) && nrow(summary_table) > 0L) {
    summary_table$Sim_Cores <- scal$sim_cores_val
  }
  if (verbose) message("Unimodal accuracy evaluation complete.")

  list(
    metrics_summary_table = summary_table,
    cell_metrics = univ$cell_results,
    feature_metrics = univ$feat_results,
    bivariate_metrics = bivariate_results,
    zero_probability_curve = sp_m$zero_prob_eval,
    accessibility_sparsity_curve = sp_m$acc_sparsity_eval,
    excess_zero_weights = list(ref = sp_m$ref_excess_w, sim = sp_m$sim_excess_w),
    generative_manifold = list(mmd = sp_m$mmd_eval, precision_recall = sp_m$gen_pr),
    ref_properties = list(cell = props$ref_cell, feature = props$ref_feat),
    sim_properties = list(cell = props$sim_cell, feature = props$sim_feat),
    batch_effects = downstream$batch_effects,
    biological_signal = downstream$biological_signal,
    trajectory = traj_res,
    resource_usage = scal$resource_usage
  )
}

#' Extract Layer-Specific Metadata Vectors from Multiomics Input
#' @keywords internal
#' @noRd
.extract_modality_label <- function(x, n_expected, mod_name, fallback_idx = 1L) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) {
    v <- if (!is.null(x$CellType)) x$CellType else if (!is.null(x$Batch)) x$Batch else if (ncol(x) >= 2L) x[[2L]] else x[[1L]]
    if (length(v) == n_expected) return(as.character(v))
  }
  if (is.list(x) && !is.data.frame(x) && !is.factor(x)) {
    if (!is.null(x[[mod_name]])) return(.extract_modality_label(x[[mod_name]], n_expected, mod_name, fallback_idx))
    alias_names <- if (grepl("atac|peak|chromatin", tolower(mod_name))) c("atac", "peaks", "chromatin") else c("rna", "gene")
    for (anm in alias_names) {
      m_idx <- match(anm, tolower(names(x)))
      if (!is.na(m_idx)) return(.extract_modality_label(x[[m_idx]], n_expected, mod_name, fallback_idx))
    }
    if (length(x) >= fallback_idx) return(.extract_modality_label(x[[fallback_idx]], n_expected, mod_name, fallback_idx))
    return(.extract_modality_label(x[[1L]], n_expected, mod_name, fallback_idx))
  }
  if (length(x) == n_expected) return(as.character(x))
  NULL
}

#' Compute Category III Clustering Concordance across Both Modalities
#' @keywords internal
#' @noRd
.evaluate_multiomics_clustering <- function(ref_multi, sim_multi, ct_m1, ct_m2, bt_m1, bt_m2,
                                            m1_name, m2_name, threads, verbose) {
  if (verbose) message(sprintf("=== [Step 3] Evaluating Clustering Preservation (%s & %s) ===", m1_name, m2_name))
  run_one <- function(sim_o, ref_o, ct, bt, tag) {
    tryCatch(
      evaluate_clustering_metrics(
        data = sim_o, cell_types = ct, ref_data = ref_o,
        ref_cell_types = ct, batch_info = bt, threads = threads
      ),
      error = function(e) {
        if (verbose) message(sprintf("Clustering evaluation note (%s): %s", tag, conditionMessage(e)))
        NULL
      }
    )
  }
  clu_m1 <- run_one(sim_multi[[1L]], ref_multi[[1L]], ct_m1, bt_m1, "M1")
  clu_m2 <- run_one(sim_multi[[2L]], ref_multi[[2L]], ct_m2, bt_m2, "M2")
  fields <- c("silhouette_sim", "cdi", "neighborhood_purity", "ari", "nmi")
  rows <- c(
    .append_rows_from_fields(clu_m1, fields, "Cellular Manifold & Clustering Concordance", "clustering", m1_name),
    .append_rows_from_fields(clu_m2, fields, "Cellular Manifold & Clustering Concordance", "clustering", m2_name)
  )
  list(clustering = list(m1 = clu_m1, m2 = clu_m2), rows = rows)
}

#' Compute Category IV Batch Mixing across Both Modalities
#' @keywords internal
#' @noRd
.evaluate_multiomics_batch <- function(sim_multi, ct_m1, ct_m2, bt_m1, bt_m2,
                                       m1_name, m2_name, threads, verbose) {
  batch_res <- list()
  rows <- list()
  fields <- c("batch_silhouette", "pcr_r2", "cms", "isi", "cross_batch_F1")

  eval_layer <- function(sim_o, bt, ct, mod_nm) {
    if (is.null(bt) || length(unique(bt[!is.na(bt)])) < 2L) return(NULL)
    if (verbose) message(sprintf("=== [Step 4] Evaluating Batch Mixing (%s) ===", mod_nm))
    tryCatch(
      evaluate_batch_metrics(data = sim_o, batch_info = bt, cell_types = ct, threads = threads, verbose = verbose),
      error = function(e) NULL
    )
  }

  res1 <- eval_layer(sim_multi[[1L]], bt_m1, ct_m1, m1_name)
  if (!is.null(res1)) {
    batch_res[[m1_name]] <- res1
    rows <- c(rows, .append_rows_from_fields(res1, fields, "Technical Variation & Batch Mixing Dynamics", "batch_mixing", m1_name))
  }
  res2 <- eval_layer(sim_multi[[2L]], bt_m2, ct_m2, m2_name)
  if (!is.null(res2)) {
    batch_res[[m2_name]] <- res2
    rows <- c(rows, .append_rows_from_fields(res2, fields, "Technical Variation & Batch Mixing Dynamics", "batch_mixing", m2_name))
  }
  list(batch_res = batch_res, rows = rows)
}

#' Compute Category V Differential Signal Fidelity across RNA and ATAC Modalities
#' @keywords internal
#' @noRd
.evaluate_multiomics_signal <- function(ref_multi, sim_multi, ct_m1, ct_m2, bt_m1, bt_m2,
                                        m1_name, m2_name, threads, verbose) {
  signal_res <- list()
  rows <- list()

  if (!is.null(ct_m1) && !is.null(bt_m1)) {
    if (verbose) message(sprintf("=== [Step 5] Evaluating %s Biological Signal Fidelity ===", m1_name))
    sim_ct_m1 <- if (length(ct_m1) == ncol(sim_multi[[1L]])) ct_m1 else rep(ct_m1, length.out = ncol(sim_multi[[1L]]))
    sim_bt_m1 <- if (length(bt_m1) == ncol(sim_multi[[1L]])) bt_m1 else rep(bt_m1, length.out = ncol(sim_multi[[1L]]))
    deg_eval <- tryCatch(
      evaluate_deg_fidelity(
        ref_data = ref_multi[[1L]], sim_data = sim_multi[[1L]],
        ref_celltypes = ct_m1, sim_celltypes = sim_ct_m1,
        batch_ref = bt_m1, batch_sim = sim_bt_m1, threads = threads
      ),
      error = function(e) {
        if (verbose) message("scRNA-seq pseudobulk DE evaluation note: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(deg_eval) && nrow(deg_eval$deg_summary_table) > 0L) {
      signal_res[[m1_name]] <- deg_eval
      rows <- c(rows, .append_rows_from_df(
        deg_eval$deg_summary_table, "Biological Signal & Differential Fidelity",
        "scRNA_differential_expression", m1_name
      ))
    }
  } else if (!is.null(ct_m1) && verbose) {
    message(sprintf("Batch annotations not provided: skipping Category V pseudobulk analysis for %s.", m1_name))
  }

  is_atac_m2 <- .is_atac_modality(ref_multi[[2L]], m2_name)
  if (is_atac_m2) {
    if (!is.null(ct_m2) && !is.null(bt_m2)) {
      if (verbose) message(sprintf("=== [Step 5] Evaluating %s Biological Signal Fidelity ===", m2_name))
      sim_ct_m2 <- if (length(ct_m2) == ncol(sim_multi[[2L]])) ct_m2 else rep(ct_m2, length.out = ncol(sim_multi[[2L]]))
      sim_bt_m2 <- if (length(bt_m2) == ncol(sim_multi[[2L]])) bt_m2 else rep(bt_m2, length.out = ncol(sim_multi[[2L]]))
      da_eval <- tryCatch(
        evaluate_atac_da_fidelity(
          ref_data = ref_multi[[2L]], sim_data = sim_multi[[2L]],
          ref_celltypes = ct_m2, sim_celltypes = sim_ct_m2,
          batch_ref = bt_m2, batch_sim = sim_bt_m2, threads = threads
        ),
        error = function(e) {
          if (verbose) message("scATAC-seq DA evaluation note in multiomics: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(da_eval) && nrow(da_eval$da_summary_table) > 0L) {
        signal_res[[m2_name]] <- da_eval
        rows <- c(rows, .append_rows_from_df(
          da_eval$da_summary_table, "Biological Signal & Differential Fidelity",
          "scATAC_differential_accessibility", m2_name
        ))
      }
    } else if (verbose) {
      message("Batch annotations not provided: skipping Category V pseudobulk DA analysis for scATAC-seq.")
    }
  }

  list(signal_res = signal_res, rows = rows)
}

#' Compute Combined and Separate Multiomics Scalability Metrics (Category VIII)
#' @keywords internal
#' @noRd
.compute_multiomics_scalability <- function(elapsed_time, memory_mb, time_unit, memory_unit,
                                            sim_cores, separate_resources, verbose) {
  t_norm <- .normalize_time_unit(time_unit)
  m_norm <- .normalize_memory_unit(memory_unit)

  parse_sep <- function(spec) {
    if (is.null(spec) || !is.list(spec)) {
      return(list(el_sec = NULL, mem_mib = NULL, cores = NULL, tu = t_norm, mu = m_norm,
                  el_in = NULL, mem_in = NULL))
    }
    tu <- if (!is.null(spec$time_unit)) .normalize_time_unit(spec$time_unit) else t_norm
    mu <- if (!is.null(spec$memory_unit)) .normalize_memory_unit(spec$memory_unit) else m_norm
    list(
      el_sec = .convert_time_to_seconds(spec$elapsed_time, tu),
      mem_mib = .convert_memory_to_mib(spec$memory_mb, mu),
      cores = .normalize_sim_cores(spec$sim_cores),
      tu = tu, mu = mu,
      el_in = spec$elapsed_time, mem_in = spec$memory_mb
    )
  }

  r_sep <- parse_sep(if (is.list(separate_resources)) separate_resources$rna else NULL)
  a_sep <- parse_sep(if (is.list(separate_resources)) separate_resources$atac else NULL)

  has_sep_time <- (!is.null(r_sep$el_sec) && !is.na(r_sep$el_sec)) || (!is.null(a_sep$el_sec) && !is.na(a_sep$el_sec))
  has_sep_mem  <- (!is.null(r_sep$mem_mib) && !is.na(r_sep$mem_mib)) || (!is.null(a_sep$mem_mib) && !is.na(a_sep$mem_mib))

  el_sec <- .convert_time_to_seconds(elapsed_time, t_norm)
  if ((is.null(el_sec) || is.na(el_sec)) && has_sep_time) {
    tp <- stats::na.omit(c(r_sep$el_sec %||% NA_real_, a_sep$el_sec %||% NA_real_))
    if (length(tp) > 0L) el_sec <- sum(tp)
  }

  mem_mib <- .convert_memory_to_mib(memory_mb, m_norm)
  if ((is.null(mem_mib) || is.na(mem_mib)) && has_sep_mem) {
    mp <- stats::na.omit(c(r_sep$mem_mib %||% NA_real_, a_sep$mem_mib %||% NA_real_))
    if (length(mp) > 0L) mem_mib <- max(mp)
  }

  sim_cores_val <- .normalize_sim_cores(sim_cores)
  if (is.null(sim_cores_val) && (!is.null(r_sep$cores) || !is.null(a_sep$cores))) {
    cp <- stats::na.omit(c(r_sep$cores %||% NA_integer_, a_sep$cores %||% NA_integer_))
    if (length(cp) > 0L) sim_cores_val <- max(cp)
  }

  comb_c <- if (!is.null(sim_cores_val)) sim_cores_val else 1L
  c_rna <- if (!is.null(r_sep$cores)) r_sep$cores else comb_c
  c_atac <- if (!is.null(a_sep$cores)) a_sep$cores else comb_c

  cs_rna <- if (!is.null(r_sep$el_sec) && !is.na(r_sep$el_sec)) r_sep$el_sec * c_rna else NULL
  cs_atac <- if (!is.null(a_sep$el_sec) && !is.na(a_sep$el_sec)) a_sep$el_sec * c_atac else NULL
  core_sec_total <- if (has_sep_time) sum(c(cs_rna %||% 0, cs_atac %||% 0))
  else if (!is.null(el_sec) && !is.na(el_sec)) el_sec * comb_c
  else NULL

  mpc_rna <- if (!is.null(r_sep$mem_mib) && !is.na(r_sep$mem_mib)) r_sep$mem_mib / c_rna else NULL
  mpc_atac <- if (!is.null(a_sep$mem_mib) && !is.na(a_sep$mem_mib)) a_sep$mem_mib / c_atac else NULL
  mem_core_total <- if (has_sep_mem) {
    parts <- stats::na.omit(c(mpc_rna %||% NA_real_, mpc_atac %||% NA_real_))
    if (length(parts) > 0L) max(parts) else NULL
  } else if (!is.null(mem_mib) && !is.na(mem_mib)) {
    mem_mib / comb_c
  } else {
    NULL
  }

  rows <- list()
  if (!is.null(core_sec_total) || !is.null(mem_core_total)) {
    if (verbose) message("=== [Step 8] Logging Computational Scalability & Efficiency ===")
    rows <- .append_rows_from_fields(
      list(core_seconds = core_sec_total, memory_per_core_mib = mem_core_total),
      c("core_seconds", "memory_per_core_mib"),
      "(VIII) Computational Scalability & Efficiency", "resource_usage", "Computation"
    )
  }

  list(
    rows = rows,
    sim_cores_val = sim_cores_val,
    resource_usage = list(
      mode = if (has_sep_time || has_sep_mem) "separate" else "combined",
      core_seconds = core_sec_total,
      memory_per_core_mib = mem_core_total,
      memory_per_core_gib = if (!is.null(mem_core_total)) mem_core_total / 1024 else NULL,
      memory_mb = mem_mib, memory_input = memory_mb, memory_unit = m_norm,
      elapsed_time = el_sec, elapsed_time_input = elapsed_time, time_unit = t_norm,
      sim_cores = sim_cores_val,
      rna = list(core_seconds = cs_rna, memory_per_core_mib = mpc_rna,
                 elapsed_time = r_sep$el_sec, elapsed_time_input = r_sep$el_in, time_unit = r_sep$tu,
                 memory_mb = r_sep$mem_mib, memory_input = r_sep$mem_in, memory_unit = r_sep$mu,
                 sim_cores = r_sep$cores),
      atac = list(core_seconds = cs_atac, memory_per_core_mib = mpc_atac,
                  elapsed_time = a_sep$el_sec, elapsed_time_input = a_sep$el_in, time_unit = a_sep$tu,
                  memory_mb = a_sep$mem_mib, memory_input = a_sep$mem_in, memory_unit = a_sep$mu,
                  sim_cores = a_sep$cores)
    )
  )
}

#' Master Multiomics Simulation Benchmarking Pipeline
#'
#' Orchestrates comprehensive benchmarking across all 8 evaluation categories
#' using \code{SingleCellExperiment} inputs (either a multi-modal SCE containing
#' an \code{altExp} for chromatin accessibility, or a named list of two SCE
#' objects such as \code{list(RNA = rna_sce, ATAC = atac_sce)}). Cell-type and
#' batch labels are read directly from \code{colData} via \code{cell_type_col}
#' and \code{batch_col}. Users starting from \code{Seurat} objects should
#' convert via \code{Seurat::as.SingleCellExperiment()}.
#'
#' \if{html}{\figure{workflow_diagram.png}{options: width="800" alt="scSimEval Unified Benchmarking Workflow Diagram"}}
#'
#' @details
#' \strong{Multiomics Metric Measurement Capacity (56 Evaluation Measures Suite):}
#' \itemize{
#'   \item \strong{Paired Multiomics (56 Metrics):} Dual unimodal properties (Cat I--II),
#'     dual clustering (Cat III), dual batch mixing (Cat IV), dual biological signal (Cat V),
#'     cross-modal integration & pairing (Cat VI), scRNA trajectory (Cat VII), and scalability (Cat VIII).
#'   \item \strong{Unpaired Multiomics (51 Metrics):} Evaluates both modalities independently
#'     with population-level integration & alignment, omitting cell-matched pairing metrics.
#' }
#'
#' @param ref_multi Reference multiomics \code{SingleCellExperiment} (with \code{altExp}) or named list of two SCE objects.
#' @param sim_multi Simulated multiomics \code{SingleCellExperiment} (with \code{altExp}) or named list of two SCE objects.
#' @param cell_type_col Column name in \code{colData} containing cell-type labels (default \code{"cell_type"}).
#' @param batch_col Column name in \code{colData} containing batch/replicate labels (default \code{"batch"}).
#' @param pairing Character string specifying \code{"paired"} (default) or \code{"unpaired"} assay structure.
#' @param memory_mb Optional numeric peak memory allocation (in \code{memory_unit}).
#' @param elapsed_time Optional numeric wall-clock runtime (in \code{time_unit}).
#' @param time_unit Unit of \code{elapsed_time}: \code{"seconds"} (default) or \code{"minutes"}.
#' @param memory_unit Unit of \code{memory_mb}: \code{"MiB"} (default) or \code{"MB"}.
#' @param sim_cores Optional positive integer specifying CPU cores used by the simulator.
#' @param separate_resources Optional named list (e.g., \code{list(rna = list(...), atac = list(...))})
#'   giving per-modality \code{elapsed_time}, \code{memory_mb}, \code{time_unit}, \code{memory_unit},
#'   and \code{sim_cores} when RNA and ATAC layers are simulated separately.
#' @param feature_pairs Optional 2-column data.frame of linked feature pairs.
#' @param compute_bivariate Logical, whether to compute 2D bivariate tests. Default \code{FALSE}.
#' @param threads Integer number of CPU threads. Default is \code{1L}.
#' @param trajectory_num_trees Number of Random Forest trees for Category VII \code{wcor_features} (default \code{10000L}).
#' @param verbose Logical, whether to print progress messages. Default is \code{TRUE}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{benchmark_summary_table}{Master tidy data frame with columns \code{Category}, \code{Property}, \code{Metric}, \code{Value}, \code{Modality} (and \code{Sim_Cores} if provided).}
#'   \item{metrics_summary_table}{Alias of \code{benchmark_summary_table}.}
#'   \item{summary_table}{Alias of \code{benchmark_summary_table}.}
#'   \item{pairing}{Resolved pairing mode (\code{"paired"} or \code{"unpaired"}).}
#'   \item{mod1_unimodal}{Full unimodal evaluation output for modality 1.}
#'   \item{mod2_unimodal}{Full unimodal evaluation output for modality 2.}
#'   \item{clustering}{List with elements \code{m1} and \code{m2} containing Category III results.}
#'   \item{batch}{Named list of Category IV batch mixing results per modality.}
#'   \item{signal}{Named list of Category V differential fidelity results per modality.}
#'   \item{trajectory}{Category VII \emph{dynverse} trajectory results inferred from the scRNA-seq layer.}
#'   \item{cross_modality}{List storing Category VI integration output under \code{correlation}.}
#'   \item{resource_usage}{Normalized Category VIII scalability list supporting both \code{"combined"} and \code{"separate"} modes.}
#' }
#' @examples
#' data(example_multiomics, package = "scSimEval")
#' res <- evaluate_multiomics_accuracy(
#'   example_multiomics$ref_multi,
#'   example_multiomics$sim_multi,
#'   verbose = FALSE
#' )
#' head(res$benchmark_summary_table, 10)
#' @export
evaluate_multiomics_accuracy <- function(
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
) {
  pairing <- match.arg(pairing)
  val_ref <- validate_multiomics_input(ref_multi, "ref_multi",
                                       cell_type_col = cell_type_col, batch_col = batch_col)
  val_sim <- validate_multiomics_input(sim_multi, "sim_multi",
                                       cell_type_col = cell_type_col, batch_col = batch_col)

  cell_types <- val_ref$cell_types %||% val_sim$cell_types
  batch_info <- val_ref$batch_info %||% val_sim$batch_info
  m1_name <- val_ref$m1_name %||% "RNA"
  m2_name <- val_ref$m2_name %||% "ATAC"

  ref_list <- stats::setNames(list(val_ref$m1_obj, val_ref$m2_obj), c(m1_name, m2_name))
  sim_list <- stats::setNames(list(val_sim$m1_obj, val_sim$m2_obj), c(m1_name, m2_name))

  c1 <- colnames(ref_list[[1L]])
  c2 <- colnames(ref_list[[2L]])
  if (pairing == "paired" && !is.null(c1) && !is.null(c2) &&
      (length(c1) != length(c2) || !identical(c1, c2))) {
    if (verbose) {
      message(sprintf(
        "scSimEval Note: Modality barcodes differ (%s: %d vs %s: %d). Enforcing unpaired mode.",
        m1_name, length(c1), m2_name, length(c2)
      ))
    }
    pairing <- "unpaired"
  }

  threads <- .set_parallel_threads(threads)
  run_unimodal_layer <- function(idx, mod_nm, step_lbl) {
    if (verbose) message(sprintf("=== [Step %s] Evaluating %s Layer Accuracy ===", step_lbl, mod_nm))
    ev <- evaluate_simulation_accuracy(
      ref_sce = ref_list[[idx]], sim_sce = sim_list[[idx]],
      cell_type_col = cell_type_col, batch_col = batch_col,
      modality = mod_nm, compute_bivariate = compute_bivariate,
      threads = threads, verbose = verbose,
      include_trajectory = FALSE, include_downstream = FALSE
    )
    tbl <- ev$metrics_summary_table
    tbl <- tbl[tbl$Property != "clustering" & !grepl("Trajectory", tbl$Category, ignore.case = TRUE), , drop = FALSE]
    tbl$Modality <- mod_nm
    list(ev = ev, tbl = tbl)
  }

  u1 <- run_unimodal_layer(1L, m1_name, "1")
  u2 <- run_unimodal_layer(2L, m2_name, "2")
  m1_eval <- u1$ev
  m2_eval <- u2$ev
  master_rows <- list(u1$tbl, u2$tbl)

  ct_m1 <- .extract_modality_label(cell_types, ncol(ref_list[[1L]]), m1_name, 1L)
  ct_m2 <- .extract_modality_label(cell_types, ncol(ref_list[[2L]]), m2_name, 2L)
  bt_m1 <- .extract_modality_label(batch_info, ncol(ref_list[[1L]]), m1_name, 1L)
  bt_m2 <- .extract_modality_label(batch_info, ncol(ref_list[[2L]]), m2_name, 2L)

  clu_out <- .evaluate_multiomics_clustering(ref_list, sim_list, ct_m1, ct_m2, bt_m1, bt_m2, m1_name, m2_name, threads, verbose)
  bat_out <- .evaluate_multiomics_batch(sim_list, ct_m1, ct_m2, bt_m1, bt_m2, m1_name, m2_name, threads, verbose)
  sig_out <- .evaluate_multiomics_signal(ref_list, sim_list, ct_m1, ct_m2, bt_m1, bt_m2, m1_name, m2_name, threads, verbose)
  master_rows <- c(master_rows, clu_out$rows, bat_out$rows, sig_out$rows)

  if (!is.null(bat_out$batch_res[[m1_name]])) m1_eval$batch_effects <- bat_out$batch_res[[m1_name]]
  if (!is.null(bat_out$batch_res[[m2_name]])) m2_eval$batch_effects <- bat_out$batch_res[[m2_name]]
  if (!is.null(sig_out$signal_res[[m1_name]])) m1_eval$biological_signal <- sig_out$signal_res[[m1_name]]
  if (!is.null(sig_out$signal_res[[m2_name]])) m2_eval$biological_signal <- sig_out$signal_res[[m2_name]]

  # Category VII: Trajectory & Lineage Topology (scRNA-seq layer)
  m1_is_atac <- .is_atac_modality(ref_list[[1L]], m1_name) || .is_atac_modality(sim_list[[1L]], m1_name)
  m2_is_atac <- .is_atac_modality(ref_list[[2L]], m2_name) || .is_atac_modality(sim_list[[2L]], m2_name)
  rna_idx <- if (!m1_is_atac) 1L else if (!m2_is_atac) 2L else NA_integer_

  traj_eval <- NULL
  if (!is.na(rna_idx)) {
    rna_nm <- if (rna_idx == 1L) m1_name else m2_name
    rna_ct <- if (rna_idx == 1L) ct_m1 else ct_m2
    if (verbose) message(sprintf("=== [Step 6] Evaluating Trajectory & Lineage Topology Inferred from %s ===", rna_nm))
    sim_ct_traj <- if (!is.null(rna_ct)) rep(rna_ct, length.out = ncol(sim_list[[rna_idx]])) else NULL
    traj_eval <- tryCatch(
      evaluate_trajectory_metrics(
        ref_data = ref_list[[rna_idx]], sim_data = sim_list[[rna_idx]],
        cell_types_ref = rna_ct, cell_types_sim = sim_ct_traj,
        num_trees = trajectory_num_trees, threads = threads, verbose = FALSE
      ),
      error = function(e) NULL
    )
    master_rows <- c(master_rows, .append_trajectory_rows(traj_eval, modality = rna_nm))
  }

  # Category VI: Multiomics Integration & Alignment
  if (verbose) message("=== [Step 7] Evaluating Multiomics Integration & Alignment (Category VI) ===")
  int_eval <- tryCatch(
    evaluate_multiomics_integration(
      ref_multi = ref_multi, sim_multi = sim_multi,
      cell_types = if (!is.null(ct_m1) && !is.null(ct_m2)) list(ct_m1, ct_m2) else ct_m1,
      batch_info = if (!is.null(bt_m1) && !is.null(bt_m2)) list(bt_m1, bt_m2) else bt_m1,
      n_pcs = 15L, threads = threads, verbose = verbose
    ),
    error = function(e) {
      if (verbose) message("Cross-modal integration evaluation note: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(int_eval)) {
    master_rows <- c(master_rows, .append_rows_from_df(
      int_eval$integration_summary_table,
      "(VI) Multiomics Integration & Alignment", "multiomics_integration", "Joint (Cross-Modal)"
    ))
  }

  # Category VIII: Scalability & Efficiency
  scal <- .compute_multiomics_scalability(elapsed_time, memory_mb, time_unit, memory_unit,
                                          sim_cores, separate_resources, verbose)
  master_rows <- c(master_rows, scal$rows)

  master_table <- do.call(rbind, master_rows)
  if (!is.null(scal$sim_cores_val) && !is.null(master_table) && nrow(master_table) > 0L) {
    master_table$Sim_Cores <- scal$sim_cores_val
  }
  if (verbose) message("Master multiomics benchmarking complete.")

  list(
    benchmark_summary_table = master_table,
    metrics_summary_table = master_table,
    summary_table = master_table,
    pairing = pairing,
    mod1_unimodal = m1_eval,
    mod2_unimodal = m2_eval,
    clustering = clu_out$clustering,
    batch = bat_out$batch_res,
    signal = sig_out$signal_res,
    trajectory = traj_eval,
    cross_modality = list(
      correlation = int_eval, foscttm = NULL, generation = NULL,
      coupling = NULL, coaccessibility = NULL, coexpression = NULL, profile = NULL
    ),
    resource_usage = scal$resource_usage
  )
}

#' Attach Cell-Type and Batch Metadata to SingleCellExperiment colData When Missing
#' @keywords internal
#' @noRd
.attach_sce_metadata <- function(sce, cell_types = NULL, batch_info = NULL,
                                 cell_type_col = "cell_type", batch_col = "batch") {
  if (!inherits(sce, "SingleCellExperiment")) return(sce)
  cd <- SummarizedExperiment::colData(sce)
  if (!is.null(cell_types) && !cell_type_col %in% colnames(cd) && length(cell_types) == ncol(sce)) {
    cd[[cell_type_col]] <- as.character(cell_types)
  }
  if (!is.null(batch_info) && !batch_col %in% colnames(cd) && length(batch_info) == ncol(sce)) {
    cd[[batch_col]] <- as.character(batch_info)
  }
  SummarizedExperiment::colData(sce) <- cd
  sce
}

#' Coerce Matrix, Path, or SCE Input to SingleCellExperiment
#' @keywords internal
#' @noRd
.ensure_sce <- function(m, auto_sparse = TRUE) {
  if (is.character(m) && length(m) == 1L && file.exists(m)) {
    ext <- tolower(tools::file_ext(m))
    m <- if (ext == "rds") readRDS(m)
    else if (ext == "csv") as.matrix(utils::read.csv(m, row.names = 1L, check.names = FALSE))
    else as.matrix(utils::read.table(m, header = TRUE, row.names = 1L, check.names = FALSE))
  }
  if (inherits(m, "Seurat")) {
    stop("Seurat objects are not directly supported; convert via Seurat::as.SingleCellExperiment().", call. = FALSE)
  }
  if (inherits(m, "SingleCellExperiment")) return(m)
  if (is.data.frame(m)) m <- as.matrix(m)
  if (auto_sparse && !inherits(m, "dgCMatrix") && (is.matrix(m) || inherits(m, "Matrix"))) {
    m <- Matrix::Matrix(m, sparse = TRUE)
  }
  SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
}

#' Parse Multi-Dataset Input Specifications into Normalized Experiment Entries
#' @keywords internal
#' @noRd
.parse_multi_dataset_spec <- function(datasets, pair_by_prefix) {
  parsed <- list()
  infer_mod <- function(nm, explicit = NULL) {
    if (!is.null(explicit)) return(explicit)
    if (grepl("ATAC", nm, ignore.case = TRUE)) "scATAC-seq"
    else if (grepl("RNA", nm, ignore.case = TRUE)) "scRNA-seq"
    else "Single-Cell"
  }
  build_uni <- function(item, nm, pfx) {
    ref_m <- item$ref %||% item$ref_data
    sim_m <- item$sim %||% item$sim_data
    list(
      type = "unimodal", name = pfx, data_name = nm, modality = infer_mod(nm, item$modality),
      ref = ref_m, sim = sim_m, cell_types = item$cell_types, batch_info = item$batch_info,
      memory_mb = item$memory_mb, elapsed_time = item$elapsed_time,
      time_unit = item$time_unit %||% "seconds", memory_unit = item$memory_unit %||% "MiB",
      sim_cores = item$sim_cores %||% item$cores %||% item$n_cores
    )
  }
  build_multi_direct <- function(item, nm, ref_mat, sim_mat) {
    dns <- if (!is.null(names(ref_mat))) paste0(nm, "-", names(ref_mat)) else paste0(nm, c("-Mod1", "-Mod2"))
    sep_res <- if (!is.null(item$elapsed_time_rna) || !is.null(item$elapsed_time_atac) ||
                   !is.null(item$memory_mb_rna) || !is.null(item$memory_mb_atac)) {
      list(
        rna = list(elapsed_time = item$elapsed_time_rna, memory_mb = item$memory_mb_rna,
                   time_unit = item$time_unit_rna %||% "seconds",
                   memory_unit = item$memory_unit_rna %||% "MiB", sim_cores = item$sim_cores_rna),
        atac = list(elapsed_time = item$elapsed_time_atac, memory_mb = item$memory_mb_atac,
                    time_unit = item$time_unit_atac %||% "seconds",
                    memory_unit = item$memory_unit_atac %||% "MiB", sim_cores = item$sim_cores_atac)
      )
    } else {
      item$separate_resources
    }
    list(
      type = "multiomics", name = nm, data_names = dns, ref = ref_mat, sim = sim_mat,
      cell_types = item$cell_types, batch_info = item$batch_info,
      pairing = item$pairing %||% (if (!is.null(item$is_paired)) (if (isTRUE(item$is_paired)) "paired" else "unpaired") else "paired"),
      memory_mb = item$memory_mb, elapsed_time = item$elapsed_time,
      time_unit = item$time_unit %||% "seconds", memory_unit = item$memory_unit %||% "MiB",
      sim_cores = item$sim_cores %||% item$cores %||% item$n_cores,
      separate_resources = sep_res, feature_pairs = item$feature_pairs
    )
  }

  if (!pair_by_prefix) {
    for (d_nm in names(datasets)) {
      item <- datasets[[d_nm]]
      ref_mat <- item$ref %||% item$ref_multi %||% item$ref_data
      sim_mat <- item$sim %||% item$sim_multi %||% item$sim_data
      if (is.list(ref_mat) && length(ref_mat) >= 2L && !is.data.frame(ref_mat)) {
        parsed[[d_nm]] <- build_multi_direct(item, d_nm, ref_mat, sim_mat)
      } else {
        parsed[[d_nm]] <- build_uni(item, d_nm, d_nm)
      }
    }
    return(parsed)
  }

  mod_pat <- "[-_ ]*(scRNA[-_]?seq|scATAC[-_]?seq|scRNA|scATAC|RNA|ATAC|CITE[-_]?seq|CITE|Methyl)"
  prefixes <- sub(paste0(mod_pat, ".*$"), "", names(datasets), ignore.case = TRUE)
  prefixes <- sub("[-_ ]+$", "", prefixes)
  prefixes[prefixes == ""] <- names(datasets)[prefixes == ""]

  for (pfx in unique(prefixes)) {
    m_idx <- which(prefixes == pfx)
    if (length(m_idx) == 1L) {
      item <- datasets[[m_idx]]
      d_nm <- names(datasets)[m_idx]
      ref_mat <- item$ref %||% item$ref_multi %||% item$ref_data
      sim_mat <- item$sim %||% item$sim_multi %||% item$sim_data
      if (is.list(ref_mat) && length(ref_mat) >= 2L && !is.data.frame(ref_mat)) {
        parsed[[pfx]] <- build_multi_direct(item, pfx, ref_mat, sim_mat)
      } else {
        parsed[[d_nm]] <- build_uni(item, d_nm, pfx)
      }
    } else {
      sub_items <- datasets[m_idx]
      sub_names <- names(datasets)[m_idx]
      is_atac <- grepl("ATAC", sub_names, ignore.case = TRUE)
      is_rna <- grepl("RNA", sub_names, ignore.case = TRUE)
      if (any(is_rna) && any(is_atac)) {
        ri <- which(is_rna)[1L]
        ai <- which(is_atac)[1L]
        r_item <- sub_items[[ri]]
        a_item <- sub_items[[ai]]
        parsed[[pfx]] <- list(
          type = "multiomics", name = pfx, data_names = c(sub_names[ri], sub_names[ai]),
          ref = list("scRNA-seq" = r_item$ref %||% r_item$ref_data,
                     "scATAC-seq" = a_item$ref %||% a_item$ref_data),
          sim = list("scRNA-seq" = r_item$sim %||% r_item$sim_data,
                     "scATAC-seq" = a_item$sim %||% a_item$sim_data),
          cell_types = r_item$cell_types %||% a_item$cell_types,
          batch_info = r_item$batch_info %||% a_item$batch_info,
          pairing = r_item$pairing %||% (if (!is.null(r_item$is_paired)) (if (isTRUE(r_item$is_paired)) "paired" else "unpaired") else "paired"),
          separate_resources = list(
            rna = list(elapsed_time = r_item$elapsed_time, memory_mb = r_item$memory_mb,
                       time_unit = r_item$time_unit %||% "seconds",
                       memory_unit = r_item$memory_unit %||% "MiB",
                       sim_cores = r_item$sim_cores %||% r_item$cores %||% r_item$n_cores),
            atac = list(elapsed_time = a_item$elapsed_time, memory_mb = a_item$memory_mb,
                        time_unit = a_item$time_unit %||% "seconds",
                        memory_unit = a_item$memory_unit %||% "MiB",
                        sim_cores = a_item$sim_cores %||% a_item$cores %||% a_item$n_cores)
          ),
          feature_pairs = r_item$feature_pairs
        )
      } else {
        for (idx_k in m_idx) {
          d_nm <- names(datasets)[idx_k]
          parsed[[d_nm]] <- build_uni(datasets[[idx_k]], d_nm, pfx)
        }
      }
    }
  }
  parsed
}

#' Evaluate Multiple Single-Cell Datasets Simultaneously
#'
#' Orchestrates simultaneous benchmarking across a collection of single-cell
#' \code{SingleCellExperiment} datasets, including paired/unpaired multiomics
#' and unimodal scRNA-seq or scATAC-seq experiments.
#'
#' @param datasets A named list of dataset specifications (\code{Format A} grouped multiomics/unimodal lists,
#'   or \code{Format B} flat modality-suffixed lists).
#' @param pair_by_prefix Logical; if \code{TRUE} (default), datasets sharing a root prefix
#'   across RNA and ATAC entries are paired for multiomics evaluation.
#' @param compute_bivariate Logical; whether to compute 2D bivariate tests (default \code{FALSE}).
#' @param threads Integer number of CPU threads. Default is \code{1L}.
#' @param verbose Logical; whether to print progress messages. Default is \code{TRUE}.
#'
#' @return An S3 object of class \code{c("scSimEval_consolidated", "list")} containing:
#' \describe{
#'   \item{consolidated_summary_table}{Master tidy data.frame with columns \code{Dataset}, \code{Data_Name}, \code{Modality}, \code{Category}, \code{Property}, \code{Metric}, and \code{Value}.}
#'   \item{consolidated_score_matrix}{Numeric matrix comparing datasets across all evaluated metric keys.}
#'   \item{dataset_overview}{Data frame summarizing metric counts and mean KS, Wasserstein, and RMSE per dataset.}
#'   \item{dataset_results}{Named list of detailed per-dataset evaluation outputs.}
#' }
#' @examples
#' data(example_scrna, package = "scSimEval")
#' ds <- list("DemoCohort" = list(ref = example_scrna$ref, sim = example_scrna$sim))
#' res <- evaluate_multiple_datasets(ds, verbose = FALSE)
#' res$dataset_overview
#' @export
evaluate_multiple_datasets <- function(
  datasets,
  pair_by_prefix = TRUE,
  compute_bivariate = FALSE,
  threads = 1L,
  verbose = TRUE
) {
  if (!is.list(datasets) || length(datasets) == 0L) {
    stop("`datasets` must be a non-empty named list of datasets.", call. = FALSE)
  }
  if (is.null(names(datasets)) || any(names(datasets) == "")) {
    names(datasets) <- paste0("Dataset_", seq_along(datasets))
  }

  parsed <- .parse_multi_dataset_spec(datasets, pair_by_prefix)
  all_master_rows <- list()
  dataset_results <- list()
  n_exp <- length(parsed)

  if (verbose) message(sprintf(">>> Initializing scSimEval benchmarking for %d experiment(s)...", n_exp))

  for (i in seq_along(parsed)) {
    exp <- parsed[[i]]
    exp_name <- exp$name
    if (verbose) message(sprintf("\n[%d/%d] Benchmarking dataset: '%s' (%s)...", i, n_exp, exp_name, exp$type))

    if (exp$type == "multiomics") {
      if (is.list(exp$ref) && !inherits(exp$ref, "SingleCellExperiment")) {
        exp$ref <- lapply(exp$ref, function(x) .attach_sce_metadata(x, exp$cell_types, exp$batch_info))
        exp$sim <- lapply(exp$sim, function(x) .attach_sce_metadata(x, exp$cell_types, exp$batch_info))
      }
      res <- evaluate_multiomics_accuracy(
        ref_multi = exp$ref, sim_multi = exp$sim,
        pairing = exp$pairing %||% "paired",
        memory_mb = exp$memory_mb, elapsed_time = exp$elapsed_time,
        time_unit = exp$time_unit %||% "seconds", memory_unit = exp$memory_unit %||% "MiB",
        sim_cores = exp$sim_cores, separate_resources = exp$separate_resources,
        feature_pairs = exp$feature_pairs, compute_bivariate = compute_bivariate,
        threads = threads, verbose = verbose
      )
      df <- res$benchmark_summary_table
      df$Dataset <- exp_name
      m1_tag <- names(exp$ref)[1L] %||% "scRNA-seq"
      m2_tag <- names(exp$ref)[2L] %||% "scATAC-seq"
      dn1 <- if (length(exp$data_names) >= 1L) exp$data_names[1L] else paste0(exp_name, "-", m1_tag)
      dn2 <- if (length(exp$data_names) >= 2L) exp$data_names[2L] else paste0(exp_name, "-", m2_tag)
      df$Data_Name <- vapply(df$Modality, function(mod_r) {
        if (mod_r == m1_tag || grepl("RNA", mod_r, ignore.case = TRUE)) dn1
        else if (mod_r == m2_tag || grepl("ATAC", mod_r, ignore.case = TRUE)) dn2
        else if (grepl("Cross|Joint", mod_r, ignore.case = TRUE)) paste0(exp_name, "-Joint")
        else paste0(exp_name, "-", mod_r)
      }, character(1))
      res_key <- exp_name
    } else {
      ref_sce <- .attach_sce_metadata(exp$ref, exp$cell_types, exp$batch_info)
      sim_sce <- .attach_sce_metadata(exp$sim, exp$cell_types, exp$batch_info)
      res <- evaluate_simulation_accuracy(
        ref_sce = ref_sce, sim_sce = sim_sce, modality = exp$modality,
        compute_bivariate = compute_bivariate, threads = threads,
        memory_mb = exp$memory_mb, elapsed_time = exp$elapsed_time,
        time_unit = exp$time_unit %||% "seconds", memory_unit = exp$memory_unit %||% "MiB",
        sim_cores = exp$sim_cores, verbose = verbose
      )
      df <- res$metrics_summary_table
      df$Dataset <- exp_name
      df$Data_Name <- exp$data_name
      df$Modality <- exp$modality
      res_key <- exp$data_name
    }

    keep_cols <- intersect(c("Dataset", "Data_Name", "Modality", "Category", "Property", "Metric", "Value", "Sim_Cores"), colnames(df))
    all_master_rows[[length(all_master_rows) + 1L]] <- df[, keep_cols, drop = FALSE]
    dataset_results[[res_key]] <- res
  }

  if (any(vapply(all_master_rows, function(d) "Sim_Cores" %in% colnames(d), logical(1)))) {
    all_master_rows <- lapply(all_master_rows, function(d) {
      if (!"Sim_Cores" %in% colnames(d)) d$Sim_Cores <- NA_integer_
      d
    })
  }
  consolidated_table <- do.call(rbind, all_master_rows)
  rownames(consolidated_table) <- NULL

  metric_keys <- paste0(consolidated_table$Category, " | ", consolidated_table$Property, " [", consolidated_table$Metric, "]")
  u_keys <- unique(metric_keys)
  u_dns <- unique(consolidated_table$Data_Name)
  score_mat <- matrix(NA_real_, nrow = length(u_keys), ncol = length(u_dns), dimnames = list(u_keys, u_dns))
  for (r in seq_len(nrow(consolidated_table))) {
    score_mat[metric_keys[r], consolidated_table$Data_Name[r]] <- consolidated_table$Value[r]
  }

  dataset_overview <- do.call(rbind, lapply(u_dns, function(dn) {
    sub_df <- consolidated_table[consolidated_table$Data_Name == dn, ]
    ks_vals <- sub_df$Value[sub_df$Metric == "KS"]
    w_vals <- sub_df$Value[sub_df$Metric == "Wasserstein"]
    r_vals <- sub_df$Value[sub_df$Metric == "RMSE"]
    data.frame(
      Dataset = sub_df$Dataset[1L], Data_Name = dn, Modality = sub_df$Modality[1L],
      Total_Metrics = nrow(sub_df),
      Mean_KS_Distance = if (length(ks_vals) > 0L) round(mean(ks_vals, na.rm = TRUE), 4) else NA_real_,
      Mean_Wasserstein = if (length(w_vals) > 0L) round(mean(w_vals, na.rm = TRUE), 4) else NA_real_,
      Mean_RMSE = if (length(r_vals) > 0L) round(mean(r_vals, na.rm = TRUE), 4) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  out <- list(
    consolidated_summary_table = consolidated_table,
    consolidated_score_matrix = score_mat,
    dataset_overview = dataset_overview,
    dataset_results = dataset_results
  )
  class(out) <- c("scSimEval_consolidated", "list")
  out
}

#' @export
print.scSimEval_consolidated <- function(x, ...) {
  cat("=====================================================================\n")
  cat("  scSimEval Consolidated Multi-Dataset Benchmark Results\n")
  cat("=====================================================================\n")
  cat(sprintf("Total Evaluated Metrics : %d across %d dataset entries\n",
              nrow(x$consolidated_summary_table),
              length(unique(x$consolidated_summary_table$Data_Name))))
  cat(sprintf("Unique Datasets        : %s\n", paste(unique(x$consolidated_summary_table$Dataset), collapse = ", ")))
  cat(sprintf("Unique Modalities      : %s\n", paste(unique(x$consolidated_summary_table$Modality), collapse = ", ")))
  cat("\nDataset Overview Summary:\n")
  print(x$dataset_overview, row.names = FALSE)
  cat("=====================================================================\n")
  invisible(x)
}

#' Aggregate Suite Scores and Cross-Dataset Summary Statistics
#' @keywords internal
#' @noRd
.aggregate_suite_stats <- function(granular_tbl) {
  if ("Category" %in% colnames(granular_tbl)) {
    mc_override <- .METRIC_CATEGORY_MAP[granular_tbl$Metric]
    idx <- !is.na(mc_override)
    granular_tbl$Category[idx] <- mc_override[idx]
  }
  if (nrow(granular_tbl) > 0L && "Score_Raw" %in% colnames(granular_tbl)) {
    granular_tbl <- normalize_bubble_scores(granular_tbl)
    granular_tbl$Score <- granular_tbl$Score_Norm
  }
  agg_norm <- stats::aggregate(Score_Norm ~ Method + Category + Metric, data = granular_tbl, FUN = mean, na.rm = TRUE)
  agg_raw  <- stats::aggregate(Score_Raw ~ Method + Category + Metric, data = granular_tbl, FUN = mean, na.rm = TRUE)
  summary_tbl <- merge(agg_norm, agg_raw, by = c("Method", "Category", "Metric"))
  summary_tbl$Score <- round(summary_tbl$Score_Norm, 4)
  summary_tbl$Score_Norm <- round(summary_tbl$Score_Norm, 4)
  summary_tbl$Score_Raw <- round(summary_tbl$Score_Raw, 4)
  summary_tbl$Simulator <- summary_tbl$Method

  stats_agg <- stats::aggregate(
    Score_Norm ~ Method + Category + Metric, data = granular_tbl,
    FUN = function(x) c(
      Mean = mean(x, na.rm = TRUE), SD = stats::sd(x, na.rm = TRUE),
      SE = stats::sd(x, na.rm = TRUE) / sqrt(length(stats::na.omit(x))),
      Min = min(x, na.rm = TRUE), Max = max(x, na.rm = TRUE),
      N = length(stats::na.omit(x))
    )
  )
  stats_df <- data.frame(
    Method = stats_agg$Method, Category = stats_agg$Category, Metric = stats_agg$Metric,
    Mean_Score = round(stats_agg$Score_Norm[, "Mean"], 4),
    SD_Score = round(stats_agg$Score_Norm[, "SD"], 4),
    SE_Score = round(stats_agg$Score_Norm[, "SE"], 4),
    Min_Score = round(stats_agg$Score_Norm[, "Min"], 4),
    Max_Score = round(stats_agg$Score_Norm[, "Max"], 4),
    N_Datasets = as.integer(stats_agg$Score_Norm[, "N"]),
    stringsAsFactors = FALSE
  )
  list(granular_tbl = granular_tbl, summary_tbl = summary_tbl, stats_df = stats_df)
}

#' Benchmark Simulation Methods Across Multiple Single-Cell Datasets
#'
#' Unified benchmarking runner for multi-dataset and multi-simulator studies.
#' Manages sparse matrix coercion, \code{SingleCellExperiment} conversion,
#' parallel metric computation, and produces both granular and cross-dataset
#' aggregated summary tables ready for visualization.
#'
#' @param datasets Named list of empirical reference datasets (\code{SingleCellExperiment},
#'   count matrix, file path, or list containing \code{ref}, \code{cell_types}, \code{batch_info}).
#' @param simulators Named list of simulator outputs corresponding to \code{datasets}.
#' @param compute_bivariate Logical; compute 2D bivariate tests. Default \code{FALSE}.
#' @param threads Integer number of parallel workers. Default \code{1L}.
#' @param auto_sparse Logical; coerce dense matrices to \code{dgCMatrix}. Default \code{TRUE}.
#' @param run_gc Logical; run \code{gc()} between datasets. Default \code{TRUE}.
#' @param save_path Optional \code{.rds} path to save the output bundle.
#' @param verbose Logical; print progress logging. Default \code{TRUE}.
#'
#' @return An S3 object of class \code{c("scSimEval_multidataset_benchmark", "scSimEval_consolidated", "list")} containing:
#' \describe{
#'   \item{granular_benchmark_table}{Tidy data.frame across all dataset-simulator evaluations.}
#'   \item{benchmark_summary_table}{Cross-dataset averaged fidelity scores formatted for \code{\link{plot_benchmark_bubble_matrix}}.}
#'   \item{cross_dataset_stats}{Mean, SD, SE, Min, Max, and \code{N_Datasets} for each simulator and metric.}
#'   \item{dataset_summary}{Summary table of reference dimensions, sparsity, and cell-type counts.}
#'   \item{dataset_summary_df}{Alias of \code{dataset_summary}.}
#'   \item{methods}{Character vector of evaluated simulator names.}
#'   \item{datasets}{Character vector of evaluated dataset names.}
#'   \item{is_multidataset}{Logical \code{TRUE}.}
#'   \item{n_datasets}{Integer count of datasets.}
#'   \item{n_simulators}{Integer count of simulators.}
#' }
#' @export
#' @examples
#' data(example_scrna, package = "scSimEval")
#' refs <- list("CohortA" = example_scrna$ref)
#' sims <- list("Splatter" = list("CohortA" = example_scrna$sim))
#' res <- benchmark_simulation_suite(refs, sims, verbose = FALSE)
#' print(res)
#' head(res$benchmark_summary_table, 8)
benchmark_simulation_suite <- function(
  datasets,
  simulators,
  compute_bivariate = FALSE,
  threads           = 1L,
  auto_sparse       = TRUE,
  run_gc            = TRUE,
  save_path         = NULL,
  verbose           = TRUE
) {
  if (!is.list(datasets) || length(datasets) == 0L) stop("`datasets` must be a non-empty named list.", call. = FALSE)
  if (!is.list(simulators) || length(simulators) == 0L) stop("`simulators` must be a non-empty named list.", call. = FALSE)
  if (is.null(names(datasets)) || any(names(datasets) == "")) names(datasets) <- paste0("Dataset_", seq_along(datasets))
  if (is.null(names(simulators)) || any(names(simulators) == "")) names(simulators) <- paste0("Simulator_", seq_along(simulators))

  d_names <- names(datasets)
  s_names <- names(simulators)
  total_evals <- length(d_names) * length(s_names)
  if (verbose) {
    message(sprintf("scSimEval Multi-Dataset Suite: %d datasets x %d simulators (%d workers)",
                    length(d_names), length(s_names), threads))
  }

  all_granular_rows <- list()
  dataset_summaries <- list()
  eval_count <- 0L

  for (i in seq_along(d_names)) {
    d_name <- d_names[i]
    d_item <- datasets[[d_name]]
    ref_raw <- if (is.list(d_item) && !inherits(d_item, "Matrix") && !is.data.frame(d_item) && !is.matrix(d_item)) {
      d_item$ref %||% d_item[[1L]]
    } else {
      d_item
    }
    ct_labels <- if (is.list(d_item)) d_item$cell_types else NULL
    bt_labels <- if (is.list(d_item)) d_item$batch_info else NULL
    ref_sce <- .attach_sce_metadata(.ensure_sce(ref_raw, auto_sparse), ct_labels, bt_labels)

    cnts_ref <- SingleCellExperiment::counts(ref_sce)
    n_c <- ncol(cnts_ref)
    n_f <- nrow(cnts_ref)
    sparsity_val <- round(1 - (as.numeric(Matrix::nnzero(cnts_ref)) / (as.numeric(n_c) * as.numeric(n_f))), 4)
    dataset_summaries[[d_name]] <- data.frame(
      "Dataset / Simulator" = d_name, "Role" = "Biological Reference",
      "Modality" = "Single-Cell (Counts)", "Cells (N)" = n_c, "Features (p)" = n_f,
      "Sparsity (%)" = sprintf("%.1f%%", sparsity_val * 100),
      "Cell Types" = if (!is.null(ct_labels)) length(unique(ct_labels)) else NA_integer_,
      check.names = FALSE, stringsAsFactors = FALSE
    )

    for (s_name in s_names) {
      eval_count <- eval_count + 1L
      sc <- simulators[[s_name]]
      sim_raw <- if (is.list(sc) && d_name %in% names(sc)) sc[[d_name]] else if (is.list(sc) && length(sc) >= i) sc[[i]] else if (is.character(sc) && length(sc) >= i) sc[i] else sc
      if (is.null(sim_raw)) {
        warning(sprintf("Simulator '%s' has no simulation output for '%s'. Skipping.", s_name, d_name))
        next
      }
      sim_sce <- .attach_sce_metadata(.ensure_sce(sim_raw, auto_sparse), ct_labels, bt_labels)
      if (verbose) message(sprintf("  [%d/%d] Evaluating '%s' on '%s'...", eval_count, total_evals, s_name, d_name))

      res_ij <- tryCatch(
        evaluate_simulation_accuracy(ref_sce = ref_sce, sim_sce = sim_sce,
                                     compute_bivariate = compute_bivariate, threads = threads, verbose = FALSE),
        error = function(e) {
          warning(sprintf("Failed evaluating '%s' on '%s': %s", s_name, d_name, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(res_ij)) {
        tbl_ij <- res_ij$metrics_summary_table
        tbl_ij$Dataset <- d_name
        tbl_ij$Method <- s_name
        tbl_ij$Simulator <- s_name
        if (!"Score_Raw" %in% colnames(tbl_ij) && "Value" %in% colnames(tbl_ij)) tbl_ij$Score_Raw <- tbl_ij$Value
        if (!"Score" %in% colnames(tbl_ij)) tbl_ij$Score <- tbl_ij$Score_Raw
        if (!"Score_Norm" %in% colnames(tbl_ij)) tbl_ij <- normalize_bubble_scores(tbl_ij)
        all_granular_rows[[length(all_granular_rows) + 1L]] <- tbl_ij
      }
    }
    if (run_gc) gc(verbose = FALSE)
  }

  if (length(all_granular_rows) == 0L) stop("No evaluations succeeded across the suite.", call. = FALSE)
  agg <- .aggregate_suite_stats(do.call(rbind, all_granular_rows))
  d_summary_df <- do.call(rbind, dataset_summaries)

  out <- list(
    granular_benchmark_table = agg$granular_tbl,
    benchmark_summary_table = agg$summary_tbl,
    cross_dataset_stats = agg$stats_df,
    dataset_summary = d_summary_df,
    dataset_summary_df = d_summary_df,
    methods = unique(agg$summary_tbl$Method),
    datasets = d_names,
    is_multidataset = TRUE,
    n_datasets = length(d_names),
    n_simulators = length(s_names)
  )
  class(out) <- c("scSimEval_multidataset_benchmark", "scSimEval_consolidated", "list")
  if (!is.null(save_path)) saveRDS(out, save_path)
  out
}

#' @export
print.scSimEval_multidataset_benchmark <- function(x, ...) {
  cat("=====================================================================\n")
  cat("  scSimEval Multi-Dataset Benchmark Results\n")
  cat("=====================================================================\n")
  cat(sprintf("Datasets Evaluated     : %d (%s)\n", length(x$datasets), paste(utils::head(x$datasets, 4), collapse = ", ")))
  cat(sprintf("Simulators Evaluated   : %d (%s)\n", length(x$methods), paste(x$methods, collapse = ", ")))
  cat(sprintf("Total Evaluated Scores : %d across all datasets & metrics\n", nrow(x$granular_benchmark_table)))
  cat(sprintf("Cross-Dataset Average  : %d averaged metric entries in summary table\n", nrow(x$benchmark_summary_table)))
  cat("=====================================================================\n")
  invisible(x)
}
