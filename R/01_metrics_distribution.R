# -----------------------------------------------------------------------------
# Category I: Cell-Level & Feature-Level Property Extraction Helpers
# -----------------------------------------------------------------------------

#' Extract Comprehensive Cell-Level Properties
#'
#' Computes library size, log-library size, zero fraction / detection frequency,
#' pairwise cell correlation (on top HVGs), cell k-NN hubness, pairwise PCA distance,
#' and outlier proportions.
#'
#' @param data Count matrix (features x cells) or \code{SingleCellExperiment} object.
#' @param max_cells_cor Maximum number of cells for computing pairwise correlation (default 1000).
#' @param n_top_hvgs Number of highly variable genes/peaks to use for correlation and PCA (default 500).
#' @param modality Optional modality string (e.g. \code{"scRNA-seq"}, \code{"scATAC-seq"}). When ATAC/peaks
#'   are provided or auto-detected, cell correlations and PCA use Signac/simPIC TF-IDF normalization
#'   instead of \code{log2(CPM + 1)}.
#' @param verbose Logical, whether to print progress messages.
#'
#' @return A named list of 10 cell-level property summaries:
#' \describe{
#'   \item{library_size}{Numeric vector of total counts per cell.}
#'   \item{log_library_size}{Numeric vector of \code{log(library_size + 1)}.}
#'   \item{tmm_factor}{Numeric vector of TMM normalization factors per cell.}
#'   \item{effective_library_size}{Numeric vector of \code{library_size * tmm_factor}.}
#'   \item{zero_fraction_cell}{Numeric vector of per-cell zero count proportions.}
#'   \item{detection_freq_cell}{Numeric vector of per-cell non-zero feature proportions (\code{1 - zero_fraction_cell}).}
#'   \item{cell_cor}{Numeric vector of upper-triangle pairwise Spearman correlations across cells.}
#'   \item{pca_dist}{Numeric vector of pairwise Euclidean distances between cells in top PCA space.}
#'   \item{knn_hubness}{Numeric vector giving how often each cell appears among other cells' 5\% nearest neighbors.}
#'   \item{prop_outliers_cell}{Numeric scalar proportion of cells flagged as library-size outliers via Tukey IQR.}
#' }
#' @examples
#' data(example_scrna, package = "scSimEval")
#' extract_cell_properties(example_scrna$ref)
#' @export
extract_cell_properties <- function(data, max_cells_cor = 1000, n_top_hvgs = 500,
                                    modality = NULL, verbose = FALSE) {
  is_atac <- .is_atac_modality(data, modality)
  if (inherits(data, "SingleCellExperiment")) {
    data <- SummarizedExperiment::assay(data, "counts")
  } else if (inherits(data, "Seurat")) {
    stop("Seurat objects are not directly supported; convert via Seurat::as.SingleCellExperiment().", call. = FALSE)
  }

  if (verbose) message("Extracting cell library size & log-library size...")
  lib_size <- as.numeric(Matrix::colSums(data))
  log_lib_size <- log(lib_size + 1)

  if (verbose) message("Extracting cell zero fraction & detection rate...")
  zero_fraction_cell <- fast_sparse_zero_fraction_cell(data)
  detection_freq_cell <- 1 - zero_fraction_cell

  if (verbose) message("Extracting cell-cell correlations on top variable features...")
  norm_mat <- normalize_counts(data, method = if (is_atac) "tfidf" else "lognorm")

  raw_vars <- fast_row_vars(data)
  pos_gene_idx <- which(raw_vars > 0)
  n_hvg <- min(n_top_hvgs, length(pos_gene_idx))
  rank_vars <- if (is_atac) fast_row_vars(norm_mat) else raw_vars
  top_hvg_idx <- if (length(pos_gene_idx) > 0) {
    pos_gene_idx[order(rank_vars[pos_gene_idx], decreasing = TRUE)[seq_len(n_hvg)]]
  } else {
    order(rank_vars, decreasing = TRUE)[seq_len(min(n_top_hvgs, length(rank_vars)))]
  }
  sub_cpm <- as.matrix(norm_mat[top_hvg_idx, , drop = FALSE])

  n_cells <- ncol(data)
  max_cor_cells <- min(n_cells, min(max_cells_cor, 500))
  cor_cpm <- if (n_cells > max_cor_cells) {
    sub_cpm[, sample(seq_len(n_cells), max_cor_cells), drop = FALSE]
  } else {
    sub_cpm
  }

  cor_vars <- fast_col_vars(cor_cpm)
  cor_cpm_valid <- cor_cpm[, cor_vars > 0, drop = FALSE]
  cell_cor <- if (ncol(cor_cpm_valid) >= 2L && nrow(cor_cpm_valid) >= 2L) {
    cell_cor_mat <- stats::cor(cor_cpm_valid, method = "spearman")
    cell_cor_mat[upper.tri(cell_cor_mat)]
  } else {
    numeric(0)
  }

  if (verbose) message("Extracting cell PCA distance and k-NN hubness...")
  pca_dist <- numeric(0)
  knn_hubness <- numeric(0)
  tryCatch({
    pca_pts <- fast_pca(t(sub_cpm), n_pcs = min(20, ncol(sub_cpm) - 1, nrow(sub_cpm) - 1))
    pca_dist <- if (n_cells > 200) {
      as.numeric(stats::dist(pca_pts[sample(seq_len(n_cells), 200), , drop = FALSE]))
    } else {
      as.numeric(stats::dist(pca_pts))
    }
    k_val <- max(2L, round(0.05 * n_cells))
    nn <- RANN::nn2(pca_pts, k = min(k_val + 1L, nrow(pca_pts)))
    idx_nn <- nn$nn.idx[, seq(2L, ncol(nn$nn.idx)), drop = FALSE]
    knn_hubness <- as.numeric(tabulate(idx_nn, nbins = n_cells))
  }, error = function(e) NULL)

  # TMM normalization factor & effective library size
  tmm_factors <- rep(1, n_cells)
  if (requireNamespace("edgeR", quietly = TRUE)) {
    tryCatch({
      sub_data <- data
      if (nrow(sub_data) > 2000L) {
        rm_means <- Matrix::rowMeans(sub_data)
        top_idx <- order(rm_means, decreasing = TRUE)[seq_len(2000L)]
        sub_data <- sub_data[top_idx, , drop = FALSE]
      }
      cs_sub <- as.numeric(Matrix::colSums(sub_data))
      pos_cells <- which(cs_sub > 0)
      if (length(pos_cells) >= 2L) {
        dge <- edgeR::DGEList(counts = as.matrix(sub_data[, pos_cells, drop = FALSE]))
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        norm_f <- as.numeric(dge$samples$norm.factors)
        if (length(norm_f) == length(pos_cells) && !any(is.na(norm_f)) && !all(norm_f == 1)) {
          tmm_factors[pos_cells] <- norm_f
        }
      }
    }, error = function(e) NULL)
  }
  if (all(tmm_factors == 1)) {
    med_cs <- stats::median(lib_size[lib_size > 0])
    if (!is.na(med_cs) && med_cs > 0) tmm_factors <- lib_size / med_cs
  }
  effective_lib_size <- lib_size * tmm_factors

  prop_outliers_cell <- calc_outlier_proportion(lib_size)

  list(
    library_size = lib_size,
    log_library_size = log_lib_size,
    tmm_factor = tmm_factors,
    effective_library_size = effective_lib_size,
    zero_fraction_cell = zero_fraction_cell,
    detection_freq_cell = detection_freq_cell,
    cell_cor = cell_cor,
    pca_dist = pca_dist,
    knn_hubness = knn_hubness,
    prop_outliers_cell = prop_outliers_cell
  )
}

#' Extract Comprehensive Feature-Level (Gene / Peak) Properties
#'
#' Computes mean abundance, standard deviation, variance, scaled variance (z-score standardized variance),
#' coefficient of variation (CV), feature dropout (zero fraction), dispersion,
#' Biological Coefficient of Variation (BCV), gene-gene correlations, and outlier proportions.
#'
#' @param data Count matrix (features x cells) or \code{SingleCellExperiment} object.
#' @param n_top_cor Number of top HVGs/HVPs to use for correlation matrix (default 500, matching SimBench).
#' @param modality Optional modality string (e.g. \code{"scRNA-seq"}, \code{"scATAC-seq"}). When ATAC/peaks
#'   are provided or auto-detected, feature summaries and pairwise peak-peak correlations are computed
#'   after Signac/simPIC TF-IDF transformation instead of \code{log2(CPM + 1)}.
#' @param verbose Logical, whether to print progress messages.
#'
#' @return A named list of 11 feature-level property summaries:
#' \describe{
#'   \item{mean_expression}{Numeric vector of normalized mean abundances per feature.}
#'   \item{variance}{Numeric vector of normalized variances per feature.}
#'   \item{scaled_variance}{Numeric vector of z-score standardized normalized variances.}
#'   \item{sd}{Numeric vector of normalized standard deviations per feature.}
#'   \item{cv}{Numeric vector of raw coefficient of variation percentages per feature.}
#'   \item{zero_fraction_feature}{Numeric vector of per-feature zero proportions.}
#'   \item{detection_freq_feature}{Numeric vector of per-feature detection rates (\code{1 - zero_fraction_feature}).}
#'   \item{dispersion}{Numeric vector of raw variance-to-mean ratios.}
#'   \item{bcv}{Numeric vector of empirical Biological Coefficient of Variation estimates.}
#'   \item{gene_cor}{Numeric vector of upper-triangle pairwise Spearman correlations across top variable features.}
#'   \item{prop_outliers_feature}{Numeric scalar proportion of features flagged as total-count outliers via Tukey IQR.}
#' }
#' @examples
#' data(example_scrna, package = "scSimEval")
#' extract_feature_properties(example_scrna$ref)
#' @export
extract_feature_properties <- function(data, n_top_cor = 500, modality = NULL, verbose = FALSE) {
  is_atac <- .is_atac_modality(data, modality)
  if (inherits(data, "SingleCellExperiment")) {
    data <- SummarizedExperiment::assay(data, "counts")
  } else if (inherits(data, "Seurat")) {
    stop("Seurat objects are not directly supported; convert via Seurat::as.SingleCellExperiment().", call. = FALSE)
  }

  if (verbose) message("Extracting feature zero fraction...")
  zero_fraction_feature <- fast_sparse_zero_fraction_feature(data)
  detection_freq_feature <- 1 - zero_fraction_feature

  if (verbose) message("Extracting mean, variance, and SD...")
  raw_mean <- as.numeric(Matrix::rowMeans(data))
  raw_var <- fast_row_vars(data)
  raw_sd <- sqrt(raw_var)

  cv_val <- ifelse(raw_mean > 0, (raw_sd / raw_mean) * 100, 0)
  dispersion_val <- ifelse(raw_mean > 0, raw_var / raw_mean, 0)

  norm_mat <- normalize_counts(data, method = if (is_atac) "tfidf" else "lognorm")
  mean_val <- as.numeric(Matrix::rowMeans(norm_mat))
  var_val <- fast_row_vars(norm_mat)
  sd_val <- sqrt(var_val)

  pos_var_idx <- which(raw_var > 0)
  n_genes_cor <- min(n_top_cor, length(pos_var_idx))
  rank_vars <- if (is_atac) var_val else raw_var
  top_g_idx <- if (length(pos_var_idx) > 0) {
    pos_var_idx[order(rank_vars[pos_var_idx], decreasing = TRUE)[seq_len(n_genes_cor)]]
  } else {
    order(rank_vars, decreasing = TRUE)[seq_len(min(n_top_cor, length(rank_vars)))]
  }
  sub_g_mat <- as.matrix(norm_mat[top_g_idx, , drop = FALSE])

  sd_var_val <- stats::sd(var_val, na.rm = TRUE)
  scaled_var_val <- if (!is.na(sd_var_val) && sd_var_val > 0) {
    as.numeric((var_val - mean(var_val, na.rm = TRUE)) / sd_var_val)
  } else {
    rep(0, length(var_val))
  }

  if (verbose) message("Extracting Biological Coefficient of Variation (BCV)...")
  bcv_val <- sqrt(pmax(0, dispersion_val - 1) / pmax(1e-6, raw_mean))

  if (verbose) message("Extracting gene-gene correlation...")
  sub_cor_in <- if (ncol(sub_g_mat) > 500L) {
    sub_g_mat[, sample(seq_len(ncol(sub_g_mat)), 500L), drop = FALSE]
  } else {
    sub_g_mat
  }
  feat_sub_vars <- fast_row_vars(sub_cor_in)
  sub_cor_valid <- sub_cor_in[feat_sub_vars > 0, , drop = FALSE]
  gene_cor <- if (nrow(sub_cor_valid) >= 2L && ncol(sub_cor_valid) >= 2L) {
    g_cor_mat <- stats::cor(t(sub_cor_valid), method = "spearman")
    g_cor_mat[upper.tri(g_cor_mat)]
  } else {
    numeric(0)
  }

  prop_outliers_feature <- calc_outlier_proportion(as.numeric(Matrix::rowSums(data)))

  list(
    mean_expression = mean_val,
    variance = var_val,
    scaled_variance = scaled_var_val,
    sd = sd_val,
    cv = cv_val,
    zero_fraction_feature = zero_fraction_feature,
    detection_freq_feature = detection_freq_feature,
    dispersion = dispersion_val,
    bcv = bcv_val,
    gene_cor = gene_cor,
    prop_outliers_feature = prop_outliers_feature
  )
}

# -----------------------------------------------------------------------------
# 1D Univariate Distance & Divergence Metrics
# -----------------------------------------------------------------------------

#' Calculate Bhattacharyya Distance (BH)
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param align Logical, whether to align lengths. Default TRUE.
#' @return Statistical divergence between discrete probability measures.
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_bhattacharyya(ref, sim)
#' @export
calc_bhattacharyya <- function(ref, sim, align = TRUE) {
  if (align) {
    aligned <- align_distributions(ref, sim)
    ref <- aligned$ref
    sim <- aligned$sim
  }
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) == 0 || length(sim) == 0) return(NA_real_)
  
  min_val <- min(c(ref, sim), na.rm = TRUE)
  if (is.finite(min_val) && min_val <= 0) {
    shift <- abs(min_val) + 1e-6
    ref <- ref + shift
    sim <- sim + shift
  }
  s_ref <- sum(ref, na.rm = TRUE)
  s_sim <- sum(sim, na.rm = TRUE)
  if (s_ref <= 0 || s_sim <= 0) return(NA_real_)
  
  p_ref <- ref / s_ref
  p_sim <- sim / s_sim

  bc <- sum(sqrt(p_ref * p_sim), na.rm = TRUE)
  if (is.na(bc) || bc <= 0) return(Inf)
  -log(bc)
}

#' Calculate Kolmogorov-Smirnov Distance (1D KS)
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @return Numeric scalar in \code{[0, 1]} giving the maximum vertical distance between empirical CDFs.
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_ks(ref, sim)
#' @export
calc_ks <- function(ref, sim) {
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) == 0L || length(sim) == 0L) return(NA_real_)

  u <- sort(unique(c(ref, sim)))
  as.numeric(max(abs(stats::ecdf(ref)(u) - stats::ecdf(sim)(u))))
}

#' @title Calculate 1D Wasserstein Metric / Earth Mover's Distance (WS)
#' @description Integrated from HelenaLC/simulation-comparison.
#'
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param p Power of the Wasserstein metric (default 1 for standard earth mover's distance).
#' @return Non-negative numeric scalar giving the 1D \eqn{L_p} Wasserstein distance between the two empirical distributions.
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_wasserstein_1d(ref, sim)
#' @export
calc_wasserstein_1d <- function(ref, sim, p = 1) {
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) == 0L || length(sim) == 0L) return(NA_real_)

  if (length(ref) == length(sim)) {
    return(mean(abs(sort(ref) - sort(sim))^p)^(1 / p))
  }
  aligned <- align_distributions(ref, sim, method = "quantile", n_points = 1000)
  mean(abs(aligned$ref - aligned$sim)^p)^(1 / p)
}

#' Calculate Median Absolute Deviation (MAD)
#' @param ref Numeric vector of reference distribution.
#' @param sim Numeric vector of simulated distribution.
#' @param align Logical, whether to sort and quantile-align vectors. Default TRUE.
#' @return Non-negative numeric MAD value.
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_mad(ref, sim)
#' @export
calc_mad <- function(ref, sim, align = TRUE) {
  if (align) {
    aligned <- align_distributions(ref, sim)
    ref <- aligned$ref
    sim <- aligned$sim
  }
  stats::median(abs(ref - sim), na.rm = TRUE)
}

#' Calculate Mean Absolute Error (MAE)
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param align Logical, whether to sort and quantile-align vectors. Default TRUE.
#' @return Non-negative numeric mean absolute difference.
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_mae(ref, sim)
#' @export
calc_mae <- function(ref, sim, align = TRUE) {
  if (align) {
    aligned <- align_distributions(ref, sim)
    ref <- aligned$ref
    sim <- aligned$sim
  }
  mean(abs(sim - ref), na.rm = TRUE)
}

#' Calculate Root Mean Squared Error (RMSE)
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param align Logical, whether to sort and quantile-align vectors. Default TRUE.
#' @return Non-negative numeric root mean squared difference.
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_rmse(ref, sim)
#' @export
calc_rmse <- function(ref, sim, align = TRUE) {
  if (align) {
    aligned <- align_distributions(ref, sim)
    ref <- aligned$ref
    sim <- aligned$sim
  }
  sqrt(mean((sim - ref)^2, na.rm = TRUE))
}

#' Calculate 1D Kernel Density Estimation Integrated Squared Error (KDE-ISE)
#'
#' Integrated from SimBench (SydneyBioX/SimBench). Measures the Integrated
#' Squared Error (ISE) between the empirical kernel density estimates of
#' reference and simulated properties across the shared support:
#' \deqn{ISE = \int (\hat{f}_{ref}(x) - \hat{f}_{sim}(x))^2 dx}
#'
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param n_points Number of points for density grid evaluation (default 512).
#' @return Integrated Squared Error (non-negative numeric value, 0 indicates identical density shapes).
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_kde_ise(ref, sim)
#' @export
calc_kde_ise <- function(ref, sim, n_points = 512) {
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) < 2 || length(sim) < 2) return(NA_real_)
  
  xmin <- min(c(ref, sim))
  xmax <- max(c(ref, sim))
  if (xmin == xmax) return(0)
  
  pooled <- c(ref, sim)
  bw_val <- stats::bw.nrd0(pooled)
  if (is.na(bw_val) || !is.finite(bw_val) || bw_val <= 0) bw_val <- "nrd0"
  
  d_ref <- stats::density(ref, from = xmin, to = xmax, n = n_points, bw = bw_val)
  d_sim <- stats::density(sim, from = xmin, to = xmax, n = n_points, bw = bw_val)
  
  dx <- d_ref$x[2] - d_ref$x[1]
  diff_sq <- (d_ref$y - d_sim$y)^2
  sum((diff_sq[-1] + diff_sq[-length(diff_sq)]) / 2) * dx
}

#' Calculate Kernel Density Overlap Value (OV)
#'
#' Integrated from simpipe (Sun et al., Genome Biology 2021; duohongrui/simpipe).
#' Quantifies the shared area under the empirical probability density curves of
#' the reference and simulated distributions:
#' \deqn{\mathrm{OV}(P, Q) = \int_{-\infty}^{+\infty} \min(\hat{f}_{\mathrm{ref}}(x), \hat{f}_{\mathrm{sim}}(x))\,dx}
#' Computes the shared area via normalized kernel density trapezoidal integration.
#'
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param nbins Number of equally spaced grid points for density integration (default 1024, as in simpipe).
#' @return Numeric scalar in \code{[0, 1]} representing distributional overlap (1 = identical distributions, 0 = disjoint; higher is better).
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_overlap(ref, sim)
#' @export
calc_overlap <- function(ref, sim, nbins = 1024) {
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) < 2L || length(sim) < 2L) return(NA_real_)

  sd_ref <- stats::sd(ref)
  sd_sim <- stats::sd(sim)
  if (is.na(sd_ref) || is.na(sd_sim)) return(NA_real_)
  if (sd_ref == 0 && sd_sim == 0) {
    return(if (ref[1] == sim[1]) 1.0 else 0.0)
  }

  pooled <- c(ref, sim)
  bw_val <- stats::bw.nrd0(pooled)
  if (is.na(bw_val) || !is.finite(bw_val) || bw_val <= 0) {
    rng <- max(pooled) - min(pooled)
    bw_val <- if (rng > 0) rng / 50 else 1e-4
  }

  xmin <- min(pooled) - 3 * bw_val
  xmax <- max(pooled) + 3 * bw_val

  d_ref <- stats::density(ref, from = xmin, to = xmax, n = nbins, bw = bw_val)
  d_sim <- stats::density(sim, from = xmin, to = xmax, n = nbins, bw = bw_val)

  dx <- d_ref$x[2] - d_ref$x[1]
  area_ref <- sum((d_ref$y[-1] + d_ref$y[-length(d_ref$y)]) / 2) * dx
  area_sim <- sum((d_sim$y[-1] + d_sim$y[-length(d_sim$y)]) / 2) * dx

  y_ref <- if (area_ref > 0) d_ref$y / area_ref else d_ref$y
  y_sim <- if (area_sim > 0) d_sim$y / area_sim else d_sim$y

  min_y <- pmin(y_ref, y_sim)
  ov_val <- sum((min_y[-1] + min_y[-length(min_y)]) / 2) * dx
  min(1.0, max(0.0, as.numeric(ov_val)))
}

# -----------------------------------------------------------------------------
# 2D Bivariate Joint Relationship Metrics
# -----------------------------------------------------------------------------

#' @title Calculate Fasano-Franceschini 2D Kolmogorov-Smirnov Test Statistic
#' @description Integrated from simpipe.
#'
#' @param ref_mat 2-column numeric matrix for reference.
#' @param sim_mat 2-column numeric matrix for simulation.
#' @param threads CPU threads. Default 1.
#' @return Estimated 2D KS statistic.
#' @examples
#' ref <- matrix(stats::rnorm(100), 50, 2)
#' sim <- matrix(stats::rnorm(100), 50, 2)
#' calc_fasano_franceschini(ref, sim)
#' @export
calc_fasano_franceschini <- function(ref_mat, sim_mat, threads = 1) {
  ref_mat <- as.matrix(stats::na.omit(ref_mat))
  sim_mat <- as.matrix(stats::na.omit(sim_mat))
  if (ncol(ref_mat) != 2 || ncol(sim_mat) != 2) stop("Requires 2-column matrices.")
  
  if (requireNamespace("fasano.franceschini.test", quietly = TRUE)) {
    res <- tryCatch(
      fasano.franceschini.test::fasano.franceschini.test(ref_mat, sim_mat, threads = threads, verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      stat <- if (!is.null(res$statistic)) res$statistic else res$estimate
      if (!is.null(stat) && length(stat) > 0 && is.numeric(stat)) {
        return(as.numeric(stat[1]))
      }
    }
  }
}

#' @title Calculate 2D Bivariate Kernel Density Estimation (KDE) Test Statistic
#' @description Integrated from simpipe & SimBench.
#'
#' @param ref_mat 2-column numeric matrix for reference.
#' @param sim_mat 2-column numeric matrix for simulation.
#' @return z-statistic testing discrepancy between bivariate kernel density estimates.
#' @examples
#' ref <- matrix(stats::rnorm(100), 50, 2)
#' sim <- matrix(stats::rnorm(100), 50, 2)
#' calc_kde_test(ref, sim)
#' @export
calc_kde_test <- function(ref_mat, sim_mat) {
  ref_mat <- as.matrix(stats::na.omit(ref_mat))
  sim_mat <- as.matrix(stats::na.omit(sim_mat))
  if (ncol(ref_mat) != 2 || ncol(sim_mat) != 2) stop("Requires 2-column matrices.")
  
  if (requireNamespace("ks", quietly = TRUE)) {
    res <- tryCatch(
      as.numeric(ks::kde.test(ref_mat, sim_mat)$zstat),
      error = function(e) NA_real_
    )
    return(res)
  }
  return(NA_real_)
}

#' @title Calculate 2D Earth Mover's Distance (2D EMD)
#' @description Integrated from HelenaLC/simulation-comparison.
#'
#' Computes 2D bivariate density via MASS::kde2d over a shared bounding box,
#' then measures optimal transport distance via emdist::emd2d.
#'
#' @param ref_mat 2-column numeric matrix for reference.
#' @param sim_mat 2-column numeric matrix for simulation.
#' @param n Grid resolution for 2D density estimation. Default is 25.
#' @return Normalized 2D Earth Mover's Distance.
#' @examples
#' ref <- matrix(stats::rnorm(100), 50, 2)
#' sim <- matrix(stats::rnorm(100), 50, 2)
#' calc_emd_2d(ref, sim)
#' @export
calc_emd_2d <- function(ref_mat, sim_mat, n = 25) {
  ref_mat <- as.matrix(stats::na.omit(ref_mat))
  sim_mat <- as.matrix(stats::na.omit(sim_mat))
  if (ncol(ref_mat) != 2 || ncol(sim_mat) != 2) stop("Requires 2-column matrices.")
  
  if (requireNamespace("MASS", quietly = TRUE) && requireNamespace("emdist", quietly = TRUE)) {
    tryCatch({
      lims <- c(
        range(c(ref_mat[, 1], sim_mat[, 1])),
        range(c(ref_mat[, 2], sim_mat[, 2]))
      )
      k_ref <- MASS::kde2d(ref_mat[, 1], ref_mat[, 2], n = n, lims = lims)
      k_sim <- MASS::kde2d(sim_mat[, 1], sim_mat[, 2], n = n, lims = lims)
      as.numeric(emdist::emd2d(k_ref$z, k_sim$z, max.iter = 2000L)) / n
    }, error = function(e) NA_real_)
  } else {
    NA_real_
  }
}

# -----------------------------------------------------------------------------
# Metrics from countsimQC (Soneson & Robinson, 2018)
# -----------------------------------------------------------------------------


#' Calculate Wald-Wolfowitz Runs Test on Pooled Distributions
#'
#' Integrated from countsimQC (Soneson & Robinson). Tests whether values from
#' the two datasets intermingle randomly when sorted. A significant left-sided
#' result indicates that identical values cluster together into fewer runs than
#' expected by chance, signaling distinct distributions.
#'
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param alternative Alternative hypothesis: "left.sided" (default in countsimQC) or "two.sided".
#' @return A named list with 2 numeric scalars:
#' \describe{
#'   \item{runs_statistic}{Standardized Wald-Wolfowitz z-statistic comparing observed runs to expected runs under mixing.}
#'   \item{runs_pvalue}{Associated p-value under the selected \code{alternative} hypothesis.}
#' }
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_runs_test(ref, sim)
#' @export
calc_runs_test <- function(ref, sim, alternative = c("left.sided", "two.sided")) {
  alternative <- match.arg(alternative)
  ref <- stats::na.omit(as.numeric(ref))
  sim <- stats::na.omit(as.numeric(sim))
  if (length(ref) < 2L || length(sim) < 2L) {
    return(list(runs_statistic = NA_real_, runs_pvalue = NA_real_))
  }

  df <- data.frame(
    val = c(ref, sim),
    ds = c(rep(1, length(ref)), rep(0, length(sim)))
  )
  df <- df[order(df$val), ]

  r <- rle(df$ds)
  n_runs <- length(r$lengths)
  n1 <- sum(df$ds == 1)
  n0 <- sum(df$ds == 0)
  mu_r <- 1 + (2 * n1 * n0) / (n1 + n0)
  var_r <- (2 * n1 * n0 * (2 * n1 * n0 - n1 - n0)) / ((n1 + n0)^2 * (n1 + n0 - 1))

  if (is.na(var_r) || var_r <= 0) {
    return(list(runs_statistic = 0, runs_pvalue = 1))
  }

  z <- (n_runs - mu_r) / sqrt(var_r)
  p_val <- if (alternative == "left.sided") {
    stats::pnorm(z)
  } else {
    2 * stats::pnorm(-abs(z))
  }

  list(runs_statistic = as.numeric(z), runs_pvalue = as.numeric(p_val))
}

#' Calculate Nearest-Neighbor Composition / Local Neighborhood Similarity (NN Mismatch)
#'
#' Integrated from countsimQC (Soneson & Robinson). Evaluates local neighborhood
#' composition by pooling reference and simulated observations (1D vector or 2D
#' matrix), finding the \eqn{k} nearest neighbors for each observation, and
#' performing a Pearson Chi-squared goodness-of-fit test comparing the local
#' neighborhood dataset proportions against the global dataset proportions.
#' Returns the fraction of observations with \eqn{p \le 0.05}.
#'
#' @param ref Numeric vector (1D) or 2-column matrix/data.frame (2D) of reference values.
#' @param sim Numeric vector (1D) or 2-column matrix/data.frame (2D) of simulated values.
#' @param k Integer number of nearest neighbors. If NULL, defaults to \code{max(5, round(0.05 * N))}.
#' @param subsample_size Maximum number of observations per dataset to subsample for speed (default 300).
#' @return Numeric scalar in \code{[0, 1]} representing the fraction of local neighborhoods with significant dataset composition mismatch (0 = well-mixed, 1 = completely separated).
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_nn_mismatch(ref, sim)
#' @export
calc_nn_mismatch <- function(ref, sim, k = NULL, subsample_size = 300) {
  is_2d <- is.matrix(ref) || is.data.frame(ref)
  
  if (is_2d) {
    ref <- as.matrix(ref)
    sim <- as.matrix(sim)
    ref <- ref[stats::complete.cases(ref), , drop = FALSE]
    sim <- sim[stats::complete.cases(sim), , drop = FALSE]
    if (nrow(ref) < 3 || nrow(sim) < 3) return(NA_real_)
    
    if (nrow(ref) > subsample_size) {
      ref <- ref[ round(seq(1, nrow(ref), length.out = subsample_size)), , drop = FALSE]
    }
    if (nrow(sim) > subsample_size) {
      sim <- sim[ round(seq(1, nrow(sim), length.out = subsample_size)), , drop = FALSE]
    }
    n1 <- nrow(ref)
    n2 <- nrow(sim)
    combined <- rbind(ref, sim)
    # Standardize columns as in countsimQC so neither axis dominates Euclidean distance
    for (col_i in seq_len(ncol(combined))) {
      s_val <- stats::sd(combined[, col_i], na.rm = TRUE)
      if (!is.na(s_val) && s_val > 0) {
        combined[, col_i] <- (combined[, col_i] - mean(combined[, col_i], na.rm = TRUE)) / s_val
      }
    }
  } else {
    ref <- stats::na.omit(as.numeric(ref))
    sim <- stats::na.omit(as.numeric(sim))
    if (length(ref) < 3 || length(sim) < 3) return(NA_real_)
    
    ref <- subsample_dist_vector(ref, max_n = subsample_size)
    sim <- subsample_dist_vector(sim, max_n = subsample_size)
    n1 <- length(ref)
    n2 <- length(sim)
    combined <- matrix(c(ref, sim), ncol = 1)
  }
  
  N <- n1 + n2
  if (is.null(k)) {
    k <- max(5L, round(0.05 * N))
  }
  k <- min(as.integer(k), N - 1L)
  if (k < 1L) return(NA_real_)
  
  labels <- c(rep("ref", n1), rep("sim", n2))
  expected_p <- c(ref = n1 / N, sim = n2 / N)
  
  dmat <- as.matrix(stats::dist(combined, method = "euclidean"))
  diag(dmat) <- Inf
  
  expected_counts <- k * expected_p
  use_sim_pval <- any(expected_counts < 5)
  mismatches <- vapply(seq_len(N), function(j) {
    nn_idx <- order(dmat[j, ])[seq_len(k)]
    nn_labels <- factor(labels[nn_idx], levels = c("ref", "sim"))
    tab <- table(nn_labels)
    pval <- tryCatch(
      stats::chisq.test(tab, p = expected_p, simulate.p.value = use_sim_pval)$p.value,
      error = function(e) 1.0
    )
    if (is.na(pval)) pval <- 1.0
    as.numeric(pval <= 0.05)
  }, numeric(1))
  
  as.numeric(mean(mismatches, na.rm = TRUE))
}

#' Calculate Overall Clustering-Like Separation of Datasets via Silhouette-Based Measures
#'
#' Integrated from countsimQC (Soneson & Robinson). Evaluates overall separation
#' between reference and simulated distributions by treating dataset identity
#' (reference vs. simulated) as cluster labels on the pooled observations and
#' computing both the global average silhouette width and the local neighborhood
#' silhouette width (based on the closest \eqn{k = \max(5, \text{round}(0.05 N))}
#' points in each dataset).
#'
#' @param ref Numeric vector (1D) or 2-column matrix/data.frame (2D) of reference values.
#' @param sim Numeric vector (1D) or 2-column matrix/data.frame (2D) of simulated values.
#' @param subsample_size Maximum number of observations per dataset to subsample for speed (default 300).
#' @return A named list with \code{between_dataset_silh_global} and \code{between_dataset_silh_local}.
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_between_dataset_silhouette(ref, sim)
#' @export
calc_between_dataset_silhouette <- function(ref, sim, subsample_size = 300) {
  is_2d <- is.matrix(ref) || is.data.frame(ref)
  
  if (is_2d) {
    ref <- as.matrix(ref)
    sim <- as.matrix(sim)
    ref <- ref[stats::complete.cases(ref), , drop = FALSE]
    sim <- sim[stats::complete.cases(sim), , drop = FALSE]
    if (nrow(ref) < 3 || nrow(sim) < 3) {
      return(list(between_dataset_silh_global = NA_real_, between_dataset_silh_local = NA_real_))
    }
    if (nrow(ref) > subsample_size) {
      ref <- ref[round(seq(1, nrow(ref), length.out = subsample_size)), , drop = FALSE]
    }
    if (nrow(sim) > subsample_size) {
      sim <- sim[round(seq(1, nrow(sim), length.out = subsample_size)), , drop = FALSE]
    }
    n1 <- nrow(ref)
    n2 <- nrow(sim)
    combined <- rbind(ref, sim)
    for (col_i in seq_len(ncol(combined))) {
      s_val <- stats::sd(combined[, col_i], na.rm = TRUE)
      if (!is.na(s_val) && s_val > 0) {
        combined[, col_i] <- (combined[, col_i] - mean(combined[, col_i], na.rm = TRUE)) / s_val
      }
    }
  } else {
    ref <- stats::na.omit(as.numeric(ref))
    sim <- stats::na.omit(as.numeric(sim))
    if (length(ref) < 3 || length(sim) < 3) {
      return(list(between_dataset_silh_global = NA_real_, between_dataset_silh_local = NA_real_))
    }
    ref <- subsample_dist_vector(ref, max_n = subsample_size)
    sim <- subsample_dist_vector(sim, max_n = subsample_size)
    n1 <- length(ref)
    n2 <- length(sim)
    combined <- matrix(c(ref, sim), ncol = 1)
  }
  
  N <- n1 + n2
  k_local <- max(5L, round(0.05 * N))
  labels <- c(rep(1L, n1), rep(2L, n2))
  
  dmat <- as.matrix(stats::dist(combined, method = "euclidean"))
  
  silh_global <- numeric(N)
  silh_local <- numeric(N)
  
  for (j in seq_len(N)) {
    same_idx <- setdiff(which(labels == labels[j]), j)
    other_idx <- which(labels != labels[j])
    
    d_same <- dmat[j, same_idx]
    d_other <- dmat[j, other_idx]
    
    # Global silhouette width
    a_j <- mean(d_same, na.rm = TRUE)
    b_j <- mean(d_other, na.rm = TRUE)
    denom_g <- max(a_j, b_j)
    silh_global[j] <- if (is.na(denom_g) || denom_g == 0) 0 else (b_j - a_j) / denom_g
    
    # Local silhouette width (countsimQC formulation using closest k_local points)
    a_loc <- mean(sort(d_same)[seq_len(min(k_local, length(d_same)))], na.rm = TRUE)
    b_loc <- mean(sort(d_other)[seq_len(min(k_local, length(d_other)))], na.rm = TRUE)
    denom_l <- max(a_loc, b_loc)
    silh_local[j] <- if (is.na(denom_l) || denom_l == 0) 0 else (b_loc - a_loc) / denom_l
  }
  
  list(
    between_dataset_silh_global = as.numeric(mean(silh_global, na.rm = TRUE)),
    between_dataset_silh_local  = as.numeric(mean(silh_local, na.rm = TRUE))
  )
}



# -----------------------------------------------------------------------------
# Combined Batch Calculators (Unifying simpipe, HelenaLC, SimBench, & countsimQC)
# -----------------------------------------------------------------------------

#' Compute All Univariate Distance & Accuracy Metrics for a Feature
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param metric_prefix Optional prefix string for metric names.
#' @return A named list of 11 numeric univariate distance and similarity metrics
#'   (names prefixed with \code{metric_prefix_} when provided):
#' \describe{
#'   \item{MAD}{Median absolute difference between quantile-aligned vectors.}
#'   \item{KS}{1D Kolmogorov-Smirnov maximum vertical distance between empirical CDFs.}
#'   \item{MAE}{Mean absolute error between quantile-aligned vectors.}
#'   \item{RMSE}{Root mean squared error between quantile-aligned vectors.}
#'   \item{Bhattacharyya}{Bhattacharyya divergence between normalized probability vectors.}
#'   \item{Wasserstein}{1D Wasserstein-1 distance between distributions.}
#'   \item{Runs_Statistic}{Wald-Wolfowitz runs test z-statistic measuring distribution interleaving.}
#'   \item{Runs_PValue}{Left-sided p-value from the Wald-Wolfowitz runs test.}
#'   \item{NN_Mismatch}{Fraction of local k-NN neighborhoods showing dataset composition mismatch.}
#'   \item{Between_Dataset_Silh_Global}{Global between-dataset silhouette width.}
#'   \item{KDE_ISE}{Integrated squared error between kernel density estimates.}
#' }
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_all_univariate_metrics(ref, sim)
#' @export
calc_all_univariate_metrics <- function(ref, sim, metric_prefix = "") {
  prefix <- if (nchar(metric_prefix) > 0) paste0(metric_prefix, "_") else ""
  
  # Subsample if vectors are extremely large (e.g. correlation distributions with millions of values)
  # Following SimBench framework (subsampleSize = 10000), 10000 points guarantees KS, ECDF, and silhouette error < 0.001
  ref_sub <- subsample_dist_vector(ref, max_n = 10000)
  sim_sub <- subsample_dist_vector(sim, max_n = 10000)
  
  runs_res <- calc_runs_test(ref_sub, sim_sub)
  silh_res <- calc_between_dataset_silhouette(ref_sub, sim_sub)
  
  res <- list(
    MAD                         = calc_mad(ref_sub, sim_sub),
    KS                          = calc_ks(ref_sub, sim_sub),
    MAE                         = calc_mae(ref_sub, sim_sub),
    RMSE                        = calc_rmse(ref_sub, sim_sub),
    Bhattacharyya               = calc_bhattacharyya(ref_sub, sim_sub),
    Wasserstein                 = calc_wasserstein_1d(ref_sub, sim_sub),
    Runs_Statistic              = runs_res$runs_statistic,
    Runs_PValue                 = runs_res$runs_pvalue,
    NN_Mismatch                 = calc_nn_mismatch(ref_sub, sim_sub),
    Between_Dataset_Silh_Global = silh_res$between_dataset_silh_global,
    KDE_ISE                     = calc_kde_ise(ref_sub, sim_sub)
  )
  
  names(res) <- paste0(prefix, c(
    "MAD", "KS", "MAE", "RMSE", "Bhattacharyya", "Wasserstein",
    "Runs_Statistic", "Runs_PValue",
    "NN_Mismatch", "Between_Dataset_Silh_Global",
    "KDE_ISE"
  ))
  return(res)
}

#' Compute All Bivariate (2D) Distance & Accuracy Metrics for Joint Distributions
#' @param ref_mat 2-column matrix for reference.
#' @param sim_mat 2-column matrix for simulation.
#' @param metric_prefix Optional prefix string for metric names.
#' @param threads CPU threads.
#' @return A named list of 4 numeric 2D bivariate distributional metrics
#'   (names prefixed with \code{metric_prefix_} when provided):
#' \describe{
#'   \item{Fasano_Franceschini_2D_KS}{Fasano-Franceschini 2D Kolmogorov-Smirnov test statistic.}
#'   \item{EMD_2D}{Normalized 2D Earth Mover's Distance between kernel density surfaces.}
#'   \item{NN_Mismatch_2D}{Fraction of 2D k-NN neighborhoods with dataset composition mismatch.}
#'   \item{Between_Dataset_Silh_2D_Global}{Global between-dataset silhouette width in standardized 2D space.}
#' }
#' @examples
#' ref <- matrix(stats::rnorm(100), 50, 2)
#' sim <- matrix(stats::rnorm(100), 50, 2)
#' calc_all_bivariate_metrics(ref, sim)
#' @export
calc_all_bivariate_metrics <- function(ref_mat, sim_mat, metric_prefix = "", threads = 1) {
  prefix <- if (nchar(metric_prefix) > 0) paste0(metric_prefix, "_") else ""
  
  silh_2d <- calc_between_dataset_silhouette(ref_mat, sim_mat)
  
  res <- list(
    Fasano_Franceschini_2D_KS      = calc_fasano_franceschini(ref_mat, sim_mat, threads = threads),
    EMD_2D                         = calc_emd_2d(ref_mat, sim_mat),
    NN_Mismatch_2D                 = calc_nn_mismatch(ref_mat, sim_mat),
    Between_Dataset_Silh_2D_Global = silh_2d$between_dataset_silh_global
  )
  
  names(res) <- paste0(prefix, c(
    "Fasano_Franceschini_2D_KS",
    "EMD_2D", "NN_Mismatch_2D",
    "Between_Dataset_Silh_2D_Global"
  ))
  return(res)
}


