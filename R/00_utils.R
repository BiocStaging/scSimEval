#' Align Two Empirical Distributions
#'
#' Sorts both distributions and aligns their lengths using either quantile interpolation
#' or nearest-neighbor resampling, ensuring element-wise distance metrics (MAD, MAE, RMSE)
#' are mathematically well-defined even when the number of cells or features differs.
#'
#' @param ref Numeric vector of reference values.
#' @param sim Numeric vector of simulated values.
#' @param method Alignment method: "quantile" (default, interpolated quantiles) or "resample" (resampling).
#' @param n_points Number of points to evaluate when method = "quantile". Default is max(length(ref), length(sim)).
#'
#' @return A list with aligned numeric vectors: \code{list(ref = ..., sim = ...)}.
#' @keywords internal
#' @noRd
align_distributions <- function(ref, sim, method = c("quantile", "resample"), n_points = NULL) {
  method <- match.arg(method)
  ref <- sort(stats::na.omit(as.numeric(ref)))
  sim <- sort(stats::na.omit(as.numeric(sim)))
  
  n_ref <- length(ref)
  n_sim <- length(sim)
  
  if (n_ref == 0 || n_sim == 0) {
    return(list(ref = numeric(0), sim = numeric(0)))
  }
  
  if (n_ref == n_sim) {
    return(list(ref = ref, sim = sim))
  }
  
  if (method == "quantile") {
    if (is.null(n_points)) {
      n_points <- min(max(n_ref, n_sim), 2500)
    }
    probs <- seq(0, 1, length.out = n_points)
    aligned_ref <- stats::quantile(ref, probs = probs, names = FALSE, type = 7)
    aligned_sim <- stats::quantile(sim, probs = probs, names = FALSE, type = 7)
    return(list(ref = aligned_ref, sim = aligned_sim))
  } else {
    diff_len <- abs(n_ref - n_sim)
    if (n_ref > n_sim) {
      add_vals <- sample(sim, diff_len, replace = TRUE)
      sim <- sort(c(sim, add_vals))
    } else {
      add_vals <- sample(ref, diff_len, replace = TRUE)
      ref <- sort(c(ref, add_vals))
    }
    return(list(ref = ref, sim = sim))
  }
}

#' Calculate Proportion of Outliers in a Numeric Vector
#'
#' Uses the standard Tukey 1.5 * IQR criterion on quantiles, as used in simpipe.
#'
#' @param x Numeric vector.
#' @return Outlier proportion (between 0 and 1).
#' @examples
#' ref <- stats::rnorm(50)
#' sim <- stats::rnorm(50)
#' calc_outlier_proportion(ref)
#' @export
calc_outlier_proportion <- function(x) {
  x <- stats::na.omit(as.numeric(x))
  if (length(x) < 4) return(0)
  q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75))
  iqr <- stats::IQR(x)
  lower_bound <- q[1] - 1.5 * iqr
  upper_bound <- q[3] + 1.5 * iqr
  sum(x < lower_bound | x > upper_bound) / length(x)
}

#' Hungarian Maximum Weight Matching for Cluster Assignment
#'
#' Matches predicted cluster labels to ground truth labels using the Hungarian
#' algorithm (Kuhn-Munkres) on the contingency/F1 matrix.
#' Adapted from HelenaLC/simulation-comparison and simpipe.
#'
#' @param pred Vector of predicted cluster labels.
#' @param truth Vector of ground truth cell-type labels.
#'
#' @return A list with matched cluster pairs, precision, recall, and macro F1 score.
#' @keywords internal
#' @noRd
hungarian_match <- function(pred, truth) {
  valid <- !is.na(pred) & !is.na(truth)
  pred <- as.character(pred[valid])
  truth <- as.character(truth[valid])
  
  u_pred <- unique(pred)
  u_truth <- unique(truth)
  
  # Compute pairwise F1 matrix
  f1_mat <- matrix(0, nrow = length(u_pred), ncol = length(u_truth),
                   dimnames = list(u_pred, u_truth))
  pr_mat <- re_mat <- f1_mat
  
  for (i in seq_along(u_pred)) {
    p_mask <- pred == u_pred[i]
    n_p <- sum(p_mask)
    for (j in seq_along(u_truth)) {
      t_mask <- truth == u_truth[j]
      n_t <- sum(t_mask)
      n_overlap <- sum(p_mask & t_mask)
      
      pr <- if (n_p > 0) n_overlap / n_p else 0
      re <- if (n_t > 0) n_overlap / n_t else 0
      f1 <- if (pr + re > 0) 2 * (pr * re) / (pr + re) else 0
      
      pr_mat[i, j] <- pr
      re_mat[i, j] <- re
      f1_mat[i, j] <- f1
    }
  }
  
  # Optimal matching
  matched_pairs <- data.frame()
  if (requireNamespace("clue", quietly = TRUE)) {
    # Maximize total F1 score using solve_LSAP
    if (nrow(f1_mat) <= ncol(f1_mat)) {
      assignment <- clue::solve_LSAP(f1_mat, maximum = TRUE)
      for (i in seq_along(assignment)) {
        j <- as.numeric(assignment[i])
        matched_pairs <- rbind(matched_pairs, data.frame(
          pred = u_pred[i], truth = u_truth[j],
          precision = pr_mat[i, j], recall = re_mat[i, j], F1 = f1_mat[i, j],
          stringsAsFactors = FALSE
        ))
      }
    } else {
      assignment <- clue::solve_LSAP(t(f1_mat), maximum = TRUE)
      for (j in seq_along(assignment)) {
        i <- as.numeric(assignment[j])
        matched_pairs <- rbind(matched_pairs, data.frame(
          pred = u_pred[i], truth = u_truth[j],
          precision = pr_mat[i, j], recall = re_mat[i, j], F1 = f1_mat[i, j],
          stringsAsFactors = FALSE
        ))
      }
    }
  } else {
    # Greedy fallback if clue is not installed
    used_t <- c()
    for (i in seq_along(u_pred)) {
      avail_j <- setdiff(seq_along(u_truth), used_t)
      if (length(avail_j) == 0) break
      best_j <- avail_j[which.max(f1_mat[i, avail_j])]
      used_t <- c(used_t, best_j)
      matched_pairs <- rbind(matched_pairs, data.frame(
        pred = u_pred[i], truth = u_truth[best_j],
        precision = pr_mat[i, best_j], recall = re_mat[i, best_j], F1 = f1_mat[i, best_j],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Align predicted cluster labels to ground truth using matched pairs
  mapping <- stats::setNames(matched_pairs$truth, matched_pairs$pred)
  aligned_pred <- unname(mapping[as.character(pred)])
  acc <- mean(aligned_pred == as.character(truth), na.rm = TRUE)
  
  list(
    matched_pairs = matched_pairs,
    accuracy = acc,
    mean_precision = mean(matched_pairs$precision, na.rm = TRUE),
    mean_recall = mean(matched_pairs$recall, na.rm = TRUE),
    mean_F1 = mean(matched_pairs$F1, na.rm = TRUE),
    aligned_pred = aligned_pred
  )
}

# -----------------------------------------------------------------------------
# High-Performance Computational Acceleration Helpers
# (Vectorized Sparse Matrix Arithmetic, BLAS Distance, & Truncated SVD)
# -----------------------------------------------------------------------------

#' Fast Row-Wise Variances for Dense or Sparse Single-Cell Matrices
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of row variances.
#' @export
fast_row_vars <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    n_cols <- ncol(mat)
    if (n_cols <= 1) return(rep(0, nrow(mat)))
    rm <- Matrix::rowMeans(mat)
    mat_sq <- mat
    mat_sq@x <- mat_sq@x^2
    sq_rm <- Matrix::rowMeans(mat_sq)
    v <- (sq_rm - rm^2) * (n_cols / (n_cols - 1))
    v[v < 0] <- 0
    return(as.numeric(v))
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(as.numeric(matrixStats::rowVars(as.matrix(mat))))
  } else {
    n_cols <- ncol(mat)
    if (n_cols <= 1) return(rep(0, nrow(mat)))
    rm <- rowMeans(mat)
    sq_rm <- rowMeans(mat^2)
    v <- (sq_rm - rm^2) * (n_cols / (n_cols - 1))
    v[v < 0] <- 0
    return(as.numeric(v))
  }
}

#' Fast Column-Wise Variances for Dense or Sparse Single-Cell Matrices
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of column variances.
#' @export
fast_col_vars <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    n_rows <- nrow(mat)
    if (n_rows <= 1) return(rep(0, ncol(mat)))
    cm <- Matrix::colMeans(mat)
    mat_sq <- mat
    mat_sq@x <- mat_sq@x^2
    sq_cm <- Matrix::colMeans(mat_sq)
    v <- (sq_cm - cm^2) * (n_rows / (n_rows - 1))
    v[v < 0] <- 0
    return(as.numeric(v))
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(as.numeric(matrixStats::colVars(as.matrix(mat))))
  } else {
    n_rows <- nrow(mat)
    if (n_rows <= 1) return(rep(0, ncol(mat)))
    cm <- colMeans(mat)
    sq_cm <- colMeans(mat^2)
    v <- (sq_cm - cm^2) * (n_rows / (n_rows - 1))
    v[v < 0] <- 0
    return(as.numeric(v))
  }
}

#' Fast Zero Fraction per Cell (Column-wise)
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of zero fraction per cell.
#' @export
fast_sparse_zero_fraction_cell <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    nz_per_col <- diff(mat@p)
    return(as.numeric((nrow(mat) - nz_per_col) / nrow(mat)))
  }
  as.numeric(colMeans(mat == 0))
}

#' Fast Zero Fraction per Feature (Row-wise)
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of zero fraction per feature.
#' @export
fast_sparse_zero_fraction_feature <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    nz_idx <- which(mat@x != 0)
    nz_per_row <- tabulate(mat@i[nz_idx] + 1, nbins = nrow(mat))
    return(as.numeric((ncol(mat) - nz_per_row) / ncol(mat)))
  }
  as.numeric(rowMeans(mat == 0))
}

#' Fast BLAS-Accelerated Cross-Distance Squared Matrix
#'
#' Computes ||X_i - Y_j||^2 using BLAS matrix product: ||x||^2 + ||y||^2 - 2 <x, y>.
#' @param X Matrix (m x p).
#' @param Y Matrix (n x p).
#' @return Matrix of squared Euclidean distances (m x n).
#' @export
fast_cross_distance_sq <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  x_norm <- rowSums(X^2)
  y_norm <- rowSums(Y^2)
  d2 <- outer(x_norm, y_norm, "+") - 2 * (X %*% t(Y))
  d2[d2 < 0] <- 0
  d2
}

#' Fast Truncated PCA on Top Highly Variable Features
#' @param mat Matrix or dgCMatrix (cells x features).
#' @param n_pcs Number of principal components. Default 20.
#' @param n_hvg Number of top highly variable genes to retain. Default 2000.
#' @param scale Logical, whether to scale features. Default FALSE.
#' @param normalize Logical, whether to apply library-size log-normalization if raw counts. Default FALSE.
#' @return Matrix of cell principal component embeddings (cells x n_pcs).
#' @export
fast_pca <- function(mat, n_pcs = 20, n_hvg = 2000, scale = FALSE, normalize = FALSE) {
  if (inherits(mat, "Matrix")) {
    mat <- as.matrix(mat)
  }
  if (normalize) {
    # mat is cells x features; normalize_counts expects features x cells
    mat <- t(as.matrix(normalize_counts(t(mat))))
  }
  n_feats <- ncol(mat)
  if (n_feats > n_hvg) {
    feat_vars <- fast_col_vars(mat)
    top_feats <- order(feat_vars, decreasing = TRUE)[seq_len(n_hvg)]
    mat <- mat[, top_feats, drop = FALSE]
  }
  rank_use <- min(n_pcs, nrow(mat) - 1, ncol(mat) - 1)
  if (rank_use < 2) return(as.matrix(mat)[, seq_len(min(1, ncol(mat))), drop = FALSE])
  
  min_dim <- min(nrow(mat), ncol(mat))
  if (requireNamespace("irlba", quietly = TRUE) && rank_use >= 3 && min_dim > 50 && rank_use < 0.5 * min_dim) {
    res <- tryCatch(suppressWarnings(irlba::prcomp_irlba(as.matrix(mat), n = rank_use, scale. = scale)), error = function(e) NULL)
    if (!is.null(res)) return(res$x)
  }
  stats::prcomp(as.matrix(mat), rank. = rank_use, scale. = scale)$x
}

#' Subsample Large Vectors for Accelerated Distributional Distance Estimation
#' @param x Numeric vector.
#' @param max_n Maximum sample size (default 3000).
#' @return Subsampled numeric vector if length(x) > max_n, else x.
#' @export
subsample_dist_vector <- function(x, max_n = 3000) {
  x <- stats::na.omit(as.numeric(x))
  if (length(x) > max_n) {
    sample(x, max_n)
  } else {
    x
  }
}

#' Normalize Single-Cell Count Matrix (Library Size + Log Transformation)
#'
#' Converts raw single-cell sequencing count matrices to library-size normalized,
#' log-transformed values (e.g. log2(CPM + 1) or log(normalized counts + 1)),
#' mitigating library size confounding across distance, correlation, and manifold metrics
#' as established across benchmark literature (SimCH, scDesign3, SPsimSeq, GDSim, SymBench, cfDiffusion, scCross).
#'
#' @param mat Expression or count matrix (features x cells), dgCMatrix, SingleCellExperiment, or Seurat object.
#' @param scale_factor Numeric scaling factor for counts per cell. Default is 1e6 (CPM); can use 1e4.
#' @param log_base Numeric base for the logarithm: 2 (default, for log2(CPM+1)) or natural log \code{exp(1)}.
#' @param pseudo_count Numeric pseudo-count added before logarithm. Default is 1.
#' @return Normalized matrix of log-transformed expression values.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' norm_mat <- normalize_counts(example_scrna$ref)
#' @export
normalize_counts <- function(mat, scale_factor = 1e6, log_base = 2, pseudo_count = 1) {
  if (is.null(mat)) return(NULL)
  if (inherits(mat, "SingleCellExperiment")) {
    mat <- SummarizedExperiment::assay(mat, "counts")
  } else if (inherits(mat, "Seurat")) {
    mat <- Seurat::GetAssayData(mat, slot = "counts")
  }
  
  # Check if matrix is already log-normalized (max < 35 and has fractional values)
  sample_vals <- if (inherits(mat, "dgCMatrix")) head(mat@x, 100) else as.numeric(mat[seq_len(min(10, nrow(mat))), seq_len(min(10, ncol(mat)))])
  sample_vals <- sample_vals[!is.na(sample_vals)]
  is_already_lognorm <- length(sample_vals) > 0 && max(sample_vals) < 35 && any(sample_vals %% 1 != 0)
  if (is_already_lognorm) {
    return(mat)
  }
  
  libs <- as.numeric(Matrix::colSums(mat))
  libs[libs == 0] <- 1
  
  if (inherits(mat, "dgCMatrix")) {
    norm_mat <- mat
    cell_factors <- scale_factor / libs
    col_idx <- rep(seq_len(ncol(mat)), diff(mat@p))
    scaled_x <- norm_mat@x * cell_factors[col_idx] + pseudo_count
    norm_mat@x <- if (log_base == 2) log2(scaled_x) else log(scaled_x)
    return(norm_mat)
  } else {
    mat <- as.matrix(mat)
    scaled <- t(t(mat) / libs) * scale_factor + pseudo_count
    norm_mat <- if (log_base == 2) log2(scaled) else log(scaled)
    return(norm_mat)
  }
}



