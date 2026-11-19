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
  n_cols <- ncol(mat)
  if (n_cols <= 1L) return(rep(0, nrow(mat)))
  if (inherits(mat, "Matrix")) {
    rm <- Matrix::rowMeans(mat)
    sq_rm <- Matrix::rowMeans(mat^2)
    v <- (sq_rm - rm^2) * (n_cols / (n_cols - 1))
    v[v < 0] <- 0
    return(as.numeric(v))
  }
  as.numeric(matrixStats::rowVars(as.matrix(mat)))
}

#' Fast Column-Wise Variances for Dense or Sparse Single-Cell Matrices
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of column variances.
#' @export
fast_col_vars <- function(mat) {
  n_rows <- nrow(mat)
  if (n_rows <= 1L) return(rep(0, ncol(mat)))
  if (inherits(mat, "Matrix")) {
    cm <- Matrix::colMeans(mat)
    sq_cm <- Matrix::colMeans(mat^2)
    v <- (sq_cm - cm^2) * (n_rows / (n_rows - 1))
    v[v < 0] <- 0
    return(as.numeric(v))
  }
  as.numeric(matrixStats::colVars(as.matrix(mat)))
}

#' Fast Zero Fraction per Cell (Column-wise)
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of zero fraction per cell.
#' @export
fast_sparse_zero_fraction_cell <- function(mat) {
  if (inherits(mat, "Matrix")) {
    nz_per_col <- as.numeric(Matrix::colSums(mat != 0))
    return((nrow(mat) - nz_per_col) / nrow(mat))
  }
  as.numeric(colMeans(mat == 0))
}

#' Fast Zero Fraction per Feature (Row-wise)
#' @param mat Matrix or dgCMatrix (features x cells).
#' @return Numeric vector of zero fraction per feature.
#' @export
fast_sparse_zero_fraction_feature <- function(mat) {
  if (inherits(mat, "Matrix")) {
    nz_per_row <- as.numeric(Matrix::rowSums(mat != 0))
    return((ncol(mat) - nz_per_row) / ncol(mat))
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

.run_prcomp_irlba <- function(X, n, center = TRUE, scale. = FALSE) {
  if (!requireNamespace("irlba", quietly = TRUE)) return(NULL)
  if (n >= 0.5 * min(dim(X))) return(NULL)
  v_init <- rep(1 / sqrt(max(1L, ncol(X))), ncol(X))
  tryCatch(
    irlba::prcomp_irlba(X, n = n, center = center, scale. = scale., v = v_init),
    error = function(e) NULL
  )
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
  if (rank_use >= 3 && min_dim > 50 && rank_use < 0.5 * min_dim) {
    res <- .run_prcomp_irlba(as.matrix(mat), n = rank_use, center = TRUE, scale. = scale)
    if (!is.null(res)) return(res$x)
  }
  stats::prcomp(as.matrix(mat), rank. = rank_use, scale. = scale)$x
}

#' Subsample Large Vectors for Accelerated Distributional Distance Estimation
#' @param x Numeric vector.
#' @param max_n Maximum sample size (default 10000, matching SimBench framework).
#' @return Subsampled numeric vector if length(x) > max_n, else x.
#' @export
subsample_dist_vector <- function(x, max_n = 10000) {
  x <- stats::na.omit(as.numeric(x))
  if (length(x) > max_n) {
    sample(x, max_n)
  } else {
    x
  }
}

#' Detect Whether Single-Cell Input Represents Chromatin Accessibility (scATAC-seq)
#'
#' Checks explicit modality strings, \code{SingleCellExperiment} experiment names,
#' or peak coordinate feature names (e.g., \code{Peak_01}, \code{chr1:100-200}).
#'
#' @param data Matrix, \code{dgCMatrix}, or \code{SingleCellExperiment} object.
#' @param modality Optional character modality hint.
#' @return Logical scalar (\code{TRUE} for scATAC-seq, \code{FALSE} otherwise).
#' @keywords internal
#' @noRd
.is_atac_modality <- function(data = NULL, modality = NULL) {
  if (!is.null(modality) && length(modality) > 0 && !is.na(modality[1]) && nzchar(as.character(modality[1]))) {
    m_low <- tolower(as.character(modality[1]))
    if (grepl("atac|peak|chromatin", m_low)) return(TRUE)
    if (grepl("rna|gene|transcript", m_low)) return(FALSE)
  }

  if (is.null(data)) return(FALSE)

  if (inherits(data, "SingleCellExperiment")) {
    exp_nm <- tryCatch(SingleCellExperiment::mainExpName(data), error = function(e) NULL)
    if (!is.null(exp_nm) && grepl("atac|peak|chromatin", tolower(exp_nm))) {
      return(TRUE)
    }
  }

  rn <- tryCatch(rownames(data), error = function(e) NULL)
  if (!is.null(rn) && length(rn) > 0) {
    is_peak_name <- grepl("^(peak|chr[0-9xyXYM]+[:-]|[A-Za-z0-9_.-]+[:-][0-9]+-[0-9]+)", rn, ignore.case = TRUE)
    if (mean(is_peak_name) > 0.5) return(TRUE)
  }

  FALSE
}

#' Perform TF-IDF Normalization for scATAC-seq Peak Matrices
#'
#' Computes Term Frequency-Inverse Document Frequency (TF-IDF) normalization
#' following the Signac (\code{RunTFIDF(..., method = 1)}) and simPIC
#' (Chugh et al., 2024) standard for scATAC-seq chromatin accessibility matrices:
#' \deqn{\text{TF-IDF}_{i,j} = \log\left(1 + \frac{X_{i,j}}{N_j} \cdot \frac{N_{\text{cells}}}{\sum_k X_{i,k}} \cdot \text{scale\_factor}\right)}
#'
#' @param mat A peak-by-cell count matrix, sparse \code{dgCMatrix}, or \code{SingleCellExperiment} object.
#' @param scale_factor Numeric scaling factor applied prior to \code{log1p}. Default is \code{1e4} (10,000).
#' @return TF-IDF transformed matrix (preserving \code{dgCMatrix} sparsity when sparse input is supplied).
#' @examples
#' set.seed(42)
#' m <- matrix(rpois(100, 0.5), nrow = 10)
#' tfidf_mat <- tfidf_transform(m)
#' @export
tfidf_transform <- function(mat, scale_factor = 1e4) {
  if (is.null(mat)) return(NULL)
  if (inherits(mat, "SingleCellExperiment")) {
    mat <- SummarizedExperiment::assay(mat, "counts")
  } else if (inherits(mat, "Seurat")) {
    stop("Seurat objects are not directly supported; convert via Seurat::as.SingleCellExperiment().", call. = FALSE)
  }

  libs <- as.numeric(Matrix::colSums(mat))
  libs[libs == 0] <- 1

  rs <- as.numeric(Matrix::rowSums(mat))
  idf <- ifelse(rs > 0, ncol(mat) / rs, 0)

  if (inherits(mat, "Matrix")) {
    trip <- Matrix::summary(mat)
    if (nrow(trip) == 0L) return(mat)
    vals <- log1p((trip$x / libs[trip$j]) * idf[trip$i] * scale_factor)
    return(Matrix::sparseMatrix(
      i = trip$i, j = trip$j, x = vals,
      dims = dim(mat), dimnames = dimnames(mat)
    ))
  }
  mat <- as.matrix(mat)
  tf <- t(t(mat) / libs)
  scaled <- sweep(tf, 1L, idf, `*`) * scale_factor
  log1p(scaled)
}

#' Fast Count Normalization (Log2-CPM or TF-IDF)
#'
#' Converts raw single-cell count matrices into log-normalized CPM (\code{"lognorm"})
#' for scRNA-seq or TF-IDF (\code{"tfidf"}) for scATAC-seq, preserving
#' \code{dgCMatrix} sparsity.
#'
#' @param mat A count matrix (genes/features x cells), \code{dgCMatrix}, or \code{SingleCellExperiment} object.
#' @param scale_factor Numeric scaling factor for counts per cell. Default is 1e6 (CPM) for \code{"lognorm"}, or 1e4 for \code{"tfidf"}.
#' @param log_base Numeric base for the logarithm in \code{"lognorm"} mode: 2 (default, for log2(CPM+1)) or natural log \code{exp(1)}.
#' @param pseudo_count Numeric pseudo-count added before logarithm in \code{"lognorm"} mode. Default is 1.
#' @param method Normalization method: \code{"lognorm"} (default) or \code{"tfidf"} (Signac/simPIC scATAC-seq standard).
#' @return Normalized matrix of transformed values.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' norm_mat <- normalize_counts(example_scrna$ref)
#' @export
normalize_counts <- function(mat, scale_factor = 1e6, log_base = 2, pseudo_count = 1,
                             method = c("lognorm", "tfidf")) {
  method <- match.arg(method)
  if (is.null(mat)) return(NULL)
  if (method == "tfidf") {
    sf <- if (missing(scale_factor)) 1e4 else scale_factor
    return(tfidf_transform(mat, scale_factor = sf))
  }
  if (inherits(mat, "SingleCellExperiment")) {
    mat <- SummarizedExperiment::assay(mat, "counts")
  } else if (inherits(mat, "Seurat")) {
    stop("Seurat objects are not directly supported; convert via Seurat::as.SingleCellExperiment().", call. = FALSE)
  }

  if (inherits(mat, "Matrix")) {
    trip <- Matrix::summary(mat)
    sample_vals <- head(trip$x, 100L)
    sample_vals <- sample_vals[!is.na(sample_vals)]
    if (length(sample_vals) > 0L && max(sample_vals) < 35 && any(sample_vals %% 1 != 0)) {
      return(mat)
    }
    if (nrow(trip) == 0L) return(mat)
    libs <- as.numeric(Matrix::colSums(mat))
    libs[libs == 0] <- 1
    cell_factors <- scale_factor / libs
    scaled_x <- trip$x * cell_factors[trip$j] + pseudo_count
    vals <- if (log_base == 2) log2(scaled_x) else log(scaled_x)
    return(Matrix::sparseMatrix(
      i = trip$i, j = trip$j, x = vals,
      dims = dim(mat), dimnames = dimnames(mat)
    ))
  }

  sample_vals <- as.numeric(mat[seq_len(min(10L, nrow(mat))), seq_len(min(10L, ncol(mat)))])
  sample_vals <- sample_vals[!is.na(sample_vals)]
  if (length(sample_vals) > 0L && max(sample_vals) < 35 && any(sample_vals %% 1 != 0)) {
    return(mat)
  }

  libs <- as.numeric(colSums(mat))
  libs[libs == 0] <- 1
  mat <- as.matrix(mat)
  scaled <- t(t(mat) / libs) * scale_factor + pseudo_count
  if (log_base == 2) log2(scaled) else log(scaled)
}

#' Extract Dataset Dimensional Properties and Summary Statistics
#'
#' Computes key dimensional and biological properties from single-cell or multiomics count matrices,
#' including cell counts, feature counts, sparsity percentage, biological group/cell type counts,
#' batch counts, median library size, median detected features per cell, and mean expression level.
#'
#' @param mat A count matrix (genes/features x cells), dgCMatrix, SingleCellExperiment, or Seurat object.
#' @param role Character string describing the role of the dataset (e.g. \code{"Biological Reference"} or \code{"Simulated"}). Default is \code{"Reference"}.
#' @param method_name Character string indicating the method or dataset name (e.g. \code{"Empirical Reference"}, \code{"Splatter"}). Default is \code{"Dataset"}.
#' @param modality Character string describing the molecular modality (e.g. \code{"scRNA-seq"}, \code{"scATAC-seq"}). Default is \code{"scRNA-seq"}.
#' @param cell_types Optional factor or character vector of cell type annotations.
#' @param batch_info Optional factor or character vector of batch annotations.
#'
#' @return A data.frame with 1 row summarizing the dataset properties.
#' @examples
#' data(example_scrna, package = "scSimEval")
#' extract_dataset_summary(example_scrna$ref, role = "Reference", method_name = "Empirical")
#' @export
extract_dataset_summary <- function(mat, role = "Reference", method_name = "Dataset",
                                     modality = "scRNA-seq", cell_types = NULL, batch_info = NULL) {
  if (is.null(mat)) return(NULL)

  if (inherits(mat, "SingleCellExperiment")) {
    cd <- SummarizedExperiment::colData(mat)
    if (is.null(cell_types) && "cell_type" %in% colnames(cd)) {
      cell_types <- cd$cell_type
    }
    if (is.null(batch_info) && "batch" %in% colnames(cd)) {
      batch_info <- cd$batch
    }
    mat <- SingleCellExperiment::counts(mat)
  } else if (inherits(mat, "Seurat")) {
    stop("Seurat objects are not directly supported; convert via Seurat::as.SingleCellExperiment().", call. = FALSE)
  }

  n_cells <- ncol(mat)
  n_features <- nrow(mat)
  total_entries <- as.numeric(n_cells) * as.numeric(n_features)

  if (inherits(mat, "Matrix")) {
    n_nz <- Matrix::nnzero(mat)
    n_zeros <- total_entries - n_nz
    sparsity_pct <- if (total_entries > 0) round((n_zeros / total_entries) * 100, 2) else 0
    lib_sizes <- Matrix::colSums(mat)
    det_feats <- Matrix::colSums(mat > 0)
    mean_expr <- Matrix::mean(mat)
  } else {
    mat_num <- as.matrix(mat)
    sparsity_pct <- round(mean(mat_num == 0) * 100, 2)
    lib_sizes <- colSums(mat_num)
    det_feats <- colSums(mat_num > 0)
    mean_expr <- mean(mat_num)
  }
  
  med_lib <- stats::median(lib_sizes)
  med_det <- stats::median(det_feats)
  
  ct_str <- if (!is.null(cell_types) && length(cell_types) == n_cells) {
    u_ct <- unique(stats::na.omit(as.character(cell_types)))
    formatC(length(u_ct), format = "d", big.mark = ",")
  } else {
    "Not specified"
  }
  
  b_str <- if (!is.null(batch_info) && length(batch_info) == n_cells) {
    u_b <- unique(stats::na.omit(as.character(batch_info)))
    formatC(length(u_b), format = "d", big.mark = ",")
  } else {
    "1"
  }
  
  data.frame(
    "Dataset / Simulator" = as.character(method_name),
    "Role" = as.character(role),
    "Modality" = as.character(modality),
    "Cells (N)" = formatC(as.integer(n_cells), format = "d", big.mark = ","),
    "Features (P)" = formatC(as.integer(n_features), format = "d", big.mark = ","),
    "Sparsity" = sprintf("%.2f%%", sparsity_pct),
    "Cell Types (Groups)" = ct_str,
    "Batches" = b_str,
    "Median Lib Size" = formatC(round(med_lib, 1), format = "f", digits = 1, big.mark = ","),
    "Median Detected Features" = formatC(round(med_det, 0), format = "d", big.mark = ","),
    "Mean Expression" = formatC(round(mean_expr, 3), format = "f", digits = 3),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

.SC_CELLTYPE_COLS <- c(
  "cell_type", "celltype", "cell_types", "celltypes", "cluster", "clusters",
  "group", "groups", "label", "labels", "cell_identity", "cell_line",
  "cellstate", "state", "annotation", "subclass", "cell_ontology_class"
)

.SC_BATCH_COLS <- c(
  "batch", "batch_info", "batches", "replicate", "rep", "donor", "sample",
  "dataset", "protocol", "run", "experiment", "plate", "library", "gem_group"
)

#' Validate and Standardize Single-Cell Input Objects
#'
#' Strictly enforces \code{SingleCellExperiment} compliance in \code{scSimEval}.
#' Validates that input datasets are \code{SingleCellExperiment} objects,
#' extracting underlying count matrices, cell annotations, batch labels, and
#' dimensional reductions. Users starting from \code{Seurat} objects should first
#' convert via \code{Seurat::as.SingleCellExperiment()}.
#'
#' @param obj A \code{SingleCellExperiment} object to validate.
#' @param arg_name Name of the argument being validated (for clear error messages).
#' @param reduction_name Optional name of precomputed dimensional reduction (e.g. \code{"PCA"}, \code{"integrated"}).
#' @param modality Optional character string specifying the modality (e.g. \code{"scRNA-seq"} or \code{"scATAC-seq"}) to guide assay/\code{altExp} selection.
#' @param cell_type_col Column name in \code{colData(obj)} containing cell-type labels (default \code{"cell_type"}).
#' @param batch_col Column name in \code{colData(obj)} containing batch/replicate labels (default \code{"batch"}).
#'
#' @return A named list with 5 elements:
#' \describe{
#'   \item{object}{The validated \code{SingleCellExperiment} object.}
#'   \item{counts}{Extracted raw count matrix (\code{dgCMatrix} or numeric matrix, features x cells).}
#'   \item{cell_types}{Character vector of extracted cell-type annotations, or \code{NULL} if absent.}
#'   \item{batch_info}{Character vector of extracted batch annotations, or \code{NULL} if absent.}
#'   \item{embedding}{Extracted dimensional reduction numeric matrix (cells x dims), or \code{NULL}.}
#' }
#'
#' @export
validate_sc_object <- function(obj,
                               arg_name = "data",
                               reduction_name = NULL,
                               modality = NULL,
                               cell_type_col = "cell_type",
                               batch_col = "batch") {
  if (is.null(obj)) return(NULL)

  if (inherits(obj, "Seurat")) {
    stop(sprintf(
      "scSimEval requires a 'SingleCellExperiment' object for '%s' (received 'Seurat').\nPlease convert your Seurat object using Seurat::as.SingleCellExperiment().",
      arg_name
    ), call. = FALSE)
  }

  if (!inherits(obj, "SingleCellExperiment")) {
    stop(sprintf(
      "scSimEval strictly requires a 'SingleCellExperiment' object for '%s' (received '%s').",
      arg_name, paste(class(obj), collapse = "/")
    ), call. = FALSE)
  }

  counts_mat <- NULL
  ct_vec <- NULL
  bt_vec <- NULL
  embed_mat <- NULL

  # 1. Counts (checking altExp if modality requests secondary assay)
  assays_list <- SummarizedExperiment::assays(obj)
  altexp_names <- tryCatch(SingleCellExperiment::altExpNames(obj), error = function(e) character(0))
  target_exp <- NULL
  if (!is.null(modality) && length(altexp_names) > 0) {
    mod_check <- tolower(trimws(as.character(modality)[1L]))
    if (grepl("atac|peak|chromatin", mod_check)) {
      hit <- grep("atac|peak|chromatin|macs2|tile", altexp_names, ignore.case = TRUE, value = TRUE)
      if (length(hit) > 0) target_exp <- SingleCellExperiment::altExp(obj, hit[1L])
    } else if (grepl("rna|gene|transcr", mod_check)) {
      hit <- grep("^rna$|gene|transcr", altexp_names, ignore.case = TRUE, value = TRUE)
      if (length(hit) == 0) hit <- grep("rna", altexp_names, ignore.case = TRUE, value = TRUE)
      if (length(hit) > 0) target_exp <- SingleCellExperiment::altExp(obj, hit[1L])
    }
  }
  if (!is.null(target_exp)) {
    exp_assays <- SummarizedExperiment::assays(target_exp)
    counts_mat <- if ("counts" %in% names(exp_assays)) exp_assays[["counts"]] else if (length(exp_assays) > 0) exp_assays[[1L]] else NULL
  } else {
    counts_mat <- if ("counts" %in% names(assays_list)) assays_list[["counts"]] else if (length(assays_list) > 0) assays_list[[1L]] else NULL
  }

  # 2. Metadata (cell_type and batch)
  cd <- SummarizedExperiment::colData(obj)
  cd_names <- colnames(cd)
  cd_cols <- tolower(cd_names)
  ct_col <- if (!is.null(cell_type_col) && cell_type_col %in% cd_names) {
    match(cell_type_col, cd_names)
  } else {
    which(cd_cols %in% tolower(.SC_CELLTYPE_COLS))[1]
  }
  if (!is.na(ct_col)) ct_vec <- as.character(cd[[ct_col]])

  bt_col <- if (!is.null(batch_col) && batch_col %in% cd_names) {
    match(batch_col, cd_names)
  } else {
    which(cd_cols %in% tolower(.SC_BATCH_COLS))[1]
  }
  if (!is.na(bt_col)) bt_vec <- as.character(cd[[bt_col]])

  # 3. Dimensional Reductions
  rd_names <- SingleCellExperiment::reducedDimNames(obj)
  if (!is.null(rd_names) && length(rd_names) > 0) {
    target_rd <- if (!is.null(reduction_name) && reduction_name %in% rd_names) {
      reduction_name
    } else {
      mod_pref <- if (!is.null(modality) && grepl("atac|peak|chromatin", tolower(modality))) {
        c("atac_lsi", "lsi", "atac_pca", "atac_umap", "harmony_atac", "integrated", "joint", "pca", "PCA")
      } else if (!is.null(modality) && grepl("rna|gene|transcr", tolower(modality))) {
        c("rna_pca", "pca", "PCA", "rna_umap", "harmony_rna", "integrated", "joint")
      } else {
        c("integrated", "joint", "harmony", "wnn", "PCA", "pca", "lsi")
      }
      found <- intersect(mod_pref, rd_names)
      if (length(found) > 0) found[1] else rd_names[1]
    }
    embed_mat <- as.matrix(SingleCellExperiment::reducedDim(obj, target_rd))
  }

  list(
    object = obj,
    counts = counts_mat,
    cell_types = ct_vec,
    batch_info = bt_vec,
    embedding = embed_mat
  )
}

#' Validate and Standardize Multiomics Input Objects
#'
#' Validates that multiomics inputs conform to \code{SingleCellExperiment} structures.
#' Accepts either a single multi-modal \code{SingleCellExperiment} object with an
#' \code{altExp} for the secondary modality, a named list of two \code{SingleCellExperiment}
#' objects (e.g. \code{list(RNA = rna_sce, ATAC = atac_sce)}), or a file path to an
#' \code{.rds} file containing either structure.
#'
#' @param obj Multiomics \code{SingleCellExperiment} object, list of two \code{SingleCellExperiment} objects, or \code{.rds} path.
#' @param arg_name Name of argument for diagnostic error messages.
#' @param reduction_name Optional reduction name to extract.
#' @param cell_type_col Column name in \code{colData} containing cell-type labels (default \code{"cell_type"}).
#' @param batch_col Column name in \code{colData} containing batch labels (default \code{"batch"}).
#' @return A named list with 9 elements:
#' \describe{
#'   \item{m1}{Count matrix for modality 1 (typically RNA).}
#'   \item{m2}{Count matrix for modality 2 (typically ATAC).}
#'   \item{m1_obj}{\code{SingleCellExperiment} object for modality 1.}
#'   \item{m2_obj}{\code{SingleCellExperiment} object for modality 2.}
#'   \item{m1_name}{Character name of modality 1.}
#'   \item{m2_name}{Character name of modality 2.}
#'   \item{cell_types}{Character vector (paired) or named list (unpaired) of cell-type labels, or \code{NULL}.}
#'   \item{batch_info}{Character vector (paired) or named list (unpaired) of batch labels, or \code{NULL}.}
#'   \item{embedding}{Extracted dimensional reduction matrix (cells x dims), or \code{NULL}.}
#' }
#' @export
validate_multiomics_input <- function(obj,
                                      arg_name = "multi_data",
                                      reduction_name = NULL,
                                      cell_type_col = "cell_type",
                                      batch_col = "batch") {
  if (is.null(obj)) return(NULL)

  if (is.character(obj) && length(obj) == 1L && file.exists(obj)) {
    obj <- readRDS(obj)
  }

  # Case A: Named list of 2 SingleCellExperiment objects
  if (is.list(obj) && !inherits(obj, "Seurat") && !inherits(obj, "SingleCellExperiment")) {
    if (length(obj) == 1L) {
      obj <- obj[[1L]]
    } else if (length(obj) >= 2L) {
      o1 <- obj[[1L]]
      o2 <- obj[[2L]]
      nm1 <- names(obj)[1L] %||% "RNA"
      nm2 <- names(obj)[2L] %||% "ATAC"
      v1 <- validate_sc_object(o1, paste0(arg_name, "[[1]]"), reduction_name = reduction_name,
                               cell_type_col = cell_type_col, batch_col = batch_col)
      v2 <- validate_sc_object(o2, paste0(arg_name, "[[2]]"), reduction_name = reduction_name,
                               cell_type_col = cell_type_col, batch_col = batch_col)

      c1 <- colnames(v1$counts)
      c2 <- colnames(v2$counts)
      is_paired_cells <- !is.null(c1) && !is.null(c2) && length(c1) == length(c2) && identical(c1, c2)

      ct <- if (is_paired_cells) {
        if (!is.null(v1$cell_types)) v1$cell_types else v2$cell_types
      } else if (!is.null(v1$cell_types) || !is.null(v2$cell_types)) {
        stats::setNames(list(v1$cell_types, v2$cell_types), c(nm1, nm2))
      } else {
        NULL
      }

      bt <- if (is_paired_cells) {
        if (!is.null(v1$batch_info)) v1$batch_info else v2$batch_info
      } else if (!is.null(v1$batch_info) || !is.null(v2$batch_info)) {
        stats::setNames(list(v1$batch_info, v2$batch_info), c(nm1, nm2))
      } else {
        NULL
      }

      emb <- if (!is.null(v1$embedding)) v1$embedding else v2$embedding

      return(list(
        m1 = v1$counts,
        m2 = v2$counts,
        m1_obj = v1$object,
        m2_obj = v2$object,
        m1_name = nm1,
        m2_name = nm2,
        cell_types = ct,
        batch_info = bt,
        embedding = emb
      ))
    }
  }

  if (inherits(obj, "Seurat")) {
    stop(sprintf(
      "scSimEval requires a 'SingleCellExperiment' object (or list of two SingleCellExperiment objects) for '%s' (received 'Seurat').\nPlease convert via Seurat::as.SingleCellExperiment().",
      arg_name
    ), call. = FALSE)
  }

  # Case B: Multi-modal SingleCellExperiment with altExp
  if (inherits(obj, "SingleCellExperiment")) {
    alt_names <- SingleCellExperiment::altExpNames(obj)
    if (length(alt_names) >= 1L) {
      a2 <- SingleCellExperiment::altExp(obj, alt_names[1L])
      m1_nm <- SingleCellExperiment::mainExpName(obj)
      if (is.null(m1_nm) || !nzchar(m1_nm)) m1_nm <- "RNA"

      cd <- SummarizedExperiment::colData(obj)
      # Handle stacked unpaired SCE objects annotated with colData$modality
      if ("modality" %in% colnames(cd)) {
        mod_vals <- toupper(as.character(cd$modality))
        rna_idx <- which(grepl("RNA|GENE|TRANSCR", mod_vals))
        atac_idx <- which(grepl("ATAC|PEAK|CHROMATIN", mod_vals))
        if (length(rna_idx) > 0L && length(atac_idx) > 0L &&
            length(intersect(rna_idx, atac_idx)) == 0L) {
          m1_obj <- obj[, rna_idx]
          m2_obj <- a2[, atac_idx]
          SingleCellExperiment::altExps(m1_obj) <- NULL
          m1 <- SummarizedExperiment::assay(m1_obj, "counts")
          m2 <- SummarizedExperiment::assay(m2_obj, "counts")
          v1 <- validate_sc_object(m1_obj, paste0(arg_name, " (RNA)"), reduction_name = reduction_name,
                                   cell_type_col = cell_type_col, batch_col = batch_col)
          v2 <- validate_sc_object(m2_obj, paste0(arg_name, " (ATAC)"), reduction_name = reduction_name,
                                   cell_type_col = cell_type_col, batch_col = batch_col)
          return(list(
            m1 = m1,
            m2 = m2,
            m1_obj = m1_obj,
            m2_obj = m2_obj,
            m1_name = m1_nm,
            m2_name = alt_names[1L],
            cell_types = stats::setNames(list(v1$cell_types, v2$cell_types), c(m1_nm, alt_names[1L])),
            batch_info = stats::setNames(list(v1$batch_info, v2$batch_info), c(m1_nm, alt_names[1L])),
            embedding = v1$embedding
          ))
        }
      }

      m1 <- SummarizedExperiment::assay(obj, "counts")
      m2 <- SummarizedExperiment::assay(a2, "counts")
      v <- validate_sc_object(obj, arg_name, reduction_name = reduction_name,
                              cell_type_col = cell_type_col, batch_col = batch_col)
      return(list(
        m1 = m1,
        m2 = m2,
        m1_obj = obj,
        m2_obj = a2,
        m1_name = m1_nm,
        m2_name = alt_names[1L],
        cell_types = v$cell_types,
        batch_info = v$batch_info,
        embedding = v$embedding
      ))
    } else {
      stop(sprintf("Multiomics SingleCellExperiment object '%s' must contain at least one altExp for the secondary modality (e.g. ATAC).", arg_name), call. = FALSE)
    }
  }

  stop(sprintf(
    "scSimEval strictly requires multiomics inputs for '%s' to be 'SingleCellExperiment' objects (with altExp or as a named list). Received '%s'.",
    arg_name, paste(class(obj), collapse = "/")
  ), call. = FALSE)
}

.sc_parallel_cache <- new.env(parent = emptyenv())
.sc_parallel_cache$bp <- NULL
.sc_parallel_cache$workers <- 0L
.sc_parallel_cache$target_threads <- 1L

.get_windows_bpparam <- function(workers) {
  target_w <- max(as.integer(workers), as.integer(.sc_parallel_cache$target_threads), 2L)
  bp_cur <- .sc_parallel_cache$bp
  w_cur <- .sc_parallel_cache$workers
  if (!is.null(bp_cur) && !is.null(w_cur) && w_cur >= workers &&
      isTRUE(tryCatch(BiocParallel::bpisup(bp_cur), error = function(e) FALSE))) {
    return(bp_cur)
  }
  if (!is.null(bp_cur)) {
    try(BiocParallel::bpstop(bp_cur), silent = TRUE)
    .sc_parallel_cache$bp <- NULL
    .sc_parallel_cache$workers <- 0L
  }
  bp_new <- BiocParallel::SnowParam(workers = target_w)
  BiocParallel::bpstart(bp_new)
  .sc_parallel_cache$bp <- bp_new
  .sc_parallel_cache$workers <- target_w
  bp_new
}

.set_parallel_threads <- function(threads = 1L) {
  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) {
    1L
  } else {
    max(1L, as.integer(threads[1L]))
  }
  .sc_parallel_cache$target_threads <- threads
  options(mc.cores = threads)
  if (requireNamespace("RcppParallel", quietly = TRUE)) {
    try(RcppParallel::setThreadOptions(numThreads = threads), silent = TRUE)
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(threads), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(threads), silent = TRUE)
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    try(data.table::setDTthreads(threads), silent = TRUE)
  }
  invisible(threads)
}

.sc_parallel_lapply <- function(X, FUN, threads = 1L, ...) {
  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) {
    1L
  } else {
    max(1L, as.integer(threads[1L]))
  }
  n_tasks <- length(X)
  if (threads <= 1L || n_tasks <= 1L) {
    return(lapply(X, FUN, ...))
  }
  if (.Platform$OS.type != "windows") {
    res <- tryCatch({
      mclapply_fn <- getExportedValue("parallel", "mclapply")
      mclapply_fn(X, FUN, ..., mc.cores = min(threads, n_tasks), mc.preschedule = TRUE)
    }, error = function(e) NULL)
    if (!is.null(res) && length(res) == n_tasks &&
        !any(vapply(res, inherits, logical(1), "try-error"))) {
      return(res)
    }
    return(lapply(X, FUN, ...))
  }
  res <- tryCatch(
    BiocParallel::bplapply(
      X, FUN, ...,
      BPPARAM = .get_windows_bpparam(min(threads, n_tasks))
    ),
    error = function(e) NULL
  )
  if (!is.null(res) && length(res) == n_tasks) {
    return(res)
  }
  lapply(X, FUN, ...)
}
