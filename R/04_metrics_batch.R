.knn_coords <- function(coords, k, threads = 1L) {
  coords <- as.matrix(coords)
  n_cells <- nrow(coords)
  if (n_cells <= 1L) {
    return(list(
      idx = matrix(integer(0), nrow = n_cells, ncol = 0L),
      dists = matrix(numeric(0), nrow = n_cells, ncol = 0L)
    ))
  }
  k <- max(1L, min(as.integer(k), n_cells - 1L))
  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) 1L else max(1L, as.integer(threads[1L]))
  n_chunks <- if (threads > 1L && n_cells >= 60L) min(threads, n_cells) else 1L

  if (requireNamespace("RANN", quietly = TRUE)) {
    if (n_chunks > 1L) {
      chunks <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / ceiling(n_cells / n_chunks)))
      part_res <- .sc_parallel_lapply(chunks, function(ch) {
        nn <- RANN::nn2(data = coords, query = coords[ch, , drop = FALSE], k = k + 1L)
        list(
          idx = nn$nn.idx[, -1L, drop = FALSE],
          dists = nn$nn.dists[, -1L, drop = FALSE]
        )
      }, threads = n_chunks)
      return(list(
        idx = do.call(rbind, lapply(part_res, `[[`, "idx")),
        dists = do.call(rbind, lapply(part_res, `[[`, "dists"))
      ))
    }
    nn <- RANN::nn2(coords, k = k + 1L)
    return(list(
      idx = nn$nn.idx[, -1L, drop = FALSE],
      dists = nn$nn.dists[, -1L, drop = FALSE]
    ))
  }

  cross_dist_fn <- fast_cross_distance_sq
  eval_knn_chunk <- function(ch) {
    m <- length(ch)
    idx_part <- matrix(0L, nrow = m, ncol = k)
    dist_part <- matrix(0, nrow = m, ncol = k)
    block_starts <- seq(1L, m, by = 512L)
    for (bs in block_starts) {
      be <- min(m, bs + 511L)
      sub_pos <- bs:be
      ch_sub <- ch[sub_pos]
      d_sub <- sqrt(cross_dist_fn(coords[ch_sub, , drop = FALSE], coords))
      d_sub[cbind(seq_along(ch_sub), ch_sub)] <- Inf
      for (r in seq_along(ch_sub)) {
        nn_r <- order(d_sub[r, ])[seq_len(k)]
        idx_part[sub_pos[r], ] <- nn_r
        dist_part[sub_pos[r], ] <- d_sub[r, nn_r]
      }
    }
    list(idx = idx_part, dists = dist_part)
  }

  if (n_chunks > 1L) {
    chunks <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / ceiling(n_cells / n_chunks)))
    part_res <- .sc_parallel_lapply(chunks, eval_knn_chunk, threads = n_chunks)
    return(list(
      idx = do.call(rbind, lapply(part_res, `[[`, "idx")),
      dists = do.call(rbind, lapply(part_res, `[[`, "dists"))
    ))
  }

  eval_knn_chunk(seq_len(n_cells))
}

.fast_pca_batch <- function(hvg_mat, n_pcs_use) {
  X <- t(hvg_mat)
  min_dim <- min(nrow(X), ncol(X))
  if (n_pcs_use >= 2L && min_dim > 50L && n_pcs_use < 0.5 * min_dim) {
    res <- .run_prcomp_irlba(X, n = n_pcs_use, center = TRUE, scale. = FALSE)
    if (!is.null(res) && !is.null(res$x) && !is.null(res$sdev)) {
      return(res)
    }
  }
  stats::prcomp(X, center = TRUE, scale. = FALSE, rank. = n_pcs_use)
}

.parallel_batch_silhouette <- function(coords, batch_info, threads = 1L, return_scores = FALSE) {
  coords <- as.matrix(coords)
  n_raw <- nrow(coords)
  if (n_raw < 3L) return(if (return_scores) numeric(0) else NA_real_)
  valid <- which(!is.na(batch_info) & stats::complete.cases(coords))
  if (length(valid) < 3L) return(if (return_scores) numeric(0) else NA_real_)
  coords <- coords[valid, , drop = FALSE]
  g_fac <- droplevels(as.factor(batch_info[valid]))
  g <- as.integer(g_fac)
  K <- max(g)
  if (K < 2L) return(if (return_scores) numeric(0) else NA_real_)
  tab <- tabulate(g, nbins = K)
  if (all(tab <= 1L)) return(if (return_scores) numeric(0) else NA_real_)

  n_cells <- nrow(coords)
  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) 1L else max(1L, as.integer(threads[1L]))
  n_chunks <- if (threads > 1L && n_cells >= 60L) min(threads, n_cells) else 1L
  chunks <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / ceiling(n_cells / n_chunks)))

  cross_dist_fn <- fast_cross_distance_sq
  eval_sil_chunk <- function(chunk_idx) {
    m <- length(chunk_idx)
    res <- numeric(m)
    block_starts <- seq(1L, m, by = 512L)
    for (bs in block_starts) {
      be <- min(m, bs + 511L)
      sub_pos <- bs:be
      idx_sub <- chunk_idx[sub_pos]
      d_sub <- sqrt(cross_dist_fn(coords[idx_sub, , drop = FALSE], coords))
      d_sub[cbind(seq_along(idx_sub), idx_sub)] <- 0
      g_sums <- rowsum(t(d_sub), group = g, reorder = TRUE)
      g_means <- g_sums / tab
      for (r in seq_along(idx_sub)) {
        gi <- g[idx_sub[r]]
        if (tab[gi] <= 1L) {
          res[sub_pos[r]] <- 0
        } else {
          a_i <- g_sums[gi, r] / (tab[gi] - 1L)
          b_i <- min(g_means[-gi, r])
          denom <- max(a_i, b_i)
          res[sub_pos[r]] <- if (denom > 0) (b_i - a_i) / denom else 0
        }
      }
    }
    res
  }

  sil_vals <- unlist(.sc_parallel_lapply(chunks, eval_sil_chunk, threads = n_chunks), use.names = FALSE)
  if (return_scores) return(sil_vals)
  mean(sil_vals, na.rm = TRUE)
}

#' Cell-Specific Mixing Score (CMS)
#'
#' Tests the hypothesis that group-specific distance distributions of k-nearest
#' neighbor cells have the same underlying unspecified distribution.
#' Computes cell-specific p-values via Anderson-Darling tests (if kSamples
#' installed) or two-sample Kolmogorov-Smirnov / Kruskal-Wallis tests, with
#' optional delegation to CellMixS::cms.
#'
#' @param data Count or log-normalized expression matrix (features x cells) or NULL if coords provided.
#' @param batch_info Factor or character vector of batch assignments for each cell.
#' @param k Neighborhood size k for k-NN graph. Default is min(table(batch_info))/2 or 30.
#' @param n_pcs Number of principal components to compute if coords is NULL. Default is 30.
#' @param coords Optional precomputed embedding matrix (cells x dimensions, e.g. PCA or UMAP).
#' @param cell_min Minimum number of cells per batch in neighborhood to perform test. Default is 4.
#' @param unbalanced Logical, whether to set pure neighborhoods to NA. Default is FALSE.
#' @param threads Number of CPU threads (default 1).
#' @param knn_res Optional precomputed k-NN list (`idx`, `dists`) from `.knn_coords`.
#'
#' @return A named list:
#'   \item{cms_scores}{Numeric vector of CMS p-values for each cell (higher means better mixed)}
#'   \item{mean_cms}{Mean CMS score across cells}
#'   \item{median_cms}{Median CMS score across cells}
#' @examples
#' coords <- matrix(stats::rnorm(100), 50, 2)
#' batch <- factor(rep(c("B1", "B2"), length.out = 50))
#' calc_cms(coords, batch)
#' @keywords internal
#' @export
calc_cms <- function(
  data = NULL,
  batch_info,
  k = NULL,
  n_pcs = 30,
  coords = NULL,
  cell_min = 4,
  unbalanced = FALSE,
  threads = 1,
  knn_res = NULL
) {
  batch_info <- as.factor(batch_info)
  
  if (is.null(coords)) {
    d_use <- normalize_counts(data)
    d_use <- if (inherits(d_use, "Matrix")) as.matrix(d_use) else d_use
    coords <- fast_pca(t(d_use), n_pcs = n_pcs, n_hvg = 2000, scale = FALSE)
  }
  
  n_cells <- nrow(coords)
  if (is.null(k)) {
    k <- max(10, round(min(table(batch_info)) / 2))
  }
  k <- min(k, n_cells - 1)
  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) 1L else max(1L, as.integer(threads[1L]))

  # Delegate to CellMixS::cms only for single-threaded uncached calls;
  # multithreaded runs use the chunk-parallel implementation to avoid
  # single-threaded internal k-NN and per-cell IPC overhead.
  if (threads <= 1L && is.null(knn_res) &&
      requireNamespace("CellMixS", quietly = TRUE) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE) &&
      requireNamespace("BiocParallel", quietly = TRUE)) {
    res_cms <- tryCatch({
      bpparam <- if (threads > 1L) {
        if (.Platform$OS.type == "windows") {
          BiocParallel::SnowParam(workers = threads)
        } else {
          BiocParallel::MulticoreParam(workers = threads)
        }
      } else {
        BiocParallel::SerialParam()
      }
      dummy_counts <- matrix(0, nrow = 2L, ncol = n_cells)
      if (!is.null(rownames(coords)) && length(rownames(coords)) == n_cells) {
        colnames(dummy_counts) <- rownames(coords)
      }
      sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = dummy_counts),
        colData = data.frame(batch = batch_info)
      )
      SingleCellExperiment::reducedDim(sce, "PCA") <- coords
      eval_res <- CellMixS::cms(
        sce, k = k, group = "batch", dim_red = "PCA",
        n_dim = ncol(coords), BPPARAM = bpparam
      )
      p_vals <- eval_res$cms
      return(list(
        cms_scores = p_vals,
        mean_cms = mean(p_vals, na.rm = TRUE),
        median_cms = stats::median(p_vals, na.rm = TRUE)
      ))
    }, error = function(e) NULL)
    if (!is.null(res_cms)) return(res_cms)
  }

  # Native parallel implementation (parallel k-NN + parallel cell chunks)
  if (is.null(knn_res)) {
    knn_res <- .knn_coords(coords, k, threads = threads)
  }
  nn_idx_mat <- knn_res$idx
  nn_dist_mat <- knn_res$dists
  has_ksamples <- requireNamespace("kSamples", quietly = TRUE)

  eval_cms_cells <- function(cell_indices) {
    vapply(cell_indices, function(i) {
      nn_idx <- nn_idx_mat[i, ]
      nn_dists <- nn_dist_mat[i, ]
      nn_batches <- batch_info[nn_idx]

      dist_by_b <- split(nn_dists, nn_batches)
      dist_by_b <- dist_by_b[vapply(dist_by_b, length, integer(1)) >= cell_min]

      if (length(dist_by_b) <= 1) {
        return(if (unbalanced) NA_real_ else 0)
      }

      if (has_ksamples) {
        ad_res <- tryCatch(kSamples::ad.test(dist_by_b), error = function(e) NULL)
        if (!is.null(ad_res)) {
          return(mean(ad_res$ad[, " asympt. P-value"], na.rm = TRUE))
        }
      }

      if (length(dist_by_b) == 2) {
        ks_res <- tryCatch(stats::ks.test(dist_by_b[[1]], dist_by_b[[2]]), error = function(e) NULL)
        return(if (!is.null(ks_res)) ks_res$p.value else NA_real_)
      } else {
        flat_dists <- unlist(dist_by_b, use.names = FALSE)
        flat_groups <- factor(rep(names(dist_by_b), lengths(dist_by_b)))
        kw_res <- tryCatch(stats::kruskal.test(flat_dists ~ flat_groups), error = function(e) NULL)
        return(if (!is.null(kw_res)) kw_res$p.value else NA_real_)
      }
    }, numeric(1))
  }

  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) 1L else max(1L, as.integer(threads[1L]))
  if (threads > 1L && n_cells >= 30L) {
    n_chunks <- min(threads, n_cells)
    chunks <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / ceiling(n_cells / n_chunks)))
    p_vals <- unlist(.sc_parallel_lapply(chunks, eval_cms_cells, threads = n_chunks), use.names = FALSE)
  } else {
    p_vals <- eval_cms_cells(seq_len(n_cells))
  }
  
  list(
    cms_scores = p_vals,
    mean_cms = mean(p_vals, na.rm = TRUE),
    median_cms = stats::median(p_vals, na.rm = TRUE)
  )
}

#' Inverse Simpson Index (ISI) for Batch Mixing
#'
#' Evaluates the diversity of batch labels within each cell's k-nearest neighborhood
#' using the Inverse Simpson Index (1 / sum(p_b^2)).
#' Implements a pure-R, lightweight formulation of LISI with optional distance weighting.
#'
#' @param coords Matrix of cell coordinates (cells x dimensions).
#' @param batch_info Factor or character vector of batch assignments.
#' @param k Number of nearest neighbors. Default is 30.
#' @param weighted Logical, whether to weight neighbor contributions by inverse distance. Default is TRUE.
#' @param threads Number of CPU threads (default 1).
#' @param knn_res Optional precomputed k-NN list (`idx`, `dists`) from `.knn_coords`.
#'
#' @return A named list:
#'   \item{isi_scores}{Numeric vector of cell-level ISI values}
#'   \item{mean_isi}{Mean ISI score (ranges from 1 to number of batches; higher indicates better mixing)}
#'   \item{median_isi}{Median ISI score}
#' @examples
#' coords <- matrix(stats::rnorm(100), 50, 2)
#' batch <- factor(rep(c("B1", "B2"), length.out = 50))
#' calc_isi(coords, batch)
#' @keywords internal
#' @export
calc_isi <- function(coords, batch_info, k = 30, weighted = TRUE, threads = 1, knn_res = NULL) {
  batch_info <- as.factor(batch_info)
  coords <- as.matrix(coords)
  n_cells <- nrow(coords)
  k <- min(k, n_cells - 1)
  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) 1L else max(1L, as.integer(threads[1L]))

  if (is.null(knn_res)) {
    knn_res <- .knn_coords(coords, k, threads = threads)
  }
  nn_idx_mat <- knn_res$idx
  nn_dist_mat <- knn_res$dists
  batch_levels <- levels(batch_info)
  n_batch <- length(batch_levels)

  eval_isi_cells <- function(cell_indices) {
    vapply(cell_indices, function(i) {
      nn_idx <- nn_idx_mat[i, ]
      nn_b <- batch_info[nn_idx]
      nn_d <- nn_dist_mat[i, ]

      if (weighted) {
        weights <- 1 / (nn_d + 1)
        p_b <- vapply(batch_levels, function(b) {
          sum(weights[nn_b == b]) / sum(weights)
        }, numeric(1))
      } else {
        p_b <- vapply(batch_levels, function(b) {
          sum(nn_b == b) / k
        }, numeric(1))
      }

      simpson <- sum(p_b^2)
      if (simpson <= 0) as.numeric(n_batch) else 1 / simpson
    }, numeric(1))
  }

  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) 1L else max(1L, as.integer(threads[1L]))
  if (threads > 1L && n_cells >= 30L) {
    n_chunks <- min(threads, n_cells)
    chunks <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / ceiling(n_cells / n_chunks)))
    isi_vals <- unlist(.sc_parallel_lapply(chunks, eval_isi_cells, threads = n_chunks), use.names = FALSE)
  } else {
    isi_vals <- eval_isi_cells(seq_len(n_cells))
  }

  list(
    isi_scores = isi_vals,
    mean_isi = mean(isi_vals, na.rm = TRUE),
    median_isi = stats::median(isi_vals, na.rm = TRUE)
  )
}

#' Calculate Batch Integration Metrics
#'
#' Evaluates how realistically simulated technical batches mix or separate,
#' computing the five canonical Category IV (Technical Variation & Batch Mixing
#' Dynamics) metrics: Batch Silhouette, PC Regression R2, CMS (Cell-Specific
#' Mixing Score), Inverse Simpson Index (ISI), and cross-batch cell type
#' classification macro F1.
#'
#' @param data Count or log-normalized expression matrix (features x cells).
#' @param batch_info Factor or character vector of batch assignments for each cell.
#' @param k Neighborhood size k for k-NN graph and mixing metrics. Default is min(table(batch_info))/2.
#' @param n_pcs Number of principal components for embedding evaluation. Default is 30.
#' @param cell_types Optional factor or character vector of cell types for cross-batch evaluation.
#' @param threads Number of CPU threads (default 1).
#' @param verbose Logical, whether to print progress.
#'
#' @return A named list of batch metrics:
#'   \item{batch_silhouette}{Average silhouette width using batch labels (lower is better mixed)}
#'   \item{pcr_r2}{Total variance explained by batch in principal components (PCR)}
#'   \item{cms}{Mean Cell-specific mixing score (CellMixS)}
#'   \item{isi}{Mean Inverse Simpson Index for batch mixing (CellMixS)}
#'   \item{cross_batch_F1}{Mean cross-batch cell type transfer macro F1 (if cell_types supplied)}
#' @examples
#' data(example_scrna, package = "scSimEval")
#' evaluate_batch_metrics(example_scrna$ref, example_scrna$batch_info,
#'                        cell_types = example_scrna$cell_types)
#' @export
evaluate_batch_metrics <- function(
  data,
  batch_info,
  k = NULL,
  n_pcs = 30,
  cell_types = NULL,
  threads = 1,
  verbose = FALSE
) {
  threads <- .set_parallel_threads(threads)
  batch_info <- as.factor(batch_info)
  if (inherits(data, "SingleCellExperiment")) {
    data <- SummarizedExperiment::assay(data, "counts")
  } else if (inherits(data, "Seurat")) {
    stop("Seurat objects are not directly supported; convert via Seurat::as.SingleCellExperiment().", call. = FALSE)
  }
  n_cells <- ncol(data)
  
  if (is.null(k)) {
    k <- max(5, round(min(table(batch_info)) / 2))
  }
  k <- max(1L, min(as.integer(k), n_cells - 1L))
  
  if (verbose) message("Running PCA for batch evaluation...")
  data_norm <- normalize_counts(data)
  gene_vars <- fast_row_vars(data_norm)
  n_hvg <- min(2000L, length(gene_vars))
  hvg_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_hvg)]
  hvg_mat <- as.matrix(data_norm[hvg_idx, , drop = FALSE])
  n_pcs_use <- max(1L, min(as.integer(n_pcs), n_cells - 1L, nrow(hvg_mat) - 1L))
  pca_res <- .fast_pca_batch(hvg_mat, n_pcs_use)
  pca_coords <- pca_res$x[, seq_len(min(n_pcs_use, ncol(pca_res$x))), drop = FALSE]
  
  # 1. Batch Silhouette Width (parallel Euclidean silhouette, with kBET/cluster fallback)
  if (verbose) message("Calculating batch silhouette width...")
  batch_sil <- if (threads > 1L && n_cells >= 60L) {
    .parallel_batch_silhouette(pca_coords, batch_info, threads = threads)
  } else {
    NA_real_
  }
  pkg_kbet <- "kBET"
  if (is.na(batch_sil)) {
    if (requireNamespace(pkg_kbet, quietly = TRUE)) {
      batch_sil <- tryCatch({
        kbet_batch_sil <- getExportedValue(pkg_kbet, "batch_sil")
        kbet_batch_sil(pca.data = pca_res, batch = batch_info, nPCs = min(n_pcs, ncol(pca_coords)))
      }, error = function(e) NA_real_)
    }
    if (is.na(batch_sil)) {
      batch_sil <- .parallel_batch_silhouette(pca_coords, batch_info, threads = threads)
    }
    if (is.na(batch_sil) && requireNamespace("cluster", quietly = TRUE)) {
      dist_pca <- stats::dist(pca_coords)
      sil_res <- tryCatch(cluster::silhouette(as.numeric(as.factor(batch_info)), dist_pca), error = function(e) NULL)
      if (!is.null(sil_res) && (is.matrix(sil_res) || inherits(sil_res, "silhouette"))) {
        batch_sil <- mean(sil_res[, 3], na.rm = TRUE)
      }
    }
  }
  
  # 2. Principal Component Regression (PCR) R2
  if (verbose) message("Calculating principal component regression (PCR)...")
  pcr_r2 <- NA_real_
  if (threads <= 1L && requireNamespace(pkg_kbet, quietly = TRUE)) {
    pcr_res <- tryCatch({
      kbet_pcr <- getExportedValue(pkg_kbet, "pcRegression")
      kbet_pcr(pca.data = pca_res, batch = as.character(batch_info), n_top = min(n_pcs, ncol(pca_coords)))
    }, error = function(e) NULL)
    if (!is.null(pcr_res)) pcr_r2 <- pcr_res$R2Var
  }
  if (is.na(pcr_r2)) {
    n_pc_cols <- ncol(pca_coords)
    mm <- stats::model.matrix(~ batch_info)
    fit_mv <- stats::lm.fit(mm, pca_coords)
    y_c <- sweep(pca_coords, 2L, colMeans(pca_coords), "-")
    y_hat_c <- sweep(as.matrix(fit_mv$fitted.values), 2L, colMeans(pca_coords), "-")
    denom_ss <- colSums(y_c^2)
    r2_vals <- ifelse(denom_ss > 0, colSums(y_hat_c^2) / denom_ss, 0)
    pc_vars <- pca_res$sdev[seq_len(n_pc_cols)]^2
    pcr_r2 <- sum(r2_vals * pc_vars) / sum(pc_vars)
  }

  # Precompute k-NN graph once in parallel and share across CMS and ISI
  knn_res <- .knn_coords(pca_coords, k, threads = threads)

  # 3. Cell-Specific Mixing Score (CMS)
  if (verbose) message("Calculating CMS (CellMixS)...")
  cms_res <- calc_cms(data = NULL, batch_info = batch_info, k = k, coords = pca_coords, threads = threads, knn_res = knn_res)
  cms_score <- cms_res$mean_cms

  # 4. Inverse Simpson Index (ISI)
  if (verbose) message("Calculating Inverse Simpson Index (ISI)...")
  isi_res <- calc_isi(coords = pca_coords, batch_info = batch_info, k = k, weighted = TRUE, threads = threads, knn_res = knn_res)
  isi_score <- isi_res$mean_isi
  
  # 5. Cross-Batch Cell Type Classification Macro F1
  cross_batch_f1 <- NA_real_
  if (!is.null(cell_types)) {
    if (verbose) message("Evaluating cross-batch cell type classification...")
    cb_res <- evaluate_cross_batch_prediction(
      coords = pca_coords,
      cell_types = cell_types,
      batch_info = batch_info,
      k = min(5, k),
      threads = threads
    )
    cross_batch_f1 <- cb_res$mean_cross_batch_F1
  }
  
  list(
    batch_silhouette = batch_sil,
    pcr_r2 = pcr_r2,
    cms = cms_score,
    isi = isi_score,
    cross_batch_F1 = cross_batch_f1
  )
}

#' Evaluate Cross-Batch Cell Classification Transfer F1
#'
#' Evaluates whether simulated or integrated multi-batch single-cell data preserve
#' cell-type identity across distinct experimental batches or donors.
#' Trains a k-NN classifier on each individual batch and tests the macro F1 score
#' of predictions on all other remaining batches.
#'
#' @param coords Embedding coordinates (cells x dimensions).
#' @param cell_types Factor or character vector of ground-truth cell type labels.
#' @param batch_info Factor or character vector of batch assignments.
#' @param k Integer number of nearest neighbors for the k-NN classifier (default 5).
#' @param threads Number of CPU threads (default 1).
#'
#' @return A list containing:
#'   \item{mean_cross_batch_F1}{Overall mean macro-averaged F1 score across all directed batch pairs}
#' @examples
#' coords <- matrix(stats::rnorm(100), 50, 2)
#' batch <- factor(rep(c("B1", "B2"), length.out = 50))
#' ct <- factor(rep(c("TypeA", "TypeB"), length.out = 50))
#' evaluate_cross_batch_prediction(coords, ct, batch)
#' @keywords internal
#' @export
evaluate_cross_batch_prediction <- function(coords, cell_types, batch_info, k = 5, threads = 1) {
  cell_types <- as.factor(cell_types)
  batch_info <- as.factor(batch_info)
  coords <- as.matrix(coords)
  
  batches <- levels(batch_info)
  n_b <- length(batches)
  pair_f1 <- matrix(NA_real_, nrow = n_b, ncol = n_b)
  cls <- levels(cell_types)
  
  pair_list <- list()
  for (i in seq_along(batches)) {
    train_idx <- which(batch_info == batches[i])
    train_y <- cell_types[train_idx]
    if (length(unique(train_y)) < 2) next
    for (j in seq_along(batches)) {
      if (i == j) next
      test_idx <- which(batch_info == batches[j])
      if (length(test_idx) == 0L) next
      k_use <- max(1L, min(k, length(train_idx) - 1L))
      pair_list[[length(pair_list) + 1L]] <- list(
        i = i, j = j, train_idx = train_idx, test_idx = test_idx, k_use = k_use
      )
    }
  }
  if (length(pair_list) == 0L) {
    return(list(mean_cross_batch_F1 = NA_real_))
  }

  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) 1L else max(1L, as.integer(threads[1L]))
  chunks_per_pair <- if (threads > 1L) max(1L, as.integer(ceiling(threads / length(pair_list)))) else 1L

  tasks <- list()
  for (p_idx in seq_along(pair_list)) {
    t_idx <- pair_list[[p_idx]]$test_idx
    n_ch <- min(chunks_per_pair, length(t_idx))
    if (n_ch > 1L) {
      grps <- split(seq_along(t_idx), ceiling(seq_along(t_idx) / ceiling(length(t_idx) / n_ch)))
      for (g_pos in grps) {
        tasks[[length(tasks) + 1L]] <- list(p_idx = p_idx, pos = g_pos, test_idx = t_idx[g_pos])
      }
    } else {
      tasks[[length(tasks) + 1L]] <- list(p_idx = p_idx, pos = seq_along(t_idx), test_idx = t_idx)
    }
  }

  task_preds <- .sc_parallel_lapply(tasks, function(tk) {
    pr <- pair_list[[tk$p_idx]]
    as.character(class::knn(
      train = coords[pr$train_idx, , drop = FALSE],
      test = coords[tk$test_idx, , drop = FALSE],
      cl = cell_types[pr$train_idx],
      k = pr$k_use
    ))
  }, threads = threads)

  preds_by_pair <- lapply(pair_list, function(pr) character(length(pr$test_idx)))
  for (tk_i in seq_along(tasks)) {
    tk <- tasks[[tk_i]]
    preds_by_pair[[tk$p_idx]][tk$pos] <- task_preds[[tk_i]]
  }

  for (p_idx in seq_along(pair_list)) {
    pr <- pair_list[[p_idx]]
    pred_y <- preds_by_pair[[p_idx]]
    test_y <- cell_types[pr$test_idx]
    f1s <- vapply(cls, function(c) {
      tp <- sum(pred_y == c & test_y == c)
      fp <- sum(pred_y == c & test_y != c)
      fn <- sum(pred_y != c & test_y == c)
      if ((2 * tp + fp + fn) > 0) (2 * tp) / (2 * tp + fp + fn) else 0
    }, numeric(1))
    pair_f1[pr$i, pr$j] <- mean(f1s[cls %in% test_y], na.rm = TRUE)
  }
  
  list(
    mean_cross_batch_F1 = mean(pair_f1, na.rm = TRUE)
  )
}
