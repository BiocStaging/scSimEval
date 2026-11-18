# Seurat count extraction compatibility helpers
#
# Keep all Seurat count extraction behind these helpers. Seurat v5 Assay5
# stores expression data in layers (including split counts.* layers), not in
# the legacy @counts slot used by Assay.

.sanitize_signac_assays <- function(obj) {
  if (!isS4(obj) || !"Seurat" %in% as.character(class(obj))) return(obj)
  if (requireNamespace("Signac", quietly = TRUE)) return(obj)
  for (nm in names(obj@assays)) {
    a <- obj@assays[[nm]]
    pkg <- attr(class(a), "package")
    if (!is.null(pkg) && identical(pkg, "Signac") && methods::.hasSlot(a, "counts")) {
      cts <- methods::slot(a, "counts")
      new_a <- tryCatch(SeuratObject::CreateAssayObject(counts = cts), error = function(e) NULL)
      if (!is.null(new_a)) {
        if (methods::.hasSlot(a, "data") && length(dim(methods::slot(a, "data"))) == 2L &&
            nrow(methods::slot(a, "data")) == nrow(cts)) {
          methods::slot(new_a, "data") <- methods::slot(a, "data")
        }
        if (methods::.hasSlot(a, "var.features")) {
          methods::slot(new_a, "var.features") <- methods::slot(a, "var.features")
        }
        if (methods::.hasSlot(a, "key") && nzchar(methods::slot(a, "key"))) {
          methods::slot(new_a, "key") <- methods::slot(a, "key")
        }
        obj@assays[[nm]] <- new_a
      }
    }
  }
  obj
}

.seurat_assay_counts <- function(obj, assay_name = NULL) {
  obj <- .sanitize_signac_assays(obj)
  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object to extract assay counts.", call. = FALSE)
  }
  if (is.null(assay_name) || length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
    assay_name <- .seurat_default_assay(obj)
  }
  if (is.null(assay_name) || !assay_name %in% names(obj@assays)) {
    stop("The requested Seurat assay is not present in the object.", call. = FALSE)
  }

  assay_obj <- obj@assays[[assay_name]]
  if (is.null(assay_obj)) {
    assay_obj <- tryCatch(obj[[assay_name]], error = function(e) NULL)
  }
  if (is.null(assay_obj)) {
    stop(sprintf("Seurat assay '%s' could not be retrieved.", assay_name), call. = FALSE)
  }

  # Assay5 may contain one counts layer or multiple split layers (counts.*).
  # Join split layers on a local assay copy; do not modify the caller's object.
  if (inherits(assay_obj, "Assay5")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("SeuratObject is required to read counts from a Seurat v5 Assay5.", call. = FALSE)
    }

    layer_names <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character())
    count_layers <- grep("^counts($|\\.)", layer_names, value = TRUE)
    if (length(count_layers) == 0L) {
      stop(sprintf("Seurat assay '%s' has no raw counts layer.", assay_name), call. = FALSE)
    }

    assay_for_read <- assay_obj
    if (length(count_layers) > 1L || !identical(count_layers, "counts")) {
      assay_for_read <- tryCatch(
        SeuratObject::JoinLayers(assay_obj, layers = count_layers, new = "counts"),
        error = function(e) {
          stop(sprintf("Could not join split counts layers in Seurat assay '%s': %s", assay_name, conditionMessage(e)), call. = FALSE)
        }
      )
    }

    counts <- tryCatch(
      SeuratObject::LayerData(assay_for_read, layer = "counts"),
      error = function(e) {
        tryCatch(
          Seurat::GetAssayData(assay_for_read, layer = "counts"),
          error = function(e2) NULL
        )
      }
    )
    if (!is.null(counts) && length(dim(counts)) == 2L && nrow(counts) > 0L) {
      return(counts)
    }
    stop(sprintf("Could not retrieve the counts layer from Seurat assay '%s'.", assay_name), call. = FALSE)
  }

  # Direct S4 @counts slot extraction (handles Seurat v3/v4 Assay and Signac ChromatinAssay without requiring Signac)
  if (methods::.hasSlot(assay_obj, "counts")) {
    counts <- tryCatch(methods::slot(assay_obj, "counts"), error = function(e) NULL)
    if (!is.null(counts) && length(dim(counts)) == 2L && nrow(counts) > 0L) {
      return(counts)
    }
  }

  # Assay-level layer API fallback
  counts <- tryCatch(
    Seurat::GetAssayData(assay_obj, layer = "counts"),
    error = function(e) NULL
  )
  if (!is.null(counts) && length(dim(counts)) == 2L && nrow(counts) > 0L) {
    return(counts)
  }

  # Object-level layer/slot fallbacks
  counts <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay_name, layer = "counts"),
    error = function(e) {
      tryCatch(
        Seurat::GetAssayData(obj, assay = assay_name, slot = "counts"),
        error = function(e2) NULL
      )
    }
  )
  if (!is.null(counts) && length(dim(counts)) == 2L && nrow(counts) > 0L) {
    return(counts)
  }

  stop(sprintf("Could not retrieve raw counts from Seurat assay '%s'.", assay_name), call. = FALSE)
}

.get_seurat_counts <- function(obj) {
  .seurat_assay_counts(obj, .seurat_default_assay(obj))
}
