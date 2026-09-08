#' Diagnose patch assignments
#'
#' Report patch size, variation in `X`, spatial connectivity, assignment
#' coverage, and final membership stability when iteration logs are available.
#'
#' @param xy Numeric matrix with cells in rows and x/y coordinates in columns.
#' @param X Numeric vector or matrix with cells in rows.
#' @param patches Final patch vector or the list returned by
#'   `getPatches(log_iters = TRUE)`. `NA` denotes an unassigned cell. A vector
#'   gives `NA` stability; a logged result compares the last two iterations and
#'   assumes patch identifiers remain stable.
#' @param k Maximum number of spatial neighbors. Default 10, matching
#'   `getPatches()`. Limited to `n - 1` for small datasets.
#' @param strict_k Neighbors used for `strict_component_fraction`. Must be no
#'   greater than `k`. `NULL` uses `min(5, k)`.
#'
#' @details
#' The metrics describe different aspects of patch quality.
#'
#' \describe{
#'   \item{`n_cells`}{Number of cells assigned to the patch.}
#'   \item{`x_sd` and `x_n`}{`x_sd` measures variation in a design variable
#'     within the patch. Values near zero indicate little within-patch contrast
#'     for estimating its effect. Larger values indicate more contrast, but the
#'     magnitude depends on the scale of `X`. `x_n`
#'     is the number of finite values used.
#'     The SD is `NA` with fewer than two finite values or if the patch contains
#'     `NaN` or infinite values.}
#'   \item{`strict_component_fraction`}{Fraction of patch cells in the largest
#'     connected component at `strict_k`. One indicates full connectivity;
#'     smaller values indicate greater fragmentation.}
#'   \item{`min_connectivity_k`}{Smallest evaluated k at which the patch is
#'     fully connected. Smaller values indicate stronger local connectivity.
#'     It is `NA` if full connectivity is not reached by `k`, and zero for a
#'     single-cell patch.}
#'   \item{`connectivity_curve`}{For each tested k, reports the fraction of
#'     patch cells in the largest connected piece. One means that all patch
#'     cells are connected. If the value becomes one only as k increases, more
#'     neighbor links are needed to join the patch.}
#'   \item{`membership_stability`}{Fraction of final patch cells with the same
#'     assignment in the preceding iteration. Values closer to one indicate
#'     greater stability. It is `NA` without at least two logged iterations.}
#'   \item{`fraction_assigned` and `fraction_unassigned`}{Proportions of
#'     analyzed cells with and without a patch assignment. A higher
#'     `fraction_assigned` means broader coverage; a higher
#'     `fraction_unassigned` means more cells were excluded. }
#' }
#'
#' Connectivity uses symmetric spatial k-nearest-neighbor graphs. `strict_k`
#' sets the neighborhood size for `strict_component_fraction`; smaller values
#' use fewer links and give a stricter local check. `k` is the largest value
#' tested by `connectivity_curve` and `min_connectivity_k`. Larger k makes
#' connectivity easier, so compare strict fractions only at the same `strict_k`.
#' Effective values are limited to `n - 1` and reported in
#' `assignment_summary`.
#'
#' @return A list with three data frames:
#'   \describe{
#'     \item{patch_diagnostics}{One row per patch. A single `X` gives `x_sd`
#'       and `x_n`; multiple columns give `x_sd_<name>` and `x_n_<name>`.
#'       Other columns report patch size, strict connectivity, the first k with
#'       full connectivity, and membership stability.}
#'     \item{assignment_summary}{Cell assignment counts and effective k values.}
#'     \item{connectivity_curve}{`patch`, `k`, and `component_fraction`, where
#'       the fraction is the largest connected component divided by patch size.}
#'   }
#'   `patch_diagnostics` can be passed to `getPatchPolys()` or
#'   `runInteractivePlotter()` as `patch_data`.
#'
#' @examples
#' ## Use a subset of the example CRC data distributed with SpaceMosaic.
#' mini <- readRDS(system.file("extdata", "miniCRC.RDS", package = "SpaceMosaic"))
#' epithelial_cells <- which(mini$clust == "epi")
#' use <- epithelial_cells[seq_len(min(200L, length(epithelial_cells)))]
#' xy <- mini$xy[use, 1:2, drop = FALSE]
#' X <- mini$X[use]
#'
#' set.seed(1)
#' patch_result <- getPatches(
#'   xy = xy,
#'   X = X,
#'   npatches = 5,
#'   n_iters = 3,
#'   log_iters = TRUE,
#'   verbose = FALSE
#' )
#'
#' res <- getPatchDiagnostics(
#'   xy = xy,
#'   X = X,
#'   patches = patch_result
#' )
#' res$patch_diagnostics
#' res$assignment_summary
#'
#' @export
getPatchDiagnostics <- function(xy, X, patches, k = 10L, strict_k = NULL) {
  if (!is.matrix(xy) || !is.numeric(xy) || ncol(xy) != 2L) {
    stop("xy must be a numeric matrix with exactly two columns.", call. = FALSE)
  }
  if (nrow(xy) < 1L) {
    stop("xy must contain at least one cell.", call. = FALSE)
  }
  if (any(!is.finite(xy))) {
    stop("xy must contain only finite coordinates.", call. = FALSE)
  }

  n <- nrow(xy)
  X <- .prepareDiagnosticX(X, n = n, xy_names = rownames(xy))

  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) ||
      k < 1 || k != floor(k)) {
    stop("k must be a positive integer.", call. = FALSE)
  }
  k <- as.integer(k)
  if (is.null(strict_k)) strict_k <- min(5L, k)
  if (length(strict_k) != 1L || !is.numeric(strict_k) ||
      !is.finite(strict_k) || strict_k < 1 ||
      strict_k != floor(strict_k) || strict_k > k) {
    stop("strict_k must be a positive integer no greater than k.",
         call. = FALSE)
  }
  strict_k <- as.integer(strict_k)

  result_info <- NULL
  if (is.null(patches)) {
    stop(
      "patches must be a patch vector or the output of getPatches().",
      call. = FALSE
    )
  } else if (is.list(patches)) {
    result_info <- .prepareDiagnosticResult(
      patches,
      n = n,
      xy_names = rownames(xy)
    )
    patch <- result_info$patch
  } else {
    .validateDiagnosticPatch(
      patches,
      n = n,
      xy_names = rownames(xy),
      label = "patches"
    )
    patch <- patches
  }

  patch <- as.character(patch)
  patch_ids <- unique(patch[!is.na(patch)])
  numeric_patch_ids <- suppressWarnings(as.numeric(patch_ids))
  patch_ids <- if (all(is.finite(numeric_patch_ids))) {
    patch_ids[order(numeric_patch_ids)]
  } else {
    sort(patch_ids)
  }

  assigned <- !is.na(patch)
  cell_indices <- split(
    which(assigned),
    factor(patch[assigned], levels = patch_ids)
  )
  n_assigned <- sum(assigned)
  effective_k <- if (n >= 2L) min(k, n - 1L) else 0L
  effective_strict_k <- if (effective_k > 0L) {
    min(strict_k, effective_k)
  } else {
    0L
  }
  assignment_summary <- data.frame(
    n_analyzed_cells = n,
    n_assigned_cells = n_assigned,
    fraction_assigned = n_assigned / n,
    n_unassigned_cells = n - n_assigned,
    fraction_unassigned = (n - n_assigned) / n,
    n_nonempty_patches = length(patch_ids),
    connectivity_k = effective_k,
    strict_connectivity_k = effective_strict_k,
    stringsAsFactors = FALSE
  )

  if (ncol(X) == 1L) {
    x_sd_names <- "x_sd"
    x_n_names <- "x_n"
  } else {
    x_sd_names <- paste0("x_sd_", colnames(X))
    x_n_names <- paste0("x_n_", colnames(X))
  }

  if (length(patch_ids) == 0L) {
    patch_diagnostics <- data.frame(
      patch = character(),
      n_cells = integer(),
      stringsAsFactors = FALSE
    )
    for (column in x_sd_names) patch_diagnostics[[column]] <- numeric()
    for (column in x_n_names) {
      patch_diagnostics[[column]] <- integer()
    }
    patch_diagnostics$strict_component_fraction <- numeric()
    patch_diagnostics$min_connectivity_k <- integer()
    patch_diagnostics$membership_stability <- numeric()
    return(list(
      patch_diagnostics = patch_diagnostics,
      assignment_summary = assignment_summary,
      connectivity_curve = data.frame(
        patch = character(),
        k = integer(),
        component_fraction = numeric(),
        stringsAsFactors = FALSE
      )
    ))
  }

  connectivity <- .calculatePatchConnectivity(
    xy = xy,
    patch = patch,
    cell_indices = cell_indices,
    max_k = effective_k,
    strict_k = effective_strict_k
  )

  previous_membership <- NULL
  if (!is.null(result_info) && ncol(result_info$membership_log) >= 2L) {
    previous_membership <- as.character(
      result_info$membership_log[, ncol(result_info$membership_log) - 1L]
    )
  }

  rows <- vector("list", length(patch_ids))
  for (i in seq_along(patch_ids)) {
    patch_id <- patch_ids[i]
    cell_index <- cell_indices[[i]]

    x_sds <- numeric(ncol(X))
    x_n <- integer(ncol(X))
    for (j in seq_len(ncol(X))) {
      values <- X[cell_index, j]
      finite <- is.finite(values)
      x_n[j] <- sum(finite)
      if (any(is.nan(values) | is.infinite(values)) ||
          x_n[j] < 2L) {
        x_sds[j] <- NA_real_
      } else {
        x_sds[j] <- stats::sd(values[finite])
      }
    }

    membership_stability <- NA_real_
    if (!is.null(previous_membership)) {
      retained <- !is.na(previous_membership[cell_index]) &
        previous_membership[cell_index] == patch_id
      membership_stability <- sum(retained) / length(cell_index)
    }

    row <- data.frame(
      patch = patch_id,
      n_cells = length(cell_index),
      stringsAsFactors = FALSE
    )
    for (j in seq_along(x_sd_names)) row[[x_sd_names[j]]] <- x_sds[j]
    for (j in seq_along(x_n_names)) {
      row[[x_n_names[j]]] <- x_n[j]
    }
    row$strict_component_fraction <-
      connectivity$strict_component_fraction[i]
    row$min_connectivity_k <- connectivity$min_connectivity_k[i]
    row$membership_stability <- membership_stability
    rows[[i]] <- row
  }

  patch_diagnostics <- do.call(rbind, rows)
  rownames(patch_diagnostics) <- NULL

  list(
    patch_diagnostics = patch_diagnostics,
    assignment_summary = assignment_summary,
    connectivity_curve = connectivity$connectivity_curve
  )
}


#' Validate and format X
#'
#' @param X Numeric vector or matrix of design variables.
#' @param n Expected number of cells.
#' @param xy_names Optional cell identifiers from `xy`.
#' @return Numeric matrix with cells in rows and uniquely named variables in
#'   columns.
#' @noRd
.prepareDiagnosticX <- function(X, n, xy_names) {
  if (is.numeric(X) && is.null(dim(X))) {
    x_rownames <- names(X)
    X <- matrix(X, ncol = 1L)
    rownames(X) <- x_rownames
    colnames(X) <- "X"
  } else if (is.matrix(X) && is.numeric(X)) {
    X <- as.matrix(X)
  } else {
    stop("X must be a numeric vector or matrix.", call. = FALSE)
  }

  if (nrow(X) != n) {
    stop("nrow(X) must equal nrow(xy).", call. = FALSE)
  }
  if (ncol(X) < 1L) {
    stop("X must contain at least one variable.", call. = FALSE)
  }
  if (!is.null(xy_names) && !is.null(rownames(X)) &&
      !identical(rownames(X), xy_names)) {
    stop("rownames(X) must match rownames(xy) in the same order.",
         call. = FALSE)
  }

  x_names <- colnames(X)
  if (is.null(x_names)) x_names <- rep("", ncol(X))
  missing_names <- is.na(x_names) | !nzchar(x_names)
  x_names[missing_names] <- paste0("X", which(missing_names))
  if (anyDuplicated(x_names)) {
    stop("X column names must be unique.", call. = FALSE)
  }
  colnames(X) <- x_names
  X
}


#' Validate a logged getPatches result
#'
#' @param result Complete output of `getPatches(log_iters = TRUE)`.
#' @param n Expected number of cells.
#' @param xy_names Optional cell identifiers from `xy`.
#' @return List containing validated `patch` and `membership_log` objects.
#' @noRd
.prepareDiagnosticResult <- function(result, n, xy_names) {
  if (!all(c("patch", "membership_log") %in% names(result))) {
    stop(
      "A list supplied as patches must be a complete getPatches(log_iters = TRUE) result containing patch and membership_log.",
      call. = FALSE
    )
  }

  .validateDiagnosticPatch(
    result$patch,
    n = n,
    xy_names = xy_names,
    label = "result$patch"
  )

  membership_log <- result$membership_log
  if (!is.matrix(membership_log) || nrow(membership_log) != n) {
    stop("result$membership_log must be a matrix with nrow(xy) rows.",
         call. = FALSE)
  }
  if (!is.null(xy_names) && !is.null(rownames(membership_log)) &&
      !identical(rownames(membership_log), xy_names)) {
    stop(
      "rownames(result$membership_log) must match rownames(xy) in the same order.",
      call. = FALSE
    )
  }
  if (ncol(membership_log) < 1L) {
    stop("result$membership_log must contain at least one iteration.",
         call. = FALSE)
  }
  if (!.samePatchAssignments(
    membership_log[, ncol(membership_log)],
    result$patch
  )) {
    stop(
      "The final column of result$membership_log must agree with result$patch.",
      call. = FALSE
    )
  }

  list(patch = result$patch, membership_log = membership_log)
}


#' Validate a patch vector
#'
#' @param patch Patch-assignment vector.
#' @param n Expected number of cells.
#' @param xy_names Optional cell identifiers from `xy`.
#' @param label Input label used in error messages.
#' @return `TRUE` invisibly when validation succeeds.
#' @noRd
.validateDiagnosticPatch <- function(patch, n, xy_names, label) {
  if (!is.atomic(patch) || !is.null(dim(patch))) {
    stop(sprintf("%s must be an atomic vector.", label), call. = FALSE)
  }
  if (length(patch) != n) {
    stop(sprintf("length(%s) must equal nrow(xy).", label), call. = FALSE)
  }
  if (!is.null(xy_names) && !is.null(names(patch)) &&
      !identical(names(patch), xy_names)) {
    stop(sprintf(
      "names(%s) must match rownames(xy) in the same order.",
      label
    ), call. = FALSE)
  }
  invisible(TRUE)
}


#' Compare patch assignments, treating paired missing values as equal
#'
#' @param x,y Patch-assignment vectors.
#' @return Logical scalar indicating whether the assignments agree.
#' @noRd
.samePatchAssignments <- function(x, y) {
  if (length(x) != length(y)) return(FALSE)
  x <- as.character(x)
  y <- as.character(y)
  same_missing <- is.na(x) == is.na(y)
  all(same_missing & (is.na(x) | x == y))
}


#' Calculate connectivity from k = 1 through max_k
#'
#' @param xy Numeric spatial coordinate matrix.
#' @param patch Character patch-assignment vector.
#' @param cell_indices Named list of cell indices for each ordered patch.
#' @param max_k Effective maximum number of neighbors.
#' @param strict_k Effective k for the strict connectivity summary.
#' @return List containing strict-k component fractions, minimum k for full
#'   connectivity, and the complete patch-by-k connectivity curve.
#' @noRd
.calculatePatchConnectivity <- function(xy, patch, cell_indices,
                                        max_k, strict_k) {
  patch_ids <- names(cell_indices)
  n_patches <- length(patch_ids)
  if (max_k == 0L) {
    return(list(
      strict_component_fraction = rep(1, n_patches),
      min_connectivity_k = rep(0L, n_patches),
      connectivity_curve = data.frame(
        patch = character(),
        k = integer(),
        component_fraction = numeric(),
        stringsAsFactors = FALSE
      )
    ))
  }

  nn_index <- FNN::get.knn(xy, k = max_k)$nn.index
  component_fraction <- matrix(
    NA_real_,
    nrow = n_patches,
    ncol = max_k,
    dimnames = list(patch_ids, paste0("k", seq_len(max_k)))
  )
  n <- nrow(xy)

  for (current_k in seq_len(max_k)) {
    neighbors <- nn_index[, seq_len(current_k), drop = FALSE]
    from <- rep(seq_len(n), each = current_k)
    to <- as.vector(t(neighbors))
    same_patch <- !is.na(patch[from]) & !is.na(patch[to]) &
      patch[from] == patch[to]

    graph <- igraph::make_empty_graph(n = n, directed = FALSE)
    if (any(same_patch)) {
      edges <- as.vector(t(cbind(from[same_patch], to[same_patch])))
      graph <- igraph::add_edges(graph, edges)
    }
    membership <- igraph::components(graph)$membership

    for (i in seq_along(patch_ids)) {
      component_sizes <- table(membership[cell_indices[[i]]])
      component_fraction[i, current_k] <-
        max(component_sizes) / length(cell_indices[[i]])
    }
  }

  minimum_k <- apply(component_fraction, 1L, function(values) {
    fully_connected <- which(values == 1)
    if (length(fully_connected) == 0L) NA_integer_ else fully_connected[1L]
  })
  minimum_k[lengths(cell_indices) == 1L] <- 0L

  connectivity_curve <- data.frame(
    patch = rep(patch_ids, times = max_k),
    k = rep(seq_len(max_k), each = n_patches),
    component_fraction = as.vector(component_fraction),
    stringsAsFactors = FALSE
  )

  list(
    strict_component_fraction = component_fraction[, strict_k],
    min_connectivity_k = as.integer(minimum_k),
    connectivity_curve = connectivity_curve
  )
}
